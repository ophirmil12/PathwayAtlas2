#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PathwayAtlas2 Analysis Summarizer with GPT-5 (structured outputs + caching)

This script iterates through statistical result CSVs in the PathwayAtlas results folder,
summarizes significant pathways using GPT, and streams results to a consolidated CSV.

For each cancer type, it:
1. Takes the top 20 significant pathways (q_value <= 0.05) based on delta_means.
2. Calls GPT-5 to extract biological insights, focusing on novel findings and relevance.
3. Caches GPT responses to avoid redundant API calls.
4. Appends insights to a master CSV
"""

import os
import sys
import time
import json
import argparse
import glob
from pathlib import Path
from typing import Dict, Any, Optional, List

import pandas as pd
from openai import OpenAI
from pydantic import BaseModel

from definitions import RESULTS_DISTANCES_P, RESULTS_P, KEGG_PATHWAY_METADATA_FILE

# ------------------ Config ------------------

MODEL = os.getenv("OPENAI_MODEL", "gpt-5")
API_KEY = os.getenv("OPENAI_API_KEY")

SYSTEM_PROMPT = """You will be provided with statistical results (including Delta Means, Q-values, and gene coverage metrics) for biological pathways in a specific cancer cohort.
Your task is to extract, summarize, and contextualize the most statistically significant and biologically interesting pathways. While acknowledging known cancer hallmarks, your primary focus should be on finding novel findings.
Focus on explaining WHY a positive Delta Mean (pathogenic shift, suggesting tumor selection or pathway activation) or negative Delta Mean (benign shift, suggesting suppression or metabolic rerouting) might be relevant to that specific cancer's biology.
When formulating your insights, provide specific mechanistic hypotheses."""

# CSV schema for the summarized output
CSV_FIELDS = [
    "cancer_type", "pathway_id", "pathway_name", "q_value", "delta_means",
    "biological_relevance", "is_novel_hit", "suggested_mechanism",
    "confidence_score", "model", "cache_hit"
]


# ------------------ Structured Output Schema ------------------

class PathwayInsight(BaseModel):
    pathway_id: str
    pathway_name: str
    biological_relevance: str
    is_novel_hit: bool
    suggested_mechanism: str
    confidence_score: float


class AnalysisResults(BaseModel):
    insights: List[PathwayInsight]


# ------------------ CLI ------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Summarize PathwayAtlas Results with GPT-5")
    ap.add_argument("--q-threshold", type=float, default=0.05, help="Q-value threshold for inclusion")
    ap.add_argument("--out", default=os.path.join(RESULTS_P, "gpt_biological_insights.csv"), help="Output CSV path")
    ap.add_argument("--cache-dir", default=".cache_gpt_atlas", help="Cache directory")
    return ap.parse_args()

# ------------------ Cache & IO Helpers ------------------

def ensure_cache_dir(cache_dir: Path):
    cache_dir.mkdir(parents=True, exist_ok=True)


def load_cache(cache_dir: Path, cancer_type: str) -> Optional[dict]:
    p = cache_dir / f"{cancer_type}.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def save_cache(cache_dir: Path, cancer_type: str, data: dict):
    p = cache_dir / f"{cancer_type}.json"
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


# ------------------ GPT Interaction ------------------

def call_gpt_for_analysis(client: OpenAI, cancer_type: str, data_payload: str) -> dict:
    print(f"Requesting biological interpretation for {cancer_type}...")

    resp = client.beta.chat.completions.parse(
        model=MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"Cancer Cohort: {cancer_type}\nTop Results Data:\n{data_payload}"}
        ],
        response_format=AnalysisResults
    )

    return resp.choices[0].message.parsed.model_dump()


# ------------------ CSV Streaming ------------------

def append_and_sort_csv(path: Path, rows: List[Dict[str, Any]]):
    """Appends new rows, sorts the entire dataset by cancer_type, and overwrites the CSV."""
    if not rows: return
    
    new_df = pd.DataFrame(rows)
    new_df = new_df.reindex(columns=CSV_FIELDS)
    
    # Load existing data if the file exists and is not empty
    if path.exists() and path.stat().st_size > 0:
        try:
            existing_df = pd.read_csv(path)
            combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        except pd.errors.EmptyDataError:
            combined_df = new_df
    else:
        combined_df = new_df

    # Remove duplicated lines
    combined_df = combined_df.drop_duplicates(subset=["cancer_type", "pathway_id", "is_novel_hit"], keep="last")
        
    # Sort by cancer type
    combined_df = combined_df.sort_values(by="cancer_type")

    # move biological_relevance column to the end
    cols = combined_df.columns.tolist()
    cols.append(cols.pop(cols.index("biological_relevance")))
    combined_df = combined_df[cols]
    
    # Overwrite the CSV with the sorted dataframe
    combined_df.to_csv(path, index=False, encoding='utf-8')


# ------------------ Main Execution ------------------

def main():
    args = parse_args()

    if not API_KEY:
        print("OPENAI_API_KEY environment variable is required.")
        sys.exit(1)

    client = OpenAI(api_key=API_KEY)
    out_csv = Path(args.out)
    cache_dir = Path(args.cache_dir)
    ensure_cache_dir(cache_dir)

    # Find all result files in RESULTS_DISTANCES_P
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    print(f"Found {len(result_files)} result files to process.")

    for file_path in sorted(result_files):
        cancer_type = Path(file_path).stem
        df = pd.read_csv(file_path)

        # Take only human pathways (those starting with 'hsa')
        df = df[df['pathway'].str.startswith('hsa')]  # Keep only KEGG pathways (those starting with 'hsa')

        # Filter for significant results
        sig_df = df[df['q_value'] <= args.q_threshold].sort_values('delta_means', ascending=False).head(20)      ###### TODO: mind tail/head!

        if sig_df.empty:
            print(f"No significant results for {cancer_type}, skipping.")
            continue

        cached = load_cache(cache_dir, cancer_type)

        if cached:
            print(f"--- Cache hit for {cancer_type} [[]]")
            analysis = cached
            cache_hit = True
        else:
            # Prepare data summary for GPT
            data_summary = sig_df[['pathway', 'pathway_name', 'q_value', 'delta_means',
                                    'num_mutations', 'num_genes', 'num_covered_genes']].to_string(index=False)
            try:
                analysis = call_gpt_for_analysis(client, cancer_type, data_summary)
                save_cache(cache_dir, cancer_type, analysis)
                cache_hit = False
                time.sleep(2)               # Prevent rate limiting
            except Exception as e:
                print(f"Failed GPT analysis for {cancer_type}: {e}")
                continue

        # Prepare rows for streaming
        final_rows = []
        insights = analysis.get("insights", [])

        # Map statistical data back to GPT's insights
        for item in insights:
            stats = sig_df[sig_df['pathway'] == item['pathway_id']]
            q_val = stats['q_value'].values[0] if not stats.empty else None
            d_mean = stats['delta_means'].values[0] if not stats.empty else None

            final_rows.append({
                "cancer_type": cancer_type,
                "pathway_id": item['pathway_id'],
                "pathway_name": item['pathway_name'],
                "q_value": q_val,
                "delta_means": d_mean,
                "biological_relevance": item['biological_relevance'],
                "is_novel_hit": item['is_novel_hit'],
                "suggested_mechanism": item['suggested_mechanism'],
                "confidence_score": item['confidence_score'],
                "model": MODEL,
                "cache_hit": cache_hit
            })

        append_and_sort_csv(out_csv, final_rows)
        print(f"Successfully processed and appended insights for {cancer_type}.")

    print(f"Process complete. Biological insights saved to: {out_csv}")


if __name__ == "__main__":
    main()