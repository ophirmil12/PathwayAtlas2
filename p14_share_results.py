import os
import pandas as pd
from definitions import RESULTS_P

def extract_novel_hits():
    ### TODO: mind the file name! ###
    input_csv = os.path.join(RESULTS_P, "gpt_biological_insights.csv")
    output_excel = os.path.join(RESULTS_P, "novel_pathway_hits.xlsx")

    print(f"Loading data from {input_csv}...")
    
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"Error: Could not find {input_csv}. Please check the path.")
        return

    # Filter for novel hits
    # Coercing to string and checking lower() catches cases where pandas reads it as text instead of bool
    novel_df = df[df['is_novel_hit'].astype(str).str.lower() == 'true']

    # Define the exact column order requested
    target_columns = [
        "cancer_type", 
        "pathway_id", 
        "pathway_name", 
        "q_value", 
        "delta_means", 
        "confidence_score", 
        "biological_relevance", 
        "suggested_mechanism"
    ]

    # Select and reorder (ignoring any columns not in this list)
    final_df = novel_df[target_columns]

    # Export to Excel
    print(f"Found {len(final_df)} novel hits. Saving to Excel...")
    final_df.to_excel(output_excel, index=False, engine='openpyxl')
    
    print(f"Done! File saved to: {output_excel}")

if __name__ == "__main__":
    extract_novel_hits()