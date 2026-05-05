"""
Mappings between Bailey et al. (2018) Table S2 and the PathwayAtlas project.

Sources:
  - Table S2: "Biological processes and pathways associated with each driver gene"
    (Bailey et al. 2018, Cell, PMC6029450)
  - Cancer files: results/distances/<name>.csv
  - Pathway IDs: data/kegg/pathway_id_to_metadata.pickle (KEGG hsa IDs)

Notes on unmapped entries are included as comments below each dict.
"""

# ──────────────────────────────────────────────────────────────────────────────
# 1.  S2 cancer type  →  cancer file name (without .csv)
# ──────────────────────────────────────────────────────────────────────────────
#
# S2 uses standard TCGA codes; files use lowercased, sometimes
# consolidated names.  Matches were determined by cross-referencing
# TCGA cohort definitions with the file list.


s2_cancer_to_mine = {
    # Direct / trivially lowercased matches
    "ACC":     "acc",        # Adrenocortical carcinoma
    "BLCA":    "blca",       # Bladder urothelial carcinoma
    "BRCA":    "brca",       # Breast invasive carcinoma
    "CESC":    "cesc",       # Cervical squamous cell carcinoma & endocervical adenocarcinoma
    "CHOL":    "chol",       # Cholangiocarcinoma
    "COADREAD":"coadread",   # Colorectal (COAD + READ combined)
    "DLBC":    "dlbclnos",   # Diffuse large B-cell lymphoma  (file adds "nos")
    "ESCA":    "esca",       # Esophageal carcinoma
    "HNSC":    "hnsc",       # Head & neck squamous cell carcinoma
    "LUAD":    "luad",       # Lung adenocarcinoma
    "LUSC":    "lusc",       # Lung squamous cell carcinoma
    "PAAD":    "paad",       # Pancreatic adenocarcinoma
    "PRAD":    "prad",       # Prostate adenocarcinoma
    "SKCM":    "skcm",       # Skin cutaneous melanoma
    "STAD":    "stad",       # Stomach adenocarcinoma
    "THYM":    "thym",       # Thymoma
    "UCEC":    "ucec",       # Uterine corpus endometrial carcinoma
    "UCS":     "ucs",        # Uterine carcinosarcoma
    "PANCAN":  "pan_cancer", # Pan-cancer (all cohorts combined)

    # Renamed / consolidated in project
    "LAML":    "aml",        # Acute myeloid leukemia → aml
    "LIHC":    "hcc",        # Liver hepatocellular carcinoma → hcc
    "LGG":     "difg",       # Lower-grade glioma → difg (diffuse glioma)
    "GBM":     "glioma",     # Glioblastoma multiforme → glioma (broader; GBM-only not in set)
    "KIRC":    "ccrcc",      # Kidney renal clear cell carcinoma → ccrcc
    "KIRP":    "prcc",       # Kidney renal papillary cell carcinoma → prcc
    "KICH":    "chrcc",      # Kidney chromophobe → chrcc
    "MESO":    "plmeso",     # Mesothelioma → plmeso (pleural mesothelioma)
    "OV":      "hgsoc",      # Ovarian serous cystadenocarcinoma → hgsoc
    "TGCT":    "testicular", # Testicular germ cell tumor → testicular (also 'testis'/'nsgct')
    "THCA":    "thpa",       # Thyroid carcinoma → thpa (thyroid papillary)
    "UVM":     "um",         # Uveal melanoma → um

    # Approximate / closest match (no exact equivalent)
    "SARC":    "soft_tissue", # Sarcoma → soft_tissue (closest available)

    # No match in file list
    "PCPG":    None,          # Pheochromocytoma & paraganglioma – not in results
}


# ──────────────────────────────────────────────────────────────────────────────
# 2.  S2 pathway name  →  KEGG hsa ID(s) from pathway_id_to_metadata.pickle
# ──────────────────────────────────────────────────────────────────────────────
#
# The S2 pathways are broad biological process categories, not single KEGG
# pathways.  Each is therefore mapped to the most representative KEGG pathway(s)
# that capture the same biology.  Where the S2 category is genuinely a broad
# umbrella (e.g. "Metabolism", "Chromatin other") several KEGG IDs are given.
#
# Categories that have NO clean KEGG equivalent in pickle
# ("Other", "Other signaling", "Transcription factor") are mapped to None
# with a comment.

s2_pathway_to_hsa = {

    "Apoptosis": [
        "hsa04210",   # Apoptosis
        "hsa04215",   # Apoptosis - multiple species
        "hsa04115",   # p53 signaling pathway  (major apoptosis regulator)
    ],

    "Cell cycle": [
        "hsa04110",   # Cell cycle
        "hsa04115",   # p53 signaling pathway  (G1/S & G2/M checkpoint)
    ],

    "Chromatin SWI/SNF complex": [
        "hsa03082",   # ATP-dependent chromatin remodeling
                      # (SWI/SNF is the main ATP-dependent remodeler in KEGG)
    ],

    "Chromatin histone modifiers": [
        "hsa03082",   # ATP-dependent chromatin remodeling
                      # Note: KEGG lacks a dedicated "histone modifier" pathway;
                      # hsa03082 is the closest available entry.
    ],

    "Chromatin other": [
        "hsa03082",   # ATP-dependent chromatin remodeling
                      # Catch-all for chromatin regulators not in SWI/SNF or
                      # histone-modifier subcategories.
    ],

    "Epigenetics DNA modifiers": [
        "hsa03082",   # ATP-dependent chromatin remodeling
                      # KEGG does not have a standalone "DNA methylation" pathway;
                      # DNA-modifier enzymes (DNMT, TET) appear in broader entries.
    ],

    "Genome integrity": [
        "hsa03410",   # Base excision repair
        "hsa03420",   # Nucleotide excision repair
        "hsa03430",   # Mismatch repair
        "hsa03440",   # Homologous recombination
        "hsa03450",   # Non-homologous end-joining
        "hsa03460",   # Fanconi anemia pathway
        "hsa04115",   # p53 signaling pathway  (DNA-damage checkpoint)
    ],

    "Histone modification": [
        "hsa03082",   # ATP-dependent chromatin remodeling
                      # (same caveat as "Chromatin histone modifiers" above)
    ],

    "Immune signaling": [
        "hsa04060",   # Cytokine-cytokine receptor interaction
        "hsa04062",   # Chemokine signaling pathway
        "hsa04620",   # Toll-like receptor signaling pathway
        "hsa04630",   # JAK-STAT signaling pathway
        "hsa04660",   # T cell receptor signaling pathway
        "hsa04662",   # B cell receptor signaling pathway
    ],

    "MAPK signaling": [
        "hsa04010",   # MAPK signaling pathway
        "hsa04014",   # Ras signaling pathway  (upstream activator of MAPK)
    ],

    "Metabolism": [
        "hsa01100",   # Metabolic pathways  (master KEGG map)
        "hsa05230",   # Central carbon metabolism in cancer
    ],

    "NFKB signaling": [
        "hsa04064",   # NF-kappa B signaling pathway
    ],

    "NOTCH signaling": [
        "hsa04330",   # Notch signaling pathway
    ],

    "Other": None,
    # "Other" in S2 = genes with driver status but no clear pathway assignment.
    # No meaningful KEGG mapping exists.

    "Other signaling": None,
    # "Other signaling" = signaling genes that don't fit the named categories.
    # No single KEGG pathway is appropriate.

    "PI3K signaling": [
        "hsa04151",   # PI3K-Akt signaling pathway
        "hsa04152",   # AMPK signaling pathway  (mTORC1 branch; shares PI3K nodes)
        "hsa04068",   # FoxO signaling pathway   (downstream of PI3K-Akt)
    ],

    "Protein homeostasis/ubiquitination": [
        "hsa04120",   # Ubiquitin mediated proteolysis
        "hsa03050",   # Proteasome
    ],

    "RNA abundance": [
        "hsa03020",   # RNA polymerase        (transcription)
        "hsa03015",   # mRNA surveillance pathway
        "hsa03018",   # RNA degradation
    ],

    "RTK signaling": [
        "hsa04012",   # ErbB signaling pathway      (EGFR/HER2 family)
        "hsa04014",   # Ras signaling pathway        (shared RTK downstream)
        "hsa04370",   # VEGF signaling pathway       (VEGFR)
        "hsa01521",   # EGFR tyrosine kinase inhibitor resistance
    ],

    "Splicing": [
        "hsa03040",   # Spliceosome
    ],

    "TGFB signaling": [
        "hsa04350",   # TGF-beta signaling pathway
    ],

    "TOR signaling": [
        "hsa04150",   # mTOR signaling pathway
    ],

    "Transcription factor": None,
    # "Transcription factor" in S2 = genes encoding TFs that act as oncogenes
    # or tumor suppressors (e.g. TP53, MYC).  KEGG has no single "transcription
    # factor" pathway; these genes appear across many KEGG pathways.

    "Wnt/B-catenin signaling": [
        "hsa04310",   # Wnt signaling pathway
    ],
}
