# ==============================================================================
# KEGG Pathway Clustering & Redundancy Reduction Pipeline
#
# EXECUTION INSTRUCTIONS:
# 1. Activate your environment:
#    conda activate project_env
#
# 2. Install required packages (run once):
#    mamba install -c conda-forge r-dplyr r-tidyr r-readr r-stringr -c bioconda bioconductor-simplifyenrichment bioconductor-keggrest
#
# 3. Execute the script:
#    /cs/labs/dina/ophirmil12/miniforge3/envs/project_env/bin/Rscript p12_pathways_clustering.R
# ==============================================================================


# --- 1. Library Loading & Paths ---
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(simplifyEnrichment)
library(KEGGREST)

RESULTS_DISTANCES_P <- "/cs/labs/dina/ophirmil12/PathwayAtlas2/results/distances"
RESULTS_P <- "/cs/labs/dina/ophirmil12/PathwayAtlas2/results"

# --- 2. The Dynamic Function ---
run_clustering <- function(db_type) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("STARTING PIPELINE FOR: %s\n", toupper(db_type)))
  cat(sprintf("========================================\n"))

  # 1. Load Data
  files <- list.files(path = RESULTS_DISTANCES_P, pattern = "\\.csv$", full.names = TRUE)
  data_list <- lapply(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE) %>%
      filter(q_value < 0.05) %>%
      select(pathway, delta_means) %>%
      mutate(cancer = str_replace(basename(f), "\\.csv$", ""))
    return(df)
  })
  
  # TODO: it might be worth changing median() to max():  ~ max(.x, na.rm = TRUE) in here just once, regenerating the heatmaps, and seeing if those "pathogenic clusters" visually pop out much clearer than they do with the median
  df_main <- bind_rows(data_list)
  matrix_df <- df_main %>%
    pivot_wider(names_from = cancer, values_from = delta_means, values_fill = list(delta_means = 0), values_fn = ~ median(.x, na.rm = TRUE))
  
  pathway_ids <- matrix_df$pathway

  # 2. KEGG API Dynamic Fetch
  cat(sprintf("Fetching %s data from API...\n", db_type))
  kegg_links <- keggLink(db_type, "hsa")
  
  # Strip prefixes like "path:" or "md:"
  raw_target_ids <- str_replace(unname(kegg_links), "^[^:]+:", "")
  
  # If it's a module, strip the "hsa_" prefix so it matches your CSVs exactly
  if(db_type == "module") {
      raw_target_ids <- str_replace(raw_target_ids, "^hsa_", "")
  }
  
  gene_to_pathway <- data.frame(
      gene = names(kegg_links),
      target_id = raw_target_ids
  )
  
  kegg_list <- split(gene_to_pathway$gene, gene_to_pathway$target_id)
  kegg_list_filtered <- kegg_list[intersect(names(kegg_list), pathway_ids)]
  
  if(length(kegg_list_filtered) == 0) {
      cat(sprintf("No matching %ss found in your significant data. Skipping.\n", db_type))
      return()
  }

  # 3. Native Jaccard
  calc_jaccard <- function(list_of_sets) {
    n <- length(list_of_sets)
    mat <- matrix(1, n, n, dimnames = list(names(list_of_sets), names(list_of_sets)))
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        sim <- length(intersect(list_of_sets[[i]], list_of_sets[[j]])) / length(unique(c(list_of_sets[[i]], list_of_sets[[j]])))
        mat[i, j] <- mat[j, i] <- sim
      }
    }
    return(mat)
  }
  sim_mat <- calc_jaccard(kegg_list_filtered)
  original_ids <- rownames(sim_mat)

  # 4. Fetch Names
  # ==========================================
  cat(sprintf("Fetching %s names...\n", db_type))
  
  # KEGG API throws a 400 Bad Request if you ask for organism-specific modules.
  # We must fetch the global module list instead.
  if(db_type == "module") {
      kegg_names_raw <- keggList("module")
  } else {
      kegg_names_raw <- keggList(db_type, "hsa")
  }
  
  clean_names <- str_replace(kegg_names_raw, " - Homo sapiens \\(human\\)", "")
  
  # Strip prefixes from the names list to match our matrix
  names_keys <- str_replace(names(clean_names), "^[^:]+:", "")
  if(db_type == "module") {
      names_keys <- str_replace(names_keys, "^hsa_", "") # Fallback just in case
  }
  names(clean_names) <- names_keys
  
  target_descriptions <- clean_names[original_ids]
  target_descriptions[is.na(target_descriptions)] <- original_ids[is.na(target_descriptions)]
  rownames(sim_mat) <- colnames(sim_mat) <- target_descriptions

  # 5. Cluster & Save (Dynamically named files)
  cl <- binary_cut(sim_mat, cutoff = 0.75)          # here chose the ~# of clusters
  df_clusters <- data.frame(id = original_ids, name = target_descriptions, cluster = cl, stringsAsFactors = FALSE)
  
  write_csv(df_clusters, file.path(RESULTS_P, sprintf("p12_%s_clusters_summary.csv", db_type)))
  
  cluster_annotations <- df_clusters %>%
    group_by(cluster) %>%
    summarise(count = n(), contained = paste(name, collapse = " | ")) %>%
    arrange(desc(count))
  write_csv(cluster_annotations, file.path(RESULTS_P, sprintf("p12_%s_cluster_annotations.csv", db_type)))
  
  pdf(file.path(RESULTS_P, sprintf("p12_simplified_kegg_%s_heatmap.pdf", db_type)), width = 15, height = 12)
  ht_clusters(sim_mat, cl, column_title = sprintf("KEGG Functional Clusters (%s)", toupper(db_type)))
  dev.off()
  
  cat(sprintf("Finished %s pipeline.\n", db_type))
}

# --- 3. Execute for Both ---
run_clustering("pathway")
run_clustering("module")