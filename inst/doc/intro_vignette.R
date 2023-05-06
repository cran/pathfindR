## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  fig.width = 7, fig.height = 7, fig.align = "center"
)
suppressPackageStartupMessages(library(pathfindR))

## ----load_pkg, eval=TRUE------------------------------------------------------
library(pathfindR)
knitr::kable(head(example_pathfindR_input))

## ----run_pathfindR------------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input)

## ----change_input_thr---------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, p_val_threshold = 0.01)

## ----change_out_dir-----------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, output_dir = "this_is_my_output_directory")

## ----change_out_dir2----------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, output_dir = "~/Desktop/my_dir")

## ----change_gset1-------------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, gene_sets = "GO-MF")

## ----change_gset2-------------------------------------------------------------
#  ## Including more terms for enrichment analysis
#  output_df <- run_pathfindR(example_pathfindR_input,
#    gene_sets = "GO-MF",
#    min_gset_size = 5,
#    max_gset_size = 500
#  )

## ----change_enr_threshold-----------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input,
#    adj_method = "fdr",
#    enrichment_threshold = 0.01
#  )

## ----change_PIN1--------------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, pin_name_path = "IntAct")

## ----change_PIN2--------------------------------------------------------------
#  # to use an external PIN of your choice
#  output_df <- run_pathfindR(example_pathfindR_input, pin_name_path = "/path/to/myPIN.sif")

## ----change_method------------------------------------------------------------
#  # for simulated annealing:
#  output_df <- run_pathfindR(example_pathfindR_input, search_method = "SA")
#  # for genetic algorithm:
#  output_df <- run_pathfindR(example_pathfindR_input, search_method = "GA")

## ----change_n_iters-----------------------------------------------------------
#  output_df <- run_pathfindR(example_pathfindR_input, iterations = 25)

## ----change_n_proc------------------------------------------------------------
#  # if not set, `n_processes` defaults to (number of detected cores - 1)
#  output_df <- run_pathfindR(example_pathfindR_input, iterations = 5, n_processes = 2)

## ----example_out, eval=TRUE---------------------------------------------------
knitr::kable(head(example_pathfindR_output, 2))

## ----encrichment_plot_shown---------------------------------------------------
#  # change number of top terms plotted (default = 10)
#  enrichment_chart(
#    result_df = example_pathfindR_output,
#    top_terms = 15
#  )

## ----KEGG_vis, eval=FALSE-----------------------------------------------------
#  input_processed <- input_processing(example_pathfindR_input)
#  visualize_terms(
#    result_df = example_pathfindR_output,
#    input_processed = input_processed,
#    hsa_KEGG = TRUE
#  )

## ----nonKEGG_viss, eval=FALSE-------------------------------------------------
#  input_processed <- input_processing(example_pathfindR_input)
#  visualize_terms(
#    result_df = example_pathfindR_output,
#    input_processed = input_processed,
#    hsa_KEGG = FALSE,
#    pin_name_path = "Biogrid"
#  )

## ----hierarchical0------------------------------------------------------------
#  example_pathfindR_output_clustered <- cluster_enriched_terms(example_pathfindR_output, plot_dend = FALSE, plot_clusters_graph = FALSE)

## ----hierarchical1, eval=TRUE-------------------------------------------------
## First 2 rows of clustered data frame
knitr::kable(head(example_pathfindR_output_clustered, 2))
## The representative terms
knitr::kable(example_pathfindR_output_clustered[example_pathfindR_output_clustered$Status == "Representative", ])

## ----hierarchical2, eval=TRUE-------------------------------------------------
# plotting only selected clusters for better visualization
selected_clusters <- subset(example_pathfindR_output_clustered, Cluster %in% 5:7)
enrichment_chart(selected_clusters, plot_by_cluster = TRUE)

## ----fuzzy--------------------------------------------------------------------
#  clustered_fuzzy <- cluster_enriched_terms(example_pathfindR_output, method = "fuzzy")

## ----hmap1--------------------------------------------------------------------
#  term_gene_heatmap(result_df = example_pathfindR_output, genes_df = example_pathfindR_input)

## ----term_gene_graph----------------------------------------------------------
#  term_gene_graph(result_df = example_pathfindR_output, use_description = TRUE)

## ----upset--------------------------------------------------------------------
#  UpSet_plot(result_df = example_pathfindR_output, genes_df = example_pathfindR_input)

## ----scores, eval=TRUE, fig.height=4, fig.width=8-----------------------------
## Vector of "Case" IDs
cases <- c(
  "GSM389703", "GSM389704", "GSM389706", "GSM389708",
  "GSM389711", "GSM389714", "GSM389716", "GSM389717",
  "GSM389719", "GSM389721", "GSM389722", "GSM389724",
  "GSM389726", "GSM389727", "GSM389730", "GSM389731",
  "GSM389733", "GSM389735"
)

## Calculate scores for representative terms
## and plot heat map using term descriptions
representative_df <- example_pathfindR_output_clustered[example_pathfindR_output_clustered$Status == "Representative", ]
score_matrix <- score_terms(
  enrichment_table = representative_df,
  exp_mat = example_experiment_matrix,
  cases = cases,
  use_description = TRUE, # default FALSE
  label_samples = FALSE, # default = TRUE
  case_title = "RA", # default = "Case"
  control_title = "Healthy", # default = "Control"
  low = "#f7797d", # default = "green"
  mid = "#fffde4", # default = "black"
  high = "#1f4037" # default = "red"
)

## ----compare2res, eval=TRUE, fig.height=4, fig.width=8------------------------
combined_df <- combine_pathfindR_results(
  result_A = example_pathfindR_output,
  result_B = example_comparison_output,
  plot_common = FALSE
)

## ----custom_prep, eval=TRUE---------------------------------------------------
## CREB target genes
CREB_target_genes <- normalizePath(system.file("extdata/CREB.txt", package = "pathfindR"))
CREB_target_genes <- readLines(CREB_target_genes)[-c(1, 2)] # skip the first two lines

## MYC target genes
MYC_target_genes <- normalizePath(system.file("extdata/MYC.txt", package = "pathfindR"))
MYC_target_genes <- readLines(MYC_target_genes)[-c(1, 2)] # skip the first two lines

## Prep for use
custom_genes <- list(TF1 = CREB_target_genes, TF2 = MYC_target_genes)
custom_descriptions <- c(TF1 = "CREB target genes", TF2 = "MYC target genes")

## ----custom_input, eval=TRUE--------------------------------------------------
set.seed(123)

## Select 40 random genes from MYC gene sets and 10 from CREB gene sets
selected_genes <- sample(MYC_target_genes, 40)
selected_genes <- c(
  selected_genes,
  sample(CREB_target_genes, 10)
)

## Assign random p value between 0.001 and 0.05 for each selected gene
rand_p_vals <- sample(seq(0.001, 0.05, length.out = 5),
  size = length(selected_genes),
  replace = TRUE
)

example_pathfindR_input <- data.frame(
  Gene_symbol = selected_genes,
  p_val = rand_p_vals
)
knitr::kable(head(example_pathfindR_input))

## ----custom_run---------------------------------------------------------------
#  example_custom_genesets_result <- run_pathfindR(
#    example_pathfindR_input,
#    gene_sets = "Custom",
#    custom_genes = custom_genes,
#    custom_descriptions = custom_descriptions,
#    max_gset_size = Inf, # DO NOT LIMIT GENE SET SIZE
#    output_dir = "misc/CREB_MYC"
#  )
#  
#  knitr::kable(example_custom_genesets_result)

## ----custom_result1, eval=TRUE, echo=FALSE------------------------------------
knitr::kable(example_custom_genesets_result)

