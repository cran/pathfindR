## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      eval = FALSE,
                      fig.width = 7, fig.height = 7, fig.align = "center")
suppressPackageStartupMessages(library(pathfindR))

## ----load_pkg, eval=TRUE-------------------------------------------------
library(pathfindR)
knitr::kable(head(RA_input))

## ----run_pathfindR-------------------------------------------------------
#  output_df <- run_pathfindR(input_df)

## ----change_input_thr----------------------------------------------------
#  output_df <- run_pathfindR(input_df, p_val_threshold = 0.01)

## ----change_out_dir------------------------------------------------------
#  output_df <- run_pathfindR(input_df, output_dir = "this_is_my_output_directory")

## ----change_out_dir2-----------------------------------------------------
#  output_df <- run_pathfindR(input_df, output_dir = "~/Desktop/my_dir")

## ----change_gset1--------------------------------------------------------
#  output_df <- run_pathfindR(input_df, gene_sets = "GO-MF")

## ----change_gset2--------------------------------------------------------
#  ## Including more terms for enrichment analysis
#  output_df <- run_pathfindR(input_df,
#                             gene_sets = "GO-MF",
#                             min_gset_size = 5,
#                             max_gset_size = 500)

## ----change_enr_threshold------------------------------------------------
#  output_df <- run_pathfindR(input_df,
#                             adj_method = "fdr",
#                             enrichment_threshold = 0.01)

## ----change_PIN1---------------------------------------------------------
#  output_df <- run_pathfindR(input_df, pin_name_path = "IntAct")

## ----change_PIN2---------------------------------------------------------
#  # to use an external PIN of your choice
#  output_df <- run_pathfindR(input_df, pin_name_path = "/path/to/myPIN.sif")

## ----change_method-------------------------------------------------------
#  # for simulated annealing:
#  output_df <- run_pathfindR(input_df, search_method = "SA")
#  # for genetic algorithm:
#  output_df <- run_pathfindR(input_df, search_method = "GA")

## ----change_n_iters------------------------------------------------------
#  output_df <- run_pathfindR(input_df, iterations = 25)

## ----change_n_proc-------------------------------------------------------
#  # n_processes defaults to (number of detected cores - 1)
#  output_df <- run_pathfindR(input_df, n_processes = 2)

## ----example_out, eval=TRUE----------------------------------------------
knitr::kable(head(RA_output, 2))

## ----encrichment_plot, eval=TRUE, echo=FALSE-----------------------------
enrichment_chart(RA_output)

## ----encrichment_plot_shown, eval=FALSE----------------------------------
#  # change number of top terms plotted (default = 10)
#  enrichment_chart(result_df = RA_output,
#                   top_terms = 15)

## ----hierarchical1, eval=TRUE--------------------------------------------
RA_clustered <- cluster_enriched_terms(RA_output, plot_dend = FALSE)
## First 2 rows of clustered data frame
knitr::kable(head(RA_clustered, 2))
## The representative terms
knitr::kable(RA_clustered[RA_clustered$Status == "Representative", ])

## ----hierarchical2, eval=TRUE--------------------------------------------
# plotting only selected clusters for better visualization
RA_selected <- subset(RA_clustered, Cluster %in% 5:7)
enrichment_chart(RA_selected, plot_by_cluster = TRUE)

## ----fuzzy, eval=TRUE, fig.height=4, fig.width=8-------------------------
RA_clustered_fuzzy <- cluster_enriched_terms(RA_output, method = "fuzzy")

## ----term_gene1, eval=TRUE, fig.height=4, fig.width=8--------------------
term_gene_graph(RA_output)

## ----term_gene2, eval=TRUE, fig.height=4, fig.width=8--------------------
term_gene_graph(RA_output, num_terms = NULL)

## ----term_gene3, eval=TRUE, fig.height=4, fig.width=8--------------------
term_gene_graph(RA_output, num_terms = 3, use_description = TRUE)

## ----term_gene4, eval=TRUE, fig.height=4, fig.width=8--------------------
term_gene_graph(RA_output, num_terms = 3, node_size = "p_val")

## ----scores, eval=TRUE, fig.height=4, fig.width=8------------------------
## Vector of "Case" IDs
cases <- c("GSM389703", "GSM389704", "GSM389706", "GSM389708", 
           "GSM389711", "GSM389714", "GSM389716", "GSM389717", 
           "GSM389719", "GSM389721", "GSM389722", "GSM389724", 
           "GSM389726", "GSM389727", "GSM389730", "GSM389731", 
           "GSM389733", "GSM389735")

## Calculate scores for representative terms 
## and plot heat map using term descriptions
score_matrix <- score_terms(enrichment_table = RA_clustered[RA_clustered$Status == "Representative", ],
                            exp_mat = RA_exp_mat,
                            cases = cases,
                            use_description = TRUE, # default FALSE
                            label_samples = FALSE, # default = TRUE
                            case_title = "RA",  # default = "Case"
                            control_title = "Healthy", # default = "Control"
                            low = "#f7797d", # default = "green"
                            mid = "#fffde4", # default = "black"
                            high = "#1f4037") # default = "red"

## ----custom_prep, eval=TRUE----------------------------------------------
## CREB target genes
CREB_target_genes <- normalizePath(system.file("extdata/CREB.txt", package = "pathfindR"))
CREB_target_genes <- readLines(CREB_target_genes)[-c(1, 2)] # skip the first two lines

## MYC target genes
MYC_target_genes <- normalizePath(system.file("extdata/MYC.txt", package = "pathfindR"))
MYC_target_genes <- readLines(MYC_target_genes)[-c(1, 2)] # skip the first two lines

## Prep for use
custom_genes <- list(TF1 = CREB_target_genes, TF2 = MYC_target_genes)
custom_descriptions <- c(TF1 = "CREB target genes", TF2 = "MYC target genes")

## ----custom_input, eval=TRUE---------------------------------------------
set.seed(123)

## Select 40 random genes from MYC gene sets and 10 from CREB gene sets
selected_genes <- sample(MYC_target_genes, 40)
selected_genes <- c(selected_genes, 
                    sample(CREB_target_genes, 10))

## Assign random p value between 0.001 and 0.05 for each selected gene
rand_p_vals <- sample(seq(0.001, 0.05, length.out = 5),
                      size = length(selected_genes),
                      replace = TRUE)

input_df <- data.frame(Gene_symbol = selected_genes,
                       p_val = rand_p_vals)
knitr::kable(head(input_df))

## ----custom_run----------------------------------------------------------
#  custom_result <- run_pathfindR(input_df,
#                                 gene_sets = "Custom",
#                                 custom_genes = custom_genes,
#                                 custom_descriptions = custom_descriptions,
#                                 max_gset_size = Inf, # DO NOT LIMIT GENE SET SIZE
#                                 iterations = 1, # for demo, setting number of iterations to 1
#                                 output_dir = "misc/v1.4/CREB_MYC")
#  
#  knitr::kable(custom_result)

## ----custom_result1, eval=TRUE-------------------------------------------
knitr::kable(custom_result)

