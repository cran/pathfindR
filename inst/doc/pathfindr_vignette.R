## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(pathfindR))
data("RA_input")
knitr::kable(head(RA_input))

## ----eval=FALSE----------------------------------------------------------
#  RA_output <- run_pathfindR(RA_input)

## ----eval=FALSE----------------------------------------------------------
#  # to change the output directory
#  RA_output <- run_pathfindR(RA_input, output = "new_directory")
#  
#  # to change the PIN (default = Biogrid)
#  RA_output <- run_pathfindR(RA_input, pin_name = "IntAct")
#  # to use an external PIN of user's choice
#  RA_output <- run_pathfindR(RA_input, pin_name = "/path/to/myPIN.sif")
#  
#  # available gene sets are KEGG, Reactome, BioCarta, GO-BP, GO-CC and GO-MF
#  # default is KEGG
#  # to change the gene sets used for enrichment analysis
#  RA_output <- run_pathfindR(RA_input, gene_sets = "BioCarta")
#  
#  # to change the active subnetwork search algorithm (default = "GR", i.e. greedy algorithm)
#  # for simulated annealing:
#  RA_output <- run_pathfindR(RA_input, search_method = "SA")
#  
#  # to change the number of iterations (default = 10)
#  RA_output <- run_pathfindR(RA_input, iterations = 5)
#  
#  # to manually specify the number processes used during parallel loop by foreach
#  # defaults to the number of detected cores
#  RA_output <- run_pathfindR(RA_input, n_processes = 2)
#  
#  # to report the non-DEG active subnetwork genes
#  RA_output <- run_pathfindR(RA_input, list_active_snw_genes = TRUE)

## ------------------------------------------------------------------------
data("RA_output")
knitr::kable(head(RA_output, 2))

## ----hierarchical, fig.height=4, fig.width=8-----------------------------
data("RA_output")
RA_clustered <- cluster_pathways(RA_output)
## First 2 rows of clustered terms data frame
knitr::kable(head(RA_clustered, 2))
## The 8 representative terms
knitr::kable(RA_clustered[RA_clustered$Status == "Representative", ])

# to display the heatmap of kappa statistics
RA_clustered <- cluster_pathways(RA_output, plot_hmap = TRUE, plot_clusters_graph = FALSE)

# to display the dendrogram and optimal clusters
RA_clustered <- cluster_pathways(RA_output, plot_dend = TRUE, plot_clusters_graph = FALSE)

# to change agglomeration method (default = "average")
RA_clustered <- cluster_pathways(RA_output, hclu_method = "centroid")

## ----fuzzy, fig.height=4, fig.width=8------------------------------------
RA_clustered <- cluster_pathways(RA_output, method = "fuzzy")

## ---- fig.height=4, fig.width=8------------------------------------------
## Pathway data frame
pws_table <- pathfindR::RA_clustered
# selecting "Representative" pathways for clear visualization
pws_table <- pws_table[pws_table$Status == "Representative", ]

## Expression matrix
exp_mat <- pathfindR::RA_exp_mat

## Vector of "Case" IDs
cases <- c("GSM389703", "GSM389704", "GSM389706", "GSM389708", 
           "GSM389711", "GSM389714", "GSM389716", "GSM389717", 
           "GSM389719", "GSM389721", "GSM389722", "GSM389724", 
           "GSM389726", "GSM389727", "GSM389730", "GSM389731", 
           "GSM389733", "GSM389735")

## Calculate pathway scores and plot heatmap
score_matrix <- calculate_pw_scores(pws_table, exp_mat, cases)

## ------------------------------------------------------------------------
## CREB target genes
CREB_targets <- normalizePath(system.file("extdata/CREB.txt",
                                      package = "pathfindR"))
CREB_targets <- read.delim(CREB_targets, skip = 2, header = FALSE, stringsAsFactors = FALSE)
CREB_targets <- CREB_targets$V1

## MYC target genes
MYC_targets <- normalizePath(system.file("extdata/MYC.txt",
                                      package = "pathfindR"))
MYC_targets <- read.delim(MYC_targets, skip = 2, header = FALSE, stringsAsFactors = FALSE)
MYC_targets <- MYC_targets$V1

## Prep for use
custom_genes <- list(CREB_targets, MYC_targets)
names(custom_genes) <- c("TF1", "TF2")

custom_pathways <- c("CREB target genes", "MYC target genes")
names(custom_pathways) <- c("TF1", "TF2")

## ------------------------------------------------------------------------
set.seed(123)
selected_genes <- sample(MYC_targets, 50)
selected_genes <- c(selected_genes,
                    sample(CREB_targets, 5))
input <- data.frame(Gene.symbol = selected_genes,
                    logFC = 1.5,
                    adj.P.Val = 0.001)
head(input)

## ----eval=FALSE----------------------------------------------------------
#  custom_result <- run_pathfindR(input,
#                                 gene_sets = "Custom",
#                                 custom_genes = custom_genes,
#                                 custom_pathways = custom_pathways,
#                                 iterations = 1)

## ----include = FALSE-----------------------------------------------------
data(custom_result)

## ------------------------------------------------------------------------
knitr::kable(custom_result)

