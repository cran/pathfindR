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
#  result <- run_pathfindR(RA_input)

## ----eval=FALSE----------------------------------------------------------
#  # to change the PIN (default = Biogrid)
#  result <- run_pathfindR(RA_input, pin_name = "IntAct")
#  # to use an external PIN of user's choice
#  result <- run_pathfindR(RA_input, pin_name = "/path/to/myPIN.sif")
#  
#  # available gene sets are KEGG, Reactome, BioCarta, GO-BP, GO-CC and GO-MF
#  # default is KEGG
#  # to change the gene sets used for enrichment analysis
#  result <- run_pathfindR(RA_input, gene_sets = "BioCarta")
#  
#  # to change the active subnetwork search algorithm (default = "GR", i.e. greedy algorithm)
#  # for simulated annealing:
#  result <- run_pathfindR(RA_input, search_method = "SA")
#  
#  # to change the number of iterations (default = 10)
#  result <- run_pathfindR(RA_input, iterations = 5)
#  
#  # to manually specify the number processes used during parallel loop by foreach
#  # defaults to the number of detected cores
#  result <- run_pathfindR(RA_input, n_processes = 2)

## ------------------------------------------------------------------------
data("RA_output")
knitr::kable(head(RA_output, 2))

## ----eval=FALSE----------------------------------------------------------
#  # example run
#  choose_clusters(RA_output)
#  
#  # to display the heatmap of pathway clustering
#  choose_clusters(RA_output, plot_heatmap = TRUE)
#  
#  # to display the heatmap of pathway clustering
#  # and change agglomeration method (default = "average")
#  choose_clusters(RA_output, agg_method = "complete", plot_heatmap = TRUE)
#  

