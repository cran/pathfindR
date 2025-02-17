## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8, fig.height = 4, fig.align = "center"
)
suppressPackageStartupMessages(library(pathfindR))

## ----compare2res--------------------------------------------------------------
combined_df <- combine_pathfindR_results(
  result_A = example_pathfindR_output,
  result_B = example_comparison_output,
  plot_common = FALSE
)

## ----compare_graph1-----------------------------------------------------------
combined_results_graph(combined_df)

## ----compare_graph2, fig.width=8, fig.height=4--------------------------------
combined_results_graph(
  combined_df,
  selected_terms = c("hsa04144", "hsa04141", "hsa04140")
)

## ----compare_graph3, eval=FALSE-----------------------------------------------
# combined_results_graph(
#   combined_df,
#   use_description = TRUE,
#   selected_terms = combined_df$Term_Description[1:4]
# )

## ----compare_graph4, eval=FALSE-----------------------------------------------
# combined_results_graph(
#   combined_df,
#   selected_terms = c("hsa04144", "hsa04141", "hsa04140"),
#   node_size = "p_val"
# )

