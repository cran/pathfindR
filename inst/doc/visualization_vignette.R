## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8, fig.height = 4, fig.align = "center"
)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(pathfindR))

## ----enr_chart, eval=FALSE----------------------------------------------------
#  enrichment_chart(example_pathfindR_output)

## ----enr_chart2, eval=FALSE---------------------------------------------------
#  ## change top_terms
#  enrichment_chart(example_pathfindR_output, top_terms = 3)
#  
#  ## set null for displaying all terms
#  enrichment_chart(example_pathfindR_output, top_terms = NULL)

## ----enr_chart3, fig.height=8, fig.width=8------------------------------------
enrichment_chart(example_pathfindR_output_clustered, plot_by_cluster = TRUE)

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

## ----hmap---------------------------------------------------------------------
term_gene_heatmap(example_pathfindR_output)

## ----hmap2, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(example_pathfindR_output, num_terms = 3)
#  
#  ## set null for displaying all terms
#  term_gene_heatmap(example_pathfindR_output, num_terms = NULL)

## ----hmap3, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(example_pathfindR_output, use_description = TRUE)

## ----hmap4, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(result_df = example_pathfindR_output, genes_df = example_pathfindR_input)

## ----term_gene1---------------------------------------------------------------
term_gene_graph(example_pathfindR_output)

## ----term_gene2, eval=FALSE---------------------------------------------------
#  term_gene_graph(example_pathfindR_output, num_terms = NULL)

## ----term_gene3, eval=FALSE---------------------------------------------------
#  term_gene_graph(example_pathfindR_output, num_terms = 3, use_description = TRUE)

## ----term_gene4, eval=FALSE---------------------------------------------------
#  term_gene_graph(example_pathfindR_output, num_terms = 3, node_size = "p_val")

## ----upset1-------------------------------------------------------------------
UpSet_plot(example_pathfindR_output)

## ----upset2, eval=FALSE-------------------------------------------------------
#  UpSet_plot(example_pathfindR_output, genes_df = example_pathfindR_input)

## ----upset3, eval=FALSE-------------------------------------------------------
#  UpSet_plot(example_pathfindR_output, num_terms = 5)

## ----upset4, eval=FALSE-------------------------------------------------------
#  UpSet_plot(example_pathfindR_output, use_description = TRUE)

## ----upset5, eval=FALSE-------------------------------------------------------
#  UpSet_plot(example_pathfindR_output, method = "barplot")

## ----upset6, eval=FALSE-------------------------------------------------------
#  UpSet_plot(example_pathfindR_output, example_pathfindR_input, method = "boxplot")

