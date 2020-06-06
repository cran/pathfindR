## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      fig.width = 8, fig.height = 4, fig.align = "center")

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(pathfindR))

## ----enr_chart, eval=FALSE----------------------------------------------------
#  enrichment_chart(RA_output)

## ----enr_chart2, eval=FALSE---------------------------------------------------
#  ## change top_terms
#  enrichment_chart(RA_output, top_terms = 3)
#  
#  ## set null for displaying all terms
#  enrichment_chart(RA_output, top_terms = NULL)

## ----enr_chart3, fig.height=8, fig.width=8------------------------------------
enrichment_chart(RA_clustered, plot_by_cluster = TRUE)

## ----hmap---------------------------------------------------------------------
term_gene_heatmap(RA_output)

## ----hmap2, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(RA_output, num_terms = 3)
#  
#  ## set null for displaying all terms
#  term_gene_heatmap(RA_output, num_terms = NULL)

## ----hmap3, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(RA_output, use_description = TRUE)

## ----hmap4, eval=FALSE--------------------------------------------------------
#  term_gene_heatmap(result_df = RA_output, genes_df = RA_input)

## ----term_gene1---------------------------------------------------------------
term_gene_graph(RA_output)

## ----term_gene2, eval=FALSE---------------------------------------------------
#  term_gene_graph(RA_output, num_terms = NULL)

## ----term_gene3, eval=FALSE---------------------------------------------------
#  term_gene_graph(RA_output, num_terms = 3, use_description = TRUE)

## ----term_gene4, eval=FALSE---------------------------------------------------
#  term_gene_graph(RA_output, num_terms = 3, node_size = "p_val")

## ----upset1-------------------------------------------------------------------
UpSet_plot(RA_output)

## ----upset2, eval=FALSE-------------------------------------------------------
#  UpSet_plot(RA_output, genes_df = RA_input)

## ----upset3, eval=FALSE-------------------------------------------------------
#  UpSet_plot(RA_output, num_terms = 5)

## ----upset4, eval=FALSE-------------------------------------------------------
#  UpSet_plot(RA_output, use_description = TRUE)

## ----upset5, eval=FALSE-------------------------------------------------------
#  UpSet_plot(RA_output, method = "barplot")

## ----upset6, eval=FALSE-------------------------------------------------------
#  UpSet_plot(RA_output, RA_input, method = "boxplot")

