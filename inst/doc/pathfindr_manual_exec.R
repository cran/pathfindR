## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----init_steps, eval=FALSE----------------------------------------------
#  library(pathfindR)
#  data(RA_input)

## ----pin_path, eval=FALSE------------------------------------------------
#  pin_path <- return_pin_path(pin_name_path = "Biogrid")

## ----process, eval=FALSE-------------------------------------------------
#  RA_processed <- input_processing(input = RA_input, ## the input: in this case, differential expression results
#                                   p_val_threshold = 0.05, ## p value threshold to filter DEGs
#                                   pin_path = pin_path, ## path/to/Protein/Interaction/Network/File
#                                   human_genes = TRUE) ## boolean indicating whether the input genes are human gene symbols

## ----snw_search, eval=FALSE----------------------------------------------
#  n_iter <- 15 ## number of iterations
#  combined_res <- NULL ## to store each iteration's result
#  for (i in 1:n_iter) {
#    ###### Active Subnetwork Search
#    ## Name of output file
#    snws_file <- paste0("active_snws_", i, ".txt")
#    active_snws <- active_snw_search(RA_processed, pin_path, snws_file = snws_file, search_method = "GR")
#  
#    ###### Enrichment analyses
#    enrichment_res <- enrichment_analyses(active_snws, gene_sets = "GO-All",
#                                          pin_path = pin_path,
#                                          input_genes = RA_processed$GENE,
#                                          list_active_snw_genes = TRUE)
#    ## combine all results via rbind
#    combined_res <- rbind(combined_res, enrichment_res)
#  }

## ----post_proc, eval=FALSE-----------------------------------------------
#  ##### Summarize Combined Enrichment Results
#  final_res <- summarize_enrichment_results(combined_res, list_active_snw_genes = TRUE)
#  
#  ##### Annotate DEGs Involved in Each Pathway
#  final_res <- annotate_pathway_DEGs(final_res, input_processed = RA_processed, gene_sets = "GO-All")

## ----vis_pws, eval=FALSE-------------------------------------------------
#  visualize_pws(final_res, RA_processed, gene_sets = "GO-All", pin_name_path = "Biogrid")

## ----enr_chart, eval=FALSE-----------------------------------------------
#  png("enrichment_chart.png", width = 650, height = 800)
#  enrichment_chart(final_res)
#  dev.off()

