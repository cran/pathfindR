## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----init_steps, eval=TRUE----------------------------------------------------
suppressPackageStartupMessages(library(pathfindR))
data(example_pathfindR_input)
head(example_pathfindR_input, 3)

## ----process------------------------------------------------------------------
#  example_processed <- input_processing(
#    input = example_pathfindR_input, # the input: in this case, differential expression results
#    p_val_threshold = 0.05, # p value threshold to filter significant genes
#    pin_name_path = "Biogrid", # the name of the PIN to use for active subnetwork search
#    convert2alias = TRUE # boolean indicating whether or not to convert missing symbols to alias symbols in the PIN
#  )

## ----gene_set-----------------------------------------------------------------
#  # using "BioCarta" as our gene sets for enrichment
#  biocarta_list <- fetch_gene_set(
#    gene_sets = "BioCarta",
#    min_gset_size = 10,
#    max_gset_size = 300
#  )
#  biocarta_gsets <- biocarta_list[[1]]
#  biocarta_descriptions <- biocarta_list[[2]]

## ----snw_search---------------------------------------------------------------
#  n_iter <- 10 ## number of iterations
#  combined_res <- NULL ## to store the result of each iteration
#  
#  for (i in 1:n_iter) {
#    ###### Active Subnetwork Search
#    snws_file <- paste0("active_snws_", i) # Name of output file
#    active_snws <- active_snw_search(
#      input_for_search = example_processed,
#      pin_name_path = "Biogrid",
#      snws_file = snws_file,
#      score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
#      sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
#      search_method = "GR", # we suggest using GR
#      seedForRandom = i # setting seed to ensure reproducibility per iteration
#    )
#  
#    ###### Enrichment Analyses
#    current_res <- enrichment_analyses(
#      snws = active_snws,
#      sig_genes_vec = example_processed$GENE,
#      pin_name_path = "Biogrid",
#      genes_by_term = biocarta_gsets,
#      term_descriptions = biocarta_descriptions,
#      adj_method = "bonferroni",
#      enrichment_threshold = 0.05,
#      list_active_snw_genes = TRUE
#    ) # listing the non-input active snw genes in output
#  
#    ###### Combine results via `rbind`
#    combined_res <- rbind(combined_res, current_res)
#  }

## ----post_proc----------------------------------------------------------------
#  ###### Summarize Combined Enrichment Results
#  summarized_df <- summarize_enrichment_results(combined_res,
#    list_active_snw_genes = TRUE
#  )
#  
#  ###### Annotate Affected Genes Involved in Each Enriched Term
#  final_res <- annotate_term_genes(
#    result_df = summarized_df,
#    input_processed = example_processed,
#    genes_by_term = biocarta_gsets
#  )

## ----vis_pws------------------------------------------------------------------
#  visualize_terms(
#    result_df = final_res,
#    hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
#    pin_name_path = "Biogrid"
#  )

## ----enr_chart----------------------------------------------------------------
#  enrichment_chart(final_res[1:10, ])

