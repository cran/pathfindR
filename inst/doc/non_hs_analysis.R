## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      fig.width = 7, fig.height = 7, fig.align = "center",
                      eval = FALSE)
suppressPackageStartupMessages(library(pathfindR))

## ----mmu_kegg-----------------------------------------------------------------
#  gsets_list <- get_gene_sets_list(source = "KEGG",
#                                   org_code = "mmu")

## ----KEGG_save----------------------------------------------------------------
#  mmu_kegg_genes <- gsets_list$gene_sets
#  mmu_kegg_descriptions <- gsets_list$descriptions
#  
#  ## Save both as RDS files for later use
#  saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
#  saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")

## ----KEGG_load----------------------------------------------------------------
#  mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
#  mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")

## ----process_PIN1-------------------------------------------------------------
#  ## Downloading the STRING PIN file to tempdir
#  url <- "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz"
#  path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
#  download.file(url, path2file)
#  
#  ## read STRING pin file
#  mmu_string_df <- read.table(path2file, header = TRUE)
#  
#  ## filter using combined_score cut-off value of 800
#  mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]
#  
#  ## fix ids
#  mmu_string_pin <- data.frame(Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
#                               Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2))
#  head(mmu_string_pin, 2)

## ----process_PIN2-------------------------------------------------------------
#  library(biomaRt)
#  
#  mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#  
#  converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
#                     filters = "ensembl_peptide_id",
#                     values = unique(unlist(mmu_string_pin)),
#                     mart = mmu_ensembl)
#  mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
#  mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
#  mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
#  mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]
#  
#  head(mmu_string_pin, 2)

## ----process_PIN3-------------------------------------------------------------
#  # remove self interactions
#  self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
#  mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]
#  
#  # remove duplicated inteactions (including symmetric ones)
#  mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) # this will return a matrix object
#  
#  mmu_string_pin <- data.frame(A = mmu_string_pin[, 1],
#                               pp = "pp",
#                               B = mmu_string_pin[, 2])

## ----process_PIN4-------------------------------------------------------------
#  path2SIF <- file.path(tempdir(), "mmusculusPIN.sif")
#  write.table(mmu_string_pin,
#              file = path2SIF,
#              col.names = FALSE,
#              row.names = FALSE,
#              sep = "\t",
#              quote = FALSE)
#  path2SIF <- normalizePath(path2SIF)

## ----mmu_input_df, eval=TRUE--------------------------------------------------
knitr::kable(head(myeloma_input))

## ----run----------------------------------------------------------------------
#  myeloma_output <- run_pathfindR(input = myeloma_input,
#                                  convert2alias = FALSE,
#                                  gene_sets = "Custom",
#                                  custom_genes = mmu_kegg_genes,
#                                  custom_descriptions = mmu_kegg_descriptions,
#                                  pin_name_path = path2SIF)

## ----enr_chart, echo=FALSE, eval=TRUE-----------------------------------------
enrichment_chart(myeloma_output)

## ----output, eval=TRUE--------------------------------------------------------
knitr::kable(myeloma_output)

## ----run2---------------------------------------------------------------------
#  myeloma_output <- run_pathfindR(input = myeloma_input,
#                                  convert2alias = FALSE,
#                                  gene_sets = "mmu_KEGG",
#                                  pin_name_path = "mmu_STRING")

