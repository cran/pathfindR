## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      fig.width = 7, fig.height = 7, fig.align = "center")
suppressPackageStartupMessages(library(pathfindR))

## ----KEGG, eval=FALSE---------------------------------------------------------
#  library(KEGGREST)
#  #### Obtain list of M.musculus pathways
#  mmu_kegg_descriptions <- keggList("pathway", "mmu")
#  
#  # Shorten descriptions
#  mmu_kegg_descriptions <- sub(" - Mus musculus \\(mouse\\)", "", mmu_kegg_descriptions)
#  
#  # Turn the identifiers into KEGG-style pathway identifiers (org_id#####)
#  names(mmu_kegg_descriptions) <- sub("path:", "", names(mmu_kegg_descriptions))
#  
#  #### Obtain and parse genes per each pathway
#  mmu_kegg_genes <- sapply(names(mmu_kegg_descriptions), function(pwid){
#    pw <- keggGet(pwid)
#    pw <- pw[[1]]$GENE[c(FALSE, TRUE)] # get gene symbols, not descriptions
#    pw <- sub(";.+", "", pw) # discard any remaining description
#    pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] # remove any mistaken lines that cannot be gene symbols
#    pw <- unique(pw) # keep unique symbols
#    return(pw)
#  })
#  
#  #### Filter terms to exclude those with 0 genes (metabolic pathways)
#  mmu_kegg_genes <- mmu_kegg_genes[sapply(mmu_kegg_genes, length) != 0]
#  mmu_kegg_descriptions <- mmu_kegg_descriptions[names(mmu_kegg_descriptions) %in% names(mmu_kegg_genes)]

## ----KEGG_save, eval=FALSE----------------------------------------------------
#  ## Save both as RDS files for later use
#  saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
#  saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")

## ----KEGG_load, eval=FALSE----------------------------------------------------
#  mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")

## ----process_PIN1, eval=FALSE-------------------------------------------------
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

## ----process_PIN2, eval=FALSE-------------------------------------------------
#  library(biomaRt)
#  mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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

## ----process_PIN3, eval=FALSE-------------------------------------------------
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

## ----process_PIN4, eval=FALSE-------------------------------------------------
#  path2SIF <- file.path(tempdir(), "mmusculusPIN.sif")
#  write.table(mmu_string_pin,
#              file = path2SIF,
#              col.names = FALSE,
#              row.names = FALSE,
#              sep = "\t")
#  path2SIF <- normalizePath(path2SIF)

## ----mmu_input_df-------------------------------------------------------------
knitr::kable(head(myeloma_input))

## ----run, eval=FALSE----------------------------------------------------------
#  myeloma_output <- run_pathfindR(input = myeloma_input,
#                                  convert2alias = FALSE,
#                                  gene_sets = "Custom",
#                                  custom_genes = mmu_kegg_genes,
#                                  custom_descriptions = mmu_kegg_descriptions,
#                                  pin_name_path = path2SIF)

## ----enr_chart, echo=FALSE----------------------------------------------------
enrichment_chart(myeloma_output)

## ----output-------------------------------------------------------------------
knitr::kable(myeloma_output)

## ----run2, eval=FALSE---------------------------------------------------------
#  myeloma_output <- run_pathfindR(input = myeloma_input,
#                                  convert2alias = FALSE,
#                                  gene_sets = "mmu_KEGG",
#                                  pin_name_path = "mmu_STRING")

