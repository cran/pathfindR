## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## -----------------------------------------------------------------------------
# ## the default organism is "Homo_sapiens"
# path_to_pin_file <- get_pin_file()

## -----------------------------------------------------------------------------
# ## retrieving PIN data for "Gallus_gallus"
# path_to_pin_file <- get_pin_file(org = "Gallus_gallus")

## -----------------------------------------------------------------------------
# ## saving the "Homo_sapiens" PIN as "/path/to/PIN/file"
# path_to_pin_file <- get_pin_file(path2pin = "/path/to/PIN/file")

## -----------------------------------------------------------------------------
# ## retrieving PIN data for "Mus_musculus" from BioGRID release 3.5.179
# path_to_pin_file <- get_pin_file(
#   org = "Mus_musculus",
#   release = "3.5.179"
# )

## -----------------------------------------------------------------------------
# ## obtaining KEGG pathway gene sets for Rattus norvegicus (rno)
# gsets_list <- get_gene_sets_list(org_code = "rno")

## -----------------------------------------------------------------------------
# gsets_list <- get_gene_sets_list(source = "Reactome")

## -----------------------------------------------------------------------------
# gsets_list <- get_gene_sets_list(
#   source = "MSigDB",
#   collection = "C2"
# )

## -----------------------------------------------------------------------------
# ## obtaining C5 gene sets data for "Drosophila melanogaster"
# gsets_list <- get_gene_sets_list(
#   source = "MSigDB",
#   species = "Drosophila melanogaster",
#   collection = "C5"
# )

## ----eval=TRUE----------------------------------------------------------------
## see msigdbr::msigdbr_species() for all available organisms
msigdbr::msigdbr_species()

## -----------------------------------------------------------------------------
# ## obtaining C3 - MIR: microRNA targets
# gsets_list <- get_gene_sets_list(
#   source = "MSigDB",
#   collection = "C3",
#   subcollection = "MIR"
# )

