% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{run_pathfindR}
\alias{run_pathfindR}
\title{Wrapper Function for pathfindR - Active-Subnetwork-Oriented Enrichment Workflow}
\usage{
run_pathfindR(
  input,
  gene_sets = "KEGG",
  min_gset_size = 10,
  max_gset_size = 300,
  custom_genes = NULL,
  custom_descriptions = NULL,
  pin_name_path = "Biogrid",
  p_val_threshold = 0.05,
  enrichment_threshold = 0.05,
  convert2alias = TRUE,
  plot_enrichment_chart = TRUE,
  output_dir = NULL,
  list_active_snw_genes = FALSE,
  ...
)
}
\arguments{
\item{input}{the input data that pathfindR uses. The input must be a data
  frame with three columns: \enumerate{
  \item Gene Symbol (Gene Symbol)
  \item Change value, e.g. log(fold change) (OPTIONAL)
  \item p value, e.g. adjusted p value associated with differential expression
}}

\item{gene_sets}{Name of the gene sets to be used for enrichment analysis.
Available gene sets are 'KEGG', 'Reactome', 'BioCarta', 'GO-All',
'GO-BP', 'GO-CC', 'GO-MF', 'cell_markers', 'mmu_KEGG' or 'Custom'.
If 'Custom', the arguments \code{custom_genes} and \code{custom_descriptions}
must be specified. (Default = 'KEGG')}

\item{min_gset_size}{minimum number of genes a term must contain (default = 10)}

\item{max_gset_size}{maximum number of genes a term must contain (default = 300)}

\item{custom_genes}{a list containing the genes involved in each custom
term. Each element is a vector of gene symbols located in the given custom
term. Names should correspond to the IDs of the custom terms.}

\item{custom_descriptions}{A vector containing the descriptions for each
custom  term. Names of the vector should correspond to the IDs of the custom
terms.}

\item{pin_name_path}{Name of the chosen PIN or absolute/path/to/PIN.sif. If PIN name,
must be one of c('Biogrid', 'STRING', 'GeneMania', 'IntAct', 'KEGG', 'mmu_STRING'). If
path/to/PIN.sif, the file must comply with the PIN specifications. (Default = 'Biogrid')}

\item{p_val_threshold}{the p value threshold to use when filtering
the input data frame. Must a numeric value between 0 and 1. (default = 0.05)}

\item{enrichment_threshold}{adjusted-p value threshold used when filtering
enrichment results (default = 0.05)}

\item{convert2alias}{boolean to indicate whether or not to convert gene symbols
in the input that are not found in the PIN to an alias symbol found in the PIN
(default = TRUE) IMPORTANT NOTE: the conversion uses human gene symbols/alias symbols.}

\item{plot_enrichment_chart}{boolean value. If TRUE, a bubble chart displaying
the enrichment results is plotted. (default = TRUE)}

\item{output_dir}{the directory to be created where the output and intermediate
files are saved (default = \code{NULL}, a temporary directory is used)}

\item{list_active_snw_genes}{boolean value indicating whether or not to report
the non-significant active subnetwork genes for the active subnetwork which was enriched for
the given term with the lowest p value (default = \code{FALSE})}

\item{...}{additional arguments for \code{\link{active_snw_enrichment_wrapper}}}
}
\value{
Data frame of pathfindR enrichment results. Columns are: \describe{
  \item{ID}{ID of the enriched term}
  \item{Term_Description}{Description of the enriched term}
  \item{Fold_Enrichment}{Fold enrichment value for the enriched term (Calculated using ONLY the input genes)}
  \item{occurrence}{the number of iterations that the given term was found to enriched over all iterations}
  \item{support}{the median support (proportion of active subnetworks leading to enrichment within an iteration) over all iterations}
  \item{lowest_p}{the lowest adjusted-p value of the given term over all iterations}
  \item{highest_p}{the highest adjusted-p value of the given term over all iterations}
  \item{non_Signif_Snw_Genes (OPTIONAL)}{the non-significant active subnetwork genes, comma-separated}
  \item{Up_regulated}{the up-regulated genes (as determined by `change value` > 0, if the `change column` was provided) in the input involved in the given term's gene set, comma-separated. If change column not provided, all affected are listed here.}
  \item{Down_regulated}{the down-regulated genes (as determined by `change value` < 0, if the `change column` was provided) in the input involved in the given term's gene set, comma-separated}
}
 The function also creates an HTML report with the pathfindR enrichment
 results linked to the visualizations of the enriched terms in addition to
 the table of converted gene symbols. This report can be found in
 '\code{output_dir}/results.html' under the current working directory.

 By default, a bubble chart of top 10 enrichment results are plotted. The x-axis
 corresponds to fold enrichment values while the y-axis indicates the enriched
 terms. Sizes of the bubbles indicate the number of significant genes in the given terms.
 Color indicates the -log10(lowest-p) value; the more red it is, the more
 significant the enriched term is. See \code{\link{enrichment_chart}}.
}
\description{
\code{run_pathfindR} is the wrapper function for the pathfindR workflow
}
\details{
This function takes in a data frame consisting of Gene Symbol, log-fold-change
and adjusted-p values. After input testing, any gene symbols that are not in
the PIN are converted to alias symbols if the alias is in the PIN. Next,
active subnetwork search is performed. Enrichment analysis is
performed using the genes in each of the active subnetworks. Terms with
adjusted-p values lower than \code{enrichment_threshold} are discarded. The
lowest adjusted-p value (over all subnetworks) for each term is kept. This
process of active subnetwork search and enrichment is repeated  for a selected
number of \code{iterations}, which is done in parallel. Over all iterations,
the lowest and the highest adjusted-p values, as well as number of occurrences
are reported for each enriched term.
}
\section{Warning}{
 Especially depending on the protein interaction network,
 the algorithm and the number of iterations you choose, 'active subnetwork
 search + enrichment' component of \code{run_pathfindR} may take a long time to finish.
}

\examples{
\dontrun{
run_pathfindR(example_pathfindR_input)
}
}
\seealso{
\code{\link{input_testing}} for input testing, \code{\link{input_processing}} for input processing,
\code{\link{active_snw_search}} for active subnetwork search and subnetwork filtering,
\code{\link{enrichment_analyses}} for enrichment analysis (using the active subnetworks),
\code{\link{summarize_enrichment_results}} for summarizing the active-subnetwork-oriented enrichment results,
\code{\link{annotate_term_genes}} for annotation of affected genes in the given gene sets,
\code{\link{visualize_terms}} for visualization of enriched terms,
\code{\link{enrichment_chart}} for a visual summary of the pathfindR enrichment results,
\code{\link[foreach]{foreach}} for details on parallel execution of looping constructs,
\code{\link{cluster_enriched_terms}} for clustering the resulting enriched terms and partitioning into clusters.
}
