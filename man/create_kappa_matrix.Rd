% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{create_kappa_matrix}
\alias{create_kappa_matrix}
\title{Create Kappa Statistics Matrix}
\usage{
create_kappa_matrix(
  enrichment_res,
  use_description = FALSE,
  use_active_snw_genes = FALSE
)
}
\arguments{
\item{enrichment_res}{data frame of pathfindR enrichment results. Must-have
columns are 'Term_Description' (if \code{use_description = TRUE}) or 'ID'
(if \code{use_description = FALSE}), 'Down_regulated', and 'Up_regulated'.
If \code{use_active_snw_genes = TRUE}, 'non_Signif_Snw_Genes' must also be
provided.}

\item{use_description}{Boolean argument to indicate whether term descriptions
(in the 'Term_Description' column) should be used. (default = \code{FALSE})}

\item{use_active_snw_genes}{boolean to indicate whether or not to use
non-input active subnetwork genes in the calculation of kappa statistics
(default = FALSE, i.e. only use affected genes)}
}
\value{
a matrix of kappa statistics between each term in the
enrichment results.
}
\description{
Create Kappa Statistics Matrix
}
\examples{
sub_df <- example_pathfindR_output[1:3, ]
create_kappa_matrix(sub_df)
}
