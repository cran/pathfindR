% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{download_KGML_file}
\alias{download_KGML_file}
\title{Obtain KGML file for a KEGG pathway (hsa)}
\usage{
download_KGML_file(pw_id, pwKGML, quiet = TRUE)
}
\arguments{
\item{pw_id}{KEGG pathway ID}

\item{pwKGML}{destination file}

\item{quiet}{If \code{TRUE}, suppress status messages (if any), and the
progress bar while downloading file(s)}
}
\value{
download status (0 for success), if warning/error returns NA
}
\description{
Obtain KGML file for a KEGG pathway (hsa)
}
