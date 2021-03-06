% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{get_relevant_genes}
\alias{get_relevant_genes}
\title{Find relevant genes for gain and loss of methylation}
\usage{
get_relevant_genes(icceTree, relevant_probes, reference)
}
\arguments{
\item{icceTree}{icceTree data structure}

\item{relevant_probes}{list of differentially methylated probes for each 
node}

\item{reference}{reference file}
}
\value{
list of differentially methylated genes for each node
}
\description{
This function returns a list of relevant gene for each node that exhibit a 
gain or loss in methylation status.
}
\examples{
relevant_genes <- get_relevant_genes(icceTree, relevant_probes, reference)
}
