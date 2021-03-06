% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{generate_branch_changes}
\alias{generate_branch_changes}
\title{Generate list of relevant probes for gain, loss, and change in methylation 
status}
\usage{
generate_branch_changes(icceTree)
}
\arguments{
\item{icceTree}{icceTree data structure}

\item{probe_node_matrix}{probe-node matrix storing the methylation status 
across all probes and nodes on a tree}
}
\value{
a structure containing a list of probes with differentially 
methylated probes for each node (gain and loss stored separately), a list 
of branch labels for the number of total changes in methylation, a list of 
branch labels for the number of gains in methylation, and a list of branch 
labels for the number of losses in methylation
}
\description{
This function returns relevant probe lists for each node that exhibits 
either a gain or a loss in methylation status. It also returns three sets 
of branch labels: number of methylation changes, number of methylation 
gains, and number of methylation losses. These labels can be given to the 
show_tree function as an edge labels parameter to visualize this data.
}
\examples{
probe_information <- generate_branch_changes(icceTree)
relevant_probes <- probe_information$relevant_probes
total_changes_label <- probe_information$total_changes_label
meth_gains_label <- probe_information$meth_gains_label
meth_losses_label <- probe_information$meth_losses_label
}
