% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/editTree.R
\name{add_data}
\alias{add_data}
\title{Set data on either probe-node matrix at a new or existing node}
\usage{
add_data(icceTree, n_new, d, matrix = 0)
}
\arguments{
\item{icceTree}{icceTree data structure}

\item{n_new}{node number}

\item{d}{column matrix of methylation status for a node across many 
probes}

\item{matrix}{numerical value indicating which probe-node matrix
to add data to (0 == methyl states, 1 == p-values)}
}
\value{
probe-node matrix with data d stored at node n_new
}
\description{
At a given node n_new, this function adds data d to the either the methyl 
probe-node matrix or p-value probe-node matrix (depending on the matrix 
parameter). Make sure that the sizes of d and the probe-node matrix are  
compatible before calling this function.
}
\examples{
tree <- add_data(icceTree, 12, probe_node_matrix['6'], matrix=0)
}
