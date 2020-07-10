#' Create icce tree object
#'
#' Returns icce tree object given a probe-node matrix, a probe-node p-values 
#' and a ape::phylo tree
#' 
#' @param probe_node_matrix probe-node matrix storing the methylation status 
#' across all probes and nodes on a tree
#' @param probe_node_matrix_p probe-node matrix storing the p-value of each 
#' node across all probes on a tree
#' @param tree ape::phylo tree
#'
#' @return icce tree object
#' @examples
#' reconstruct_internal_nodes(get_root(tree), node_matrix, 1)
#' @export
icce_tree <- function(probe_node_matrix, probe_node_matrix_p, tree) {
    return(list(probe_node_matrix=probe_node_matrix, probe_node_matrix_p = probe_node_matrix_p, tree=tree))
}
