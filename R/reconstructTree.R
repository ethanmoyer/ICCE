#' Reconstruct internal nodes
#'
#' This function will recursively walk through a given tree and reconstruct 
#' all of its internal states using Fitch's algorithm
#' 
#' @param root root of non-proper subtree
#' @param node_matrix matrix obtained from create_node_matrix function
#' @param i row in node matrix, start at 1
#'
#' @return None
#' @examples
#' reconstruct_internal_nodes(get_root(tree), node_matrix, 1)
#' @export
reconstruct_internal_nodes <- function(root, node_matrix, i) {
    paths = P_descendant[P_ancestor == root]
    daughter_node = paths[1]
    son_node = paths[2]

    if (daughter_node %in% get_leaves(tree) & son_node %in% get_leaves(tree)) {

        for (j in 1:nrow(sorted_probes)) {
            probe = probe_names[j]

            daughter_node_value = sorted_probes[probe, tree$tip.label[daughter_node]]
            son_node_value = sorted_probes[probe, tree$tip.label[son_node]]

             probe_node_matrix[probe, toString(root)] <<- as.numeric(daughter_node_value == son_node_value & daughter_node_value == "g")
      
             if (probe_node_matrix[probe, toString(root)] == 0 & (daughter_node_value == "g" | son_node_value == "g")) {
                probe_node_matrix[probe, toString(root)] <<- 0.5
            }
        }
    
    } else {
        if (length(P_descendant[P_ancestor == daughter_node]) == 2) {
            node_matrix[i + 1, ] = P_edges[P_ancestor == root, ][1, ]
            reconstruct_internal_nodes(daughter_node, node_matrix, i + 1)
        }

        if (length(P_descendant[P_ancestor == son_node]) == 2) {
            node_matrix[i + 1, ] = P_edges[P_ancestor == root, ][2, ]
            reconstruct_internal_nodes(son_node, node_matrix, i + 1)
        }

        for (j in 1:nrow(sorted_probes)) {
            probe = probe_names[j]

            daughter_node_value = probe_node_matrix[probe, toString(daughter_node)]
            if (daughter_node_value == 0.5) {
                daughter_node_value = c(1,0)
            }

            son_node_value = probe_node_matrix[probe, toString(son_node)]
            if (son_node_value == 0.5) {
                son_node_value = c(1,0)
            }

            if (length(intersect(as.integer(son_node_value), as.integer(daughter_node_value))) == 1) {
                probe_node_matrix[probe, toString(root)] <<- intersect(as.integer(son_node_value), as.integer(daughter_node_value))
            } else {
                probe_node_matrix[probe, toString(root)] <<- 0.5
            }       
        }
    }
}

#' Inherit parental states to ambiguous nodes
#'
#' This function will recursively walk through a given tree and assigns 
#' parental states to ambiguous children nodes
#' 
#' @param root root of non-proper subtree
#' @param node_matrix matrix obtained from create_node_matrix function
#' @param i row in node matrix, start at 1
#'
#' @return None
#' @examples
#' inherit_parental_state(get_root(tree), node_matrix, 1)
#' @export
inherit_parental_state <- function(root, node_matrix, i) {
    paths = P_descendant[P_ancestor == root]
    daughter_node = paths[1]
    son_node = paths[2]

    for (j in 1:nprobes(consensus_vector)) {
        probe = probe_names[j]

        daughter_node_value = probe_node_matrix[probe, toString(daughter_node)]
        son_node_value = probe_node_matrix[probe, toString(son_node)]

        if (daughter_node_value == 0.5) {
            daughter_node_value = c(1,0)
        }

        if (son_node_value == 0.5) {
            son_node_value = c(1,0)
        }

        root_node_value = probe_node_matrix[probe, toString(root)]

        if (length(intersect(as.integer(son_node_value), as.integer(root_node_value))) == 1) {
            probe_node_matrix[probe, toString(son_node)] <<- intersect(as.integer(son_node_value), as.integer(root_node_value))
        }

        if (length(intersect(as.integer(daughter_node_value), as.integer(root_node_value))) == 1) {
            probe_node_matrix[probe, toString(daughter_node)] <<- intersect(as.integer(daughter_node_value), as.integer(root_node_value))
        }
    }

    if (length(P_descendant[P_ancestor == daughter_node]) == 2) {
        node_matrix[i + 1, ] = P_edges[P_ancestor == root, ][1, ]
        inherit_parental_state(daughter_node, node_matrix, i + 1)
    }

    if (length(P_descendant[P_ancestor == son_node]) == 2) {
        node_matrix[i + 1, ] = P_edges[P_ancestor == root, ][2, ]
        inherit_parental_state(son_node, node_matrix, i + 1)   
    }
}

#' Fully generate probe-node matrix with available data
#'
#' This function fully generates the probe-node matrix with a given 
#' tree. It performs ancestral tree reconstruction using Fitch's algorithm.
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#' @param tree ape::phylo tree
#'
#' @return prone-node matrix storing the methylation status across all probes 
#' and nodes on a tree
#' @examples
#' icceTree <- build_icceTree(consensus_vector, tree)
#' tree <- icceTree$tree
#' probe_node_matrix <- icceTree$probe_node_matrix
#' @export
build_icceTree <- function(consensus_vector, tree) {
    probe_node_matrix <<- initialize_probe_node_matrix(consensus_vector, tree)

    P_edges <<- tree$edge
    P_edges <<- rbind(P_edges, c(0, get_root(tree)))
    P_edges <<- P_edges[order(P_edges[,1 ], P_edges[, 2]),]
    P_ancestor <<- P_edges[, 1]
    P_descendant <<- P_edges[, 2]

    reconstruct_internal_nodes(get_root(tree), create_node_matrix(tree), 1)

    inherit_parental_state(get_root(tree), create_node_matrix(tree), 1)

    return(icce_tree(probe_node_matrix, NULL, tree))
}

#' Get UPGMA tree construction
#'
#' This function uses the UPGMA algorithm to construct a phylogenetic tree 
#' using methylation data.
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return ape::phylo tree constructed using UPGMA algorithm
#' @examples
#' tree <- upgma_tree(consensus_vector)
#' @export
upgma_tree <- function(consensus_vector) {
    library(phangorn)

    consensus_vector = t(consensus_vector)

    consensus_vector[consensus_vector == 1] = 'g'
    consensus_vector[consensus_vector == 0] = 'a'

    consensus_vector_phyDat <- phyDat(consensus_vector)

    dm <- dist.hamming(consensus_vector_phyDat) 

    tree <- upgma(dm)

    return(tree)
}

create_tree <- function() {

}

#' Initialize the probe-node matrix
#'
#' This function simply initializes a probe-node matrix across given nodes and 
#' probes with zeros.
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#' @param tree ape::phylo tree
#'
#' @return initialized probe_node_matrix
#' @examples
#' probe_node_matrix <- initialize_probe_node_matrix(consensus_vector, 
#' upgma_tree(consensus_vector))
#' @export
initialize_probe_node_matrix <- function(consensus_vector, tree) {
    internal_nodes = sort(unique(tree$edge[,1]))

    ninternal_ndoes = length(internal_nodes)

    probe_node_matrix <- data.frame()

    probe_names <- rownames(consensus_vector)
    probe_node_matrix[probe_names, 1:ninternal_ndoes] = 0

    colnames(probe_node_matrix) = internal_nodes

    for (i in 1:nleaves(tree)) {
      probe_node_matrix[toString(i)] <- as.matrix(consensus_vector[,i])
      probe_node_matrix[[toString(i)]] <- as.character(probe_node_matrix[,toString(i)])
    }

    probe_node_matrix[probe_node_matrix == "g"] = 1
    probe_node_matrix[probe_node_matrix == "a"] = 0

    return(probe_node_matrix)
}

#' Get number of probes in consensus vector
#'
#' Returns the number of probes in a list of consensus vectors
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return number of probes
#' @examples
#' number_of_probes <- nprobes(consensus_vector)
#' ... some visualization 
#' @export
nprobes <- function(consensus_vector) {
    return(nrow(consensus_vector))
}


