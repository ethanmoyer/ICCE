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
        # Replace 'node' with 'j'?
        for (j in 1:nrow(probe_node_matrix)) {
            probe = rownames(probe_node_matrix)[j]
            daughter_node_value = consensus_state[probe, tree$tip.label[daughter_node]]
            son_node_value = consensus_state[probe, tree$tip.label[son_node]]

             probe_node_matrix[probe, toString(root)] <<- as.numeric(daughter_node_value == son_node_value & daughter_node_value == 1)
             if (probe_node_matrix[probe, toString(root)] == 0 & (daughter_node_value == 1 | son_node_value == 1)) {
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

        for (j in 1:nrow(probe_node_matrix)) {
            probe = rownames(probe_node_matrix)[j]

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

    for (j in 1:nrow(consensus_state)) {
        probe = rownames(consensus_state)[j]

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

#' Run a simulation of Fitch's algorithm for all nodes
#'
#' This function runs n simulations in which is randomizes methylation calls 
#' of probes according to their p-value and calls Fitch's algorithm. It 
#' assigns the p-value of all internal nodes based on the number of 
#' observations of methylation in n simulations.
#' 
#' @param probe_node_matrix probe-node matrix of methylation calls
#' @param probe_node_matrix_p probe-node matrix of p-values
#'
#' @return m probe vector of methylation calls for all node
#' @return p probe vector of p-values for all node
#' @examples
#' simulate_fitch_info <- reconstruct_internal_nodes_p(probe_node_matrix, 
#' probe_vector <- simulate_fitch_info$m
#' probe_vector_p <- simulate_fitch_info$p 
#' @export
simulate_fitch <- function(tree, consensus_state, probe_node_matrix, probe_node_matrix_p, n = 1000) {
    probe_node_matrix_ <- probe_node_matrix
    m <- initialize_probe_node_matrix(consensus_state, tree)
    p <- probe_node_matrix_p
    for (k in 1:n) {
        for (i in 1:nrow(probe_node_matrix)) {
            for (j in get_leaves(tree)) {
                if (runif(1, 0, 1) < probe_node_matrix_p[i, toString(j)]) {
                    probe_node_matrix[i, toString(j)] <- bitwXor(as.integer(probe_node_matrix[i, toString(j)]), 1)
                } 
            }
        }

        reconstruct_internal_nodes(tree_root, create_node_matrix(tree), 1)
        inherit_parental_state(tree_root, create_node_matrix(tree), 1)

        for (i in 1:nrow(probe_node_matrix)) {
            for (j in get_internal_nodes(tree)) {
                p[i, toString(j)] <- p[i, toString(j)] + as.integer(probe_node_matrix[i, toString(j)] == 1 | probe_node_matrix[i, toString(j)] == "1.0")
            }
        }
        probe_node_matrix <- probe_node_matrix_
    }
    
    for (i in 1:nrow(probe_node_matrix)) {
        for (j in get_internal_nodes(tree)) {
            mi = p[i, toString(j)] / n

            if (mi > 0.5) {
                p[i, toString(j)] <- 1 - mi
                m[i, toString(j)] <- 1
            } else {
                p[i, toString(j)] <- mi
                m[i, toString(j)] <- 0
           }
        }
    }

    return(list(m=m, p=p))

}

#' Get UPGMA tree construction
#'
#' This function uses the UPGMA algorithm to construct a phylogenetic tree 
#' using methylation data.
#' 
#' @param consensus_state consensus vector data frame for each cell in group 
#' names
#'
#' @return ape::phylo tree constructed using UPGMA algorithm
#' @examples
#' tree <- upgma_tree(consensus_state)
#' @export
upgma_tree <- function(consensus_state) {
    library(phangorn)

    consensus_state = t(consensus_state)

    consensus_state[consensus_state == 1] = 'g'
    consensus_state[consensus_state == 0] = 'a'

    consensus_state_phyDat <- phyDat(consensus_state)

    dm <- dist.hamming(consensus_state_phyDat) 

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
#' @param consensus_state consensus vector data frame for each cell in group 
#' names
#' @param tree ape::phylo tree
#'
#' @return initialized probe_node_matrix and probe_names
#' @examples
#' probe_node_matrix <- initialize_probe_node_matrix(consensus_state, 
#' upgma_tree(consensus_state))
#' @export
initialize_probe_node_matrix <- function(consensus_state, tree) {
    internal_nodes = sort(unique(tree$edge[,1]))

    ninternal_ndoes = length(internal_nodes)

    probe_node_matrix <- data.frame()

    probe_names <- rownames(consensus_state)
    probe_node_matrix[probe_names, 1:ninternal_ndoes] = 0

    colnames(probe_node_matrix) = internal_nodes

    for (i in 1:nleaves(tree)) {
      probe_node_matrix[toString(i)] <- as.matrix(consensus_state[,i])
      probe_node_matrix[[toString(i)]] <- as.character(probe_node_matrix[,toString(i)])
    }

    return(probe_node_matrix)
}

#' Get number of probes in consensus vector
#'
#' Returns the number of probes in a list of consensus vectors
#' 
#' @param consensus_state consensus vector data frame for each cell in group 
#' names
#'
#' @return number of probes
#' @examples
#' number_of_probes <- nprobes(consensus_state)
#' ... some visualization 
#' @export
nprobes <- function(consensus_state) {
    return(nrow(consensus_state))
}

#' Fully generate probe-node matrix with available data
#'
#' This function fully generates the probe-node matrix with a given 
#' tree. It performs ancestral tree reconstruction using Fitch's algorithm.
#' 
#' @param consensus_state consensus methylation call vector data frame for each cell in group names
#' @param consensus_p consensus p-value vector data frame for each cell in group names
#' @param tree ape::phylo tree
#'
#' @return prone-node matrix storing the methylation status across all probes 
#' and nodes on a tree
#' @examples
#' icceTree <- build_icceTree(consensus_state, tree)
#' tree <- icceTree$tree
#' probe_node_matrix <- icceTree$probe_node_matrix
#' @export
build_icceTree <- function(consensus_state, consensus_p, tree) {
    probe_node_matrix <<- initialize_probe_node_matrix(consensus_state, tree)

    tree_root <<- get_root(tree)

    P_edges <<- tree$edge
    P_edges <<- rbind(P_edges, c(0, tree_root))
    P_edges <<- P_edges[order(P_edges[,1 ], P_edges[, 2]),]
    P_ancestor <<- P_edges[, 1]
    P_descendant <<- P_edges[, 2]

    reconstruct_internal_nodes(tree_root, create_node_matrix(tree), 1)

    inherit_parental_state(tree_root, create_node_matrix(tree), 1)

    probe_node_matrix_copy <<- probe_node_matrix

    probe_node_matrix_p <<- probe_node_matrix
    for (i in get_nodes(tree)) {
        if (i %in% get_leaves(tree)) {
            probe_node_matrix_p[toString(i)] = consensus_p[i]
        } else {
            probe_node_matrix_p[toString(i)] = integer(nrow(consensus_p))
        }
    }

    simulate_fitch_info <- simulate_fitch(tree, consensus_state, probe_node_matrix_copy, probe_node_matrix_p)
    probe_node_matrix <- simulate_fitch_info$m
    probe_node_matrix_p <- simulate_fitch_info$p

    return(icce_tree(probe_node_matrix, probe_node_matrix_p, tree))
}

