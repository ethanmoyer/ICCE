#' Relabel node in probe-node matrix
#'
#' Relabel node label n_old with node label n_new in either probe-node matrix 
#' (depending on the matrix parameter).
#'
#' @param icceTree icceTree data structure
#' @param n_old node number
#' @param n_new node number
#' @param matrix numerical value indicating which probe-node matrix
#' to add data to (0 == methyl states, 1 == p-values)
#' @return relabeled probe_node_matrix with node label n_old with node label
#' n_new
#' @examples
#' probe_node_matrix <- relabel_node_matrix(icceTree, 7, 15, matrix=2)
#' @export
relabel_node_matrix <- function(icceTree, n_old, n_new, matrix=2) {
	if (matrix == 0) {
		names(icceTree$probe_node_matrix)[!is.na(match(names(icceTree$probe_node_matrix), toString(n_old)))] = toString(n_new)
	} else if (icceTree$matrix == 1){
		names(icceTree$probe_node_matrix_p)[!is.na(match(names(icceTree$probe_node_matrix_p), toString(n_old)))] = toString(n_new)
	} else if (matrix == 2){
		names(icceTree$probe_node_matrix)[!is.na(match(names(icceTree$probe_node_matrix), toString(n_old)))] = toString(n_new)
		names(icceTree$probe_node_matrix_p)[!is.na(match(names(icceTree$probe_node_matrix_p), toString(n_old)))] = toString(n_new)
	}

	return(probe_node_matrix)
}

#' Relabel node on tree
#'
#' Relabel node label n_old with node label n_new on tree.
#'
#' @param icceTree icceTree data structure
#' @param n_old node number
#' @param n_new node number
#'
#' @return relabeled tree with node label n_old with node label
#' n_new
#' @examples
#' icceTree <- relabel_tree(icceTree, 7, 15)
#' @export
relabel_tree <- function(icceTree, n_old, n_new) {
	icceTree$tree$edge[icceTree$tree$edge == n_old] = n_new
	return(icceTree)
}

#' Get the leaves of a tree 
#'
#' Given a tree, this function returns its leaves.
#'
#' @param icceTree icceTree data structure
#'
#' @return a list of nodes representing the leaves of the tree
#' @examples
#' leaves <- get_leaves(icceTree)
#' @export
get_leaves <- function(icceTree) {
	ancestors = icceTree$tree$edge[, 1]
	descendants = icceTree$tree$edge[, 2]

	leaves = descendants[!(descendants %in% ancestors)]
	return(leaves)
}

#' Get number of leaves of a tree 
#'
#' Returns the number of leaves on a given tree.
#'
#' @param icceTree icceTree data structure
#'
#' @return number of leafs of a given tree
#' @examples
#' number_of_leaves = nleaves(icceTree)
#' @export
nleaves <- function(icceTree) {
	return(length(get_leaves(icceTree$tree)))
}

#' Get total number of nodes on a tree
#'
#' Returns the number of nodes on a given tree.
#'
#' @param icceTree icceTree data structure
#'
#' @return number of nodes on a tree
#' @examples
#' number_of_nodes = nnodes(icceTree)
#' @export
nnodes <- function(icceTree) {
	return(nrow(icceTree$tree$edge))
}

#' Get number of internal nodes on a tree
#'
#' Returns the number of internal nodes on the given tree.
#'
#' @param icceTree icceTree data structure
#'
#' @return number of internal nodes on a given tree
#' @examples
#' number_of_internal_nodes = ninternal_nodes(icceTree)
#' @export
ninternal_nodes <- function(icceTree) {
	return(nrow(icceTree$tree$edge) - nleaves(icceTree$tree))
}

#' Grow a tree at a given leaf node 
#'
#' Add a pair of leaves branching off of a given leaf node on a tree. 
#' Assigns the two newly added leaves labels from node_names. The old leaf 
#' where the two new branches were add is also assigned as an internal node.
#' 
#' @param icceTree icceTree data structure
#' @param leaf node number
#' @param node_names list of two node names
#'
#' @return new tree topology--two branches added to given leaf and labeled 
#' accordingly
#' @examples
#' icceTree <- grow_off_leaf(icceTree, 1, c('Fibroblasts', 'Fibrocytes'))
#' @export
grow_off_leaf <- function(icceTree, leaf, node_names) {
	tree <- icceTree$tree

	if (!(node %in% get_leaves(tree))) {
		print('Please provide a leaf as the node to which the branches will be added')
		return(tree)
	}

	max_leaf = max(tree$edge)

	tree$edge = rbind(tree$edge, c(node, max_leaf + 1), c(node, max_leaf + 2))

	tree$Nnode = tree$Nnode + 1

	tree$tip.label[number_of_get_leaves(tree)] = node_names[1]
	tree$tip.label[number_of_get_leaves(tree) + 1] = node_names[2]

	tree$tip.label = tree$tip.label[-node]

	icceTree$tree <- tree

	return(tree)
}

#' Set data on either probe-node matrix at a new or existing node
#'
#' At a given node n_new, this function adds data d to the either the methyl 
#' probe-node matrix or p-value probe-node matrix (depending on the matrix 
#' parameter). Make sure that the sizes of d and the probe-node matrix are  
#' compatible before calling this function.
#' 
#' @param icceTree icceTree data structure
#' @param n_new node number
#' @param d column matrix of methylation status for a node across many 
#' probes
#' @param matrix numerical value indicating which probe-node matrix
#' to add data to (0 == methyl states, 1 == p-values)
#'
#' @return probe-node matrix with data d stored at node n_new
#' @examples
#' tree <- add_data(icceTree, 12, probe_node_matrix['6'], matrix=0)
#' @export
add_data <- function(icceTree, n_new, d, matrix=0) {
	if (matrix == 0) {
		icceTree$probe_node_matrix[toString(n_new)] = d
	} else if (matrix == 1){
		icceTree$probe_node_matrix_p[toString(n_new)] = d
	}
	
	return(probe_node_matrix)
}

#' Get root of a tree 
#'
#' Returns the root of a given tree structure.
#' 
#' @param icceTree icceTree data structure
#'
#' @return node corresponding the the root of the tree
#' @examples
#' root <- get_root(icceTree)
#' @export
get_root <- function(icceTree) {
	tree <- icceTree$tree

	ancestor = tree$edge[, 1]
	descendant = tree$edge[, 2]

	root = ancestor[!(ancestor %in% descendant)]

	return(unique(root))
}

#' Create node matrix used for tree walking functions
#'
#' Returns a matrix which tracks where recursive functions are in a tree 
#' structure. This is passed to most recursive functions in this package.
#' 
#' @param tree ape::phylo tree
#'
#' @return node matrix with one column for ancestors and another column for 
#' descendants
#' @examples
#' node_matrix <- create_node_matrix(tree)
#' @export
create_node_matrix <- function(tree) {
	node_matrix = data.frame(ancestor=integer(), descendents=integer())

	node_matrix[1, ] = c(0, get_root(tree))

	return(node_matrix)
}


#' Recursively walk through a tree and structure remove descendants of given 
#' node.
#'
#' This function will recursively walk through a given tree and remove all of 
#' the descendants of a given node, making the given node a new leaf.
#' 
#' @param node node number
#' @param root root of non-proper subtree
#' @param node_matrix node matrix obtained from create_node_matrix function
#' @param i row in node matrix, starting at 1
#'
#' @return None
#' @examples
#' cut_tree(node, get_root(tree), node_matrix, 1)
#' @export
cut_tree <- function(node, root, node_matrix, i) {
	tree_edge = tree$edge

	ancestor = tree_edge[, 1]
	descendant = tree_edge[, 2]

  	paths = descendant[ancestor == root]

  	daughter_node = paths[1]
 	son_node = paths[2]

 	if (root %in% get_leaves(tree)) {
 		if (node %in% node_matrix[, 1]) {

  			remove_nodes = node_matrix[which(node == node_matrix[,1]):nrow(node_matrix), ]

  			for (i in 1:nrow(remove_nodes)) {
  				tree$edge <<- tree$edge[-match(remove_nodes[i, ], tree$edge)[1], ]
  				tree$tip.label <<- tree$tip.label(-root)
  				probe_node_matrix <<- probe_node_matrix[-match(toString(root), names(probe_node_matrix))]
  				probe_node_matrix_p <<- probe_node_matrix_p[-match(toString(root), names(probe_node_matrix))]
  			}
  		}
 	} else {
	    if (length(daughter_node) != 0) {

		    node_matrix[i + 1, ] = tree_edge[ancestor == root, ][1, ]

	      	cut_tree(node, daughter_node, node_matrix, i + 1)
	    }

		if (length(son_node) != 0) {

		    node_matrix[i + 1, ] = tree_edge[ancestor == root, ][2, ]

		    cut_tree(node, son_node, node_matrix, i + 1)	 
	  	}
 	}
}

#' Cut tree at specified node
#'
#' Given a tree and a selected node, this function removes the branches 
#' stemming off of that node, thus making that node a leaf. It returns this 
#' altered tree. 
#' 
#' @param icceTree icceTree data structure
#' @param node node number
#'
#' @return new tree topology--node is now a leaf and all of node's descendants 
#' are removed.
#' @examples
#' icceTree <- remove_branches_from_node(icceTree, node=8)
#' @export
remove_branches_from_node <- function(icceTree, node) {
	tree <- icceTree$tree

	node_matrix = create_node_matrix(tree)

	cut_tree(node, get_root(tree), node_matrix, 1)

	icceTree$tree <- tree

	return(icceTree)
}

