#' Display tree
#'
#' Displays constructed tree and saves it at the given filename.
#'
#' @param filename filename to which the plot will be saved
#' @param tree ape::phylo tree
#' @param title title name
#' @param edge_labels optional labels for each branch
#' @return None
#' @examples
#' show_tree(filename='tree.pdf', tree=icceTree$tree, title='Hematopoietic Tree', 
#' edge_labels=meth_gains_label) 
#' @export
show_tree <- function(filename, tree, title, edge_labels) {
	pdf(filename)
	plot(tree)
	title(title)
	if(!missing(edge_labels)) {
		edgelabels(edge_labels)
	}
	nodelabels()
	dev.off()
}

#' Track the single probe resolution change of a gene
#'
#' Displays constructed tree with node labels for the proportion of each probe 
#' state for a selected gene. This is saved at the given filename.
#'
#' @param filename filename to which the plot will be saved
#' @param tree ape::phylo tree
#' @param gene gene name
#' @param thermo_prop_node thermo proportion plot for internal nodes of all of 
#' the probe states
#' @param thermo_prop_leaves thermo proportion plot for leaves of all of
#' the probe states
#' @return None
#' @examples
#' show_gene_track('tree_CD81', tree=icceTree$tree, gene='CD81', 
#' thermo_prop_internal_nodes=thermo_prop_internal_nodes, thermo_prop_leaves=
#' thermo_prop_leaves)
#' ... some visualization 
#' @export
show_gene_track <- function(filename, tree, gene, thermo_prop_internal_nodes, thermo_prop_leaves) {
	pdf(filename)
	plot(tree)
	title(paste("Summary plot for gene: ", gene, sep=""))
	nodelabels(thermo = thermo_prop_internal_nodes, piecol=topo.colors(3))
	tiplabels(adj = c(0.59,0.5), thermo = thermo_prop_leaves, piecol=topo.colors(3))
	legend("topleft", inset=.02, title="Methylation Status", c("Methylated", "Unmethylated", "Ambiguous"), fill=topo.colors(3), horiz=TRUE, cex=0.8)
	dev.off()
}

#' Display heatmap of probes {In progress}
#'
#' ...
#'
#' @param None
#' @return 
#' @examples
#' head_map <- show_heatmap()
#' @export
show_heatmap <- function(betas, reference) {


jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

library(gplots)

pdf('temp_repo/hm12.pdf')
heatmap.2(matrix, col=jet.colors, density.info="none", trace="none", dendrogram='none')
dev.off()

}

