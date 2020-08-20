# ICCE - Iterative Construction of Cyto-epigenetic Evolution

NEED TO UPDATE

To install from github {Not currently implemented},
```R
install.packages("BiocManager")
library(BiocManager)
install("ethanmoyer/icce")
```

Just clone to install for now.

ICCE can process Illumina Infinium DNA methylation data. It currently only supports the EPIC Platform. 

## Install ICCE

Development version can be installed from github
```{r, eval=FALSE}
BiocManager::install('ethanmoyer/icce')
BiocManager::install('ethanmoyer/icceData')
```

## Creating consensus vectors with the ICE Pipeline

The ICCE data package contains a small reference file that can be obtained through:
```{r eval=FALSE, message = FALSE, warning = FALSE}
ref_file <- system.file("extdata/reference.rds", package = "icceData")
reference <- get_reference(ref_file)
```

ICCE utilizes the openSesame pipeline through:
```{r eval=FALSE, message = FALSE, warning = FALSE}
idat_dir <- system.file("extdata/", package = "sesameData")
betas <- read_IDAT(idat_dir)
```
where `idat_dir` is the directory containing all the IDAT files (they can be
present under nested sub-directories).

Some data may have groups that are heterogeneous or undesirable for tree reconstruction. These can be removed by simply including their sample names in the following variable:
```{r eval=FALSE, essage = FALSE, warning = FALSE}
remove_names = c("PCA0612", "WB1148", "A2.mix", "B6.mix", "PCA1547", "WB1147", "B3.mix", "B2.mix", "A4.mix", "B4.mix", "A3.mix", "ST2007", "B5.mix", "B1.mix", "A5.mix", "WB1145", "A6.mix", "A1.mix")
```
The betas generated earlier from openSesame are compressed into consensus vectors for each of the indicated cell types:
```{r eval=FALSE, message = FALSE, warning = FALSE}
consensus_vector_names = c("NK", "Mono", "Tc", "Neu", "Th", "B")
series_matrix_file <-  system.file("extdata/series_matrix.rds", package = "sesameData")
consensus_vector <- betas_to_consensus_vector(reference, series_matrix_file, betas, consensus_vector_names)
```

### Why are certain probes removed?

The openSesame pipeline sometimes replaces beta values with NA because of certain probe masking in the pipeline. These are masked because of their low quality during the Illumina assay. 

Before the consensus vectors are built for each of the specified groups, the beta values are fit to a different bimodal distribution for each cell type. Several considerations about this distribution can explain the removal of certain probes. First, probes located on sex chromosomes are removed prior to creating the distributions because they have a disproportionate global methylation status. Also, probes that have a beta value outside of two standard deviations away from either of the two means in the distribution are labeled intermediates and removed.

In addition to these parameters, certain probes exhibit zero parsimony. Including them in tree reconstruction would be redundant as they do not add to the tree topology. Therefore, they are also removed.

## Reconstructing a tree using the ICCE Pipeline

A tree can be constructed using the UPGMA algorithm with the following:
```{r eval=FALSE, message = FALSE, warning = FALSE}
tree <- upgma_tree(consensus_vector)
```

To load a known tree into the internal format of the pipeline, please use:
```{r eval=FALSE, message = FALSE, warning = FALSE}
tree <- load_tree()
```

Using the consensus vectors, the internal states of the tree can be reconstructed in the following way:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- build_icceTree(consensus_vector, tree)
```

### How are internal states reconstructed?

Using a recursive implementation of Fitch's algorithm, all of the internal states are inferred based off of the nodes of the tree and assigned as either a single state (methylated or unmethylated) or as an ambiguous state (methylated and unmethylated). In a second pass, those internal ambiguous states are resolved by inheriting parental unambiguous states to children.

### What is the structure of an icceTree object?

An icceTree object is composed of three different structures. The first structure is a n probe by m node data frame (prone-node matrix) that houses the methylation states of all of the nodes on the tree across all probes. A methylated status is assigned 1, unmethylated is 0, and ambiguous is 0.5. The second structure is also a data frame of the same size. This probe-node matrix houses p-values for each of the nodes across all probes. These are used for iterative reconstruction of the cellular tree when more and more data is added. the third structure is a phylogenetic tree structure from ape. 

## Editing a tree using the ICCE package

Data can be added to either one of the two probe-node matrices through:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- add_data(icceTree, n_new=12, d=icceTree icceTree$probe_node_matrix['6'], matrix=0)
```
where matrix is a single value indicating to which probe-node matrix data will be added (0 == methyl states, 1 == p-values).

Either probe-node matrix can be relabeled with the following:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- relabel_node_matrix(icceTree, n_old=7, old_new=15, matrix=2)
```
where matrix is a single value indicating to which probe-node matrix data will be added (0 == methyl states, 1 == p-values, 2 == both).

The phylo tree stored in the icceTree data structure can also be relabeled:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- relabel_tree(icceTree, n_old=7, n_new=15)
```

Certain parameters about the tree structure can also be obtained, such as a list of leaves, the number of leaves, the total number of nodes, and the number of internal nodes:
```{r eval=FALSE, message = FALSE, warning = FALSE}
leaves <- get_leaves(icceTree)
number_of_leaves = nleaves(icceTree)
number_of_nodes = nnodes(icceTree)
number_of_internal_nodes = ninternal_nodes(icceTree)
```

Also the root of the tree can be obtained with:
```{r eval=FALSE, message = FALSE, warning = FALSE}
root <- get_root(icceTree)
```

A tree can be grown from any current leaf on the tree structure with:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- grow_off_leaf(icceTree, leaf=1, node_names=c('Fibroblasts', 'Fibrocytes'))
```
where node_names contains the labels for the two new leaves add at the first node.

The tree can also be cut at any node using:
```{r eval=FALSE, message = FALSE, warning = FALSE}
icceTree <- remove_branches_from_node(icceTree, node=8)
```
where the 8th node will become a new leaf and all of its descendants will be removed from the tree structure. Descendants are also automatically removed from both the methyl probe-node matrix and the p-value probe-node matrix.

## Analyzing the reconstructed tree using ICCE pipeline

Differentially methylated probes on each branch can be obtained using:
```{r eval=FALSE, message = FALSE, warning = FALSE}
probe_information <- generate_branch_changes(icceTree)
relevant_probes <- probe_information$relevant_probes
```
The probe_information variable also has three other structures in it: a node label list for the total number of changes in methylation, number of gains in methylation, and the number of losses in methylation along each branch on the tree: 
```{r eval=FALSE, message = FALSE, warning = FALSE}
total_changes_label <- probe_information$total_changes_label
meth_gains_label <- probe_information$meth_gains_label
meth_losses_label <- probe_information$meth_losses_label
```

Using these relevant probes, corresponding differentially methylated genes can be obtained with:
```{r eval=FALSE, message = FALSE, warning = FALSE}
relevant_genes <- get_relevant_genes(icceTree, relevant_probes, reference)
```

One powerful tool of this package is the ability to track the differentiation on the tree of probes along a selected gene:
```{r eval=FALSE, message = FALSE, warning = FALSE}
gene_information <- track_gene(icceTree, reference, gene='CD81')
```
where the output are two node label lists with the proportion of methylated internal nodes and leaves at each of the three states (methylated, unmethylated, and ambiguous).
```{r eval=FALSE, message = FALSE, warning = FALSE}
thermo_prop_internal_nodes = gene_information$thermo_prop_internal_nodes
thermo_prop_leaves = gene_information$thermo_prop_leaves
```

## Visualizing trees using the ICCE package

A reconstructed tree can be visualized with the following command:
```{r eval=FALSE, message = FALSE, warning = FALSE}
show_tree(filename='tree.pdf', tree=icceTree$tree, title='Hematopoietic Tree', edge_labels=meth_gains_label)
```
where edge_labels is an optional parameter.

The edge labels generated for a specific gene from the track_gene function can be visualized as well:
```{r eval=FALSE, message = FALSE, warning = FALSE}
show_gene_track('tree_CD81', tree=icceTree$tree, gene='CD81', thermo_prop_internal_nodes=thermo_prop_internal_nodes, thermo_prop_leaves=thermo_prop_leaves)
```

Another useful plot is a heat map of probes with high variance:
```{r eval=FALSE, message = FALSE, warning = FALSE}
in progress...
head_map <- show_heatmap()
```
