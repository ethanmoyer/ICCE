#' Load betas data
#'
#' Using sesame, this function loads betas data using IDAT files stored in a 
#' given directory
#' 
#' @param IDAT_dir directory location of where given IDAT file is located
#'
#' @return betas matrix
#' @examples
#' betas <- read_IDAT('~/20191212_GEO_datasets/GSE110554')
#' @export
read_IDAT <- function(IDAT_dir) {
	library(sesame)
	
	pfx = searchIDATprefixes(IDAT_dir)
	betas = openSesame(pfx)

	return(betas)
}

#' Removes probes on sex chromosomes 
#'
#' This function removes probes located on sex chromosomes because these 
#' probes have a disproportionate global methylation status
#' 
#' @param reference reference file
#' @param betas matrix (rows are probe names, columns are cell names)
#'
#' @return betas matrix without probes located on sex chromosomes.
#' @examples
#' betas <- remove_sex_chr(reference, betas)
#' @export
remove_sex_chr <- function(reference, betas) {
	sex_probes_remove = rownames(reference)[reference[,"seqnames"] == "chrX" | reference[,"seqnames"] == "chrY" | reference[,"seqnames"] == "*"]

	betas = betas[,!(colnames(betas) %in% sex_probes_remove)]
	return(betas)
}

#' Return reference data frame 
#'
#' This function returns the contents of a given reference file as a data 
#' frame. The rows of the data frame are the probe names, allowing for quick
#' indexing.
#' 
#' @param reference_file file location of reference file
#'
#' @return reference file
#' @examples
#' reference <- get_reference('~/references/hg38/annotation/R/
#' EPIC.hg38.manifest.rds')
#' @export
get_reference <- function(reference_file) {
	reference <- readRDS(reference_file)
	probe_names <- names(reference)
	reference <- data.frame(reference)
	rownames(reference) <- probe_names

	return(reference)
}

#' Removes specific cell types from betas matrix 
#'
#' This function removes specified group names. Common names include MIX or 
#' are typically heterogeneous mixtures of cells.
#' 
#' @param betas matrix (rows are cell names, columns are probe names)
#' @param remove_names list of cell names in betas
#'
#' @return betas matrix without specified cell names
#' @examples
#' betas <- remove_cell_types(betas, c('A.mix', 'B.mix', ...))
#' @export
remove_cell_types <- function(betas, remove_names) {
	betas <- betas[!(rownames(betas) %in% remove_names),]
	return(betas)
}

#' Generate bimodal distribution on cell types
#'
#' This function generates a bimodal distribution on the betas matrix data 
#' within each cell type
#' 
#' @param betas matrix (rows are cell names, columns are probe names)
#'
#' @return matrix specifying mu0, mu1, sigma0, sigma1 across all cell types
#' @examples
#' bi <- generate_bimodal_distribution(betas)
#' @export
generate_bimodal_distribution <- function(betas) {
	library(BimodalIndex)

	b_data <- betas[, apply(betas, 2, function(x) !any(is.na(x)))]
	bi <- bimodalIndex(b_data, verbose=FALSE)
	return(bi)
}

#' Categorize betas based on bimodal distribution 
#'
#' This function uses a given bimodal distribution generated within cell types
#' in order to categorize betas as either methylated (1), unmethylated (0), or 
#' ambiguous (-1).
#' 
#' @param betas matrix (rows are cell names, columns are probe names)
#' @param bi matrix specifying mu0, mu1, sigma0, sigma1 across all cell types
#'
#' @return betas matrix categorized based on bimodal distribution
#' @examples
#' betas <- fit_betas_to_distribution(betas, bi)
#' @export
fit_betas_to_distribution <- function(betas, bi) {
	for(i in 1:nrow(betas)) {
		for(k in 1:ncol(betas)) {
	    	if (is.na(betas[i,k])) {              
	   		} else if (betas[i,k] <= (bi[i,1] + 2*bi[i,3])) { 
	    		betas[i,k] = 0
	   		} else if (betas[i,k] >= (bi[i,2] - 2*bi[i,3])) { 
	   			betas[i,k] = 1
	   		} else {
	   			betas[i,k] = -1                   
	  		}
	  	}
	}
	return(betas)
}

#' Create consensus vector of probes for each group
#'
#' This function creates consensus vector for each group using majority ruling:
#' across all probes, whichever category (methylated/unmethylated) is 
#' represented in more than 2/3's of the sample will saved. Otherwise, -1 for 
#' ambiguous is saved. A probabilistic model will supplement this in the 
#' future.
#' 
#' @param betas matrix (rows are cell names, columns are probe names)
#' @param group_names list of cell group names in betas matrix
#'
#' @return consensus vector for each cell in group names
#' @examples
#' consensus_vector <- create_consensus_vector(betas, c('NK', "Monocytes", 
#' ...))
#' @export
create_consensus_vector <- function(betas, group_names) {
	consensus_vector <- c()

	for (i in 1:length(group_names)) {
		name = group_names[i]

		select_group = betas[grepl(name, rownames(betas)), ]

		n = nrow(select_group)

		S = select_group[1,]
		S[colSums(select_group == 1) > (n - n / 3)] = 1
		S[colSums(select_group == 0) > (n - n / 3)] = 0
		S[colSums(select_group == 0) <= (n - n / 3) & colSums(select_group == 1) <= (n - n / 3)] = -1

		consensus_vector[[name]] = S
	} 

	return(data.frame(consensus_vector))
}

#' Convert betas matrix to consensus vectors of probes for each cell group 
#'
#' This function creates a consensus vector for each group using the betas 
#' matrix. Sex cells and unwanted cell groups given the function are removed 
#' first. Then the cells are fit to a bimodal distribution and finally 
#' compiled into consensus vectors for each cell group.
#' 
#' @param reference reference file
#' @param series_matrix_file directory location of where given series matrix 
#' file is located
#' @param betas matrix (rows are cell names, columns are probe names)
#' @param group_names list of cell names in betas
#' @param remove_names list of cell names in betas
#'
#' @return consensus vector for each cell in group names
#' @examples
#' consensus_vector <- betas_to_consensus_vector(reference, 'GSE110554/samples_GEOLoadSeriesMatrix.rds', 
#' betas, c('NK', "Monocytes", ...), c('A.mix', 'B.mix', ...))
#' @export
betas_to_consensus_vector <- function(reference_file, series_matrix_file, betas, group_names, remove_names) {

	betas <- remove_sex_chr(reference_file, betas)
	save_data(betas, 'betas_raw')

	series_matrix <- readRDS(series_matrix_file)

	betas = t(betas)
	rownames(betas) <- series_matrix$title

	betas <- remove_cell_types(betas, remove_names)
	save_data(betas, 'betas_selected_cells')

	bi <- generate_bimodal_distribution(betas)
	save_data(bi, 'bi')

	betas <- fit_betas_to_distribution(betas, bi)
	save_data(betas, 'betas_fit_to_distribution')

	consensus_vector <- create_consensus_vector(betas, group_names)
	names(consensus_vector) <- c('NK.cells', 'Monocytes', 'Tc.cells', 'Neutrophils','Th.cells', 'B.cells')
	save_data(consensus_vector, 'consensus_vector')

	consensus_vector <- remove_unapplicable_probes(consensus_vector)

	return(consensus_vector)
}

#' Removes NA probes
#'
#' Returns consensus vector with NA probes removed
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return consensus vector without NA probes
#' @examples
#' consensus_vector <- remove_NA_probes(consensus_vector)
#' @export
remove_NA_probes <- function(consensus_vector) {
	NA_remove <- rownames(consensus_vector)[is.na(rowSums(consensus_vector))]
	print(sprintf("Number of NAs: %d", length(NA_remove)))
	
	consensus_vector = consensus_vector[!(rownames(consensus_vector) %in% NA_remove), ]
	return(consensus_vector)
}

#' Removes intermediate probes
#'
#' Returns consensus vector with intermediate probes removed.
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return consensus vector without intermediate probes
#' @examples
#' consensus_vector <- remove_intermediate_probes(consensus_vector)
#' ... some visualization 
#' @export
remove_intermediate_probes <- function(consensus_vector) {
	nprobes = nrow(consensus_vector)

	intermediate_remove = unique(which(consensus_vector == -1) %% nprobes)
	intermediate_remove[intermediate_remove == 0] = nprobes	
	intermediate_remove <- rownames(consensus_vector)[intermediate_remove]
	print(sprintf("Number of intermediates: %d", length(intermediate_remove)))

	consensus_vector <- consensus_vector[!(rownames(consensus_vector) %in% intermediate_remove), ]

	return(consensus_vector)
}

#' Removes zero parsimony probes
#'
#' Returns consensus vector with zero parsimony probes removed
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return consensus vector without zero parsimony probes
#' @examples
#' consensus_vector <- remove_zero_parsimony_probes(consensus_vector)
#' ... some visualization 
#' @export
remove_zero_parsimony_probes <- function(consensus_vector) {
	ngroups = ncol(consensus_vector)

	zero_parsimony_remove = rownames(consensus_vector)[rowSums(consensus_vector) == ngroups | rowSums(consensus_vector) == 0]
	print(sprintf("Number of zero parsimony: %d", length(zero_parsimony_remove)))

	consensus_vector <- consensus_vector[!(rownames(consensus_vector) %in% zero_parsimony_remove), ]
	return(consensus_vector)
}

#' Removes all inapplicable probes 
#'
#' Returns consensus vector with NA, intermediate, and zero parsimony probes 
#' removed
#' 
#' @param consensus_vector consensus vector data frame for each cell in group 
#' names
#'
#' @return consensus vector without inapplicable probes
#' @examples
#' consensus_vector <- remove_unapplicable_probes(consensus_vector)
#' @export
remove_unapplicable_probes <- function(consensus_vector) {
	consensus_vector <- remove_NA_probes(consensus_vector)
	consensus_vector <- remove_intermediate_probes(consensus_vector)
	consensus_vector <- remove_zero_parsimony_probes(consensus_vector)

	return(consensus_vector)
}

#' Saves specified data with variable name to working directory 
#'
#' Given data, a variable name, and a location for that data to be stored, 
#' this function simply saves data.
#' 
#' @param data data structure
#' @param variable_name name of variable/file
#' @param dir directory to which the data should be saved
#'
#' @return None
#' @examples
#' save_data(betas, 'betas', getwd())
#' @export
save_data <- function(data, variable_name, dir) {
	if(missing(dir)) {
		dir = getwd()
	}

	data_file = paste(dir, "/", variable_name, '.rds', sep = "")
	saveRDS(data, file = data_file)
}
