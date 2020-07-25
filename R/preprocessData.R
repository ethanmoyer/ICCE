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

#' Work on this...
#'
#' Given...
#' 
#' @param x data structure
#' @param mu0 name of variable/file
#' @param mu1 directory to which the data should be saved
#' @param nsamples number of cell types
#'
#' @return None
#' @examples
#' betas_info = generate_p_values(0.2, .119, .880, .112)
#' betas_state = betas_info$betas
#' betas_p = betas_info$betas_p
#' @export
generate_p_value <- function(x, mu0=0.1, mu1=0.9, sigma=0.1) {
	m0_p = pnorm(x, mu0, sigma)
    m1_p = pnorm(x, mu1, sigma)
    if ((1 - m0_p) > m1_p) {
    	return(list(state = 0, p = m0_p))
    } else {
    	return(list(state = 1, p = 1 - m1_p))
    }
}

#' Categorize betas based on bimodal distribution and return the p-value of 
#' each categorization.
#'
#' This function uses a given bimodal distribution generated within cell types
#' in order to categorize betas as either methylated (1) or unmethylated (0) 
#' in the betas structure. It also returns the p-value of each categorization 
#' in the betas_p structure.
#' 
#' @param betas matrix (rows are cell names, columns are probe names)
#' @param bi matrix specifying mu0, mu1, sigma0, sigma1 across all cell types
#'
#' @return betas matrix categorized based on bimodal distribution
#' @examples
#' betas <- fit_betas_to_distribution(betas, bi)
#' @export
fit_betas_to_distribution_p <- function(betas, bi) {
	betas_p = betas
	for(i in 1:nrow(betas)) {
		for(k in 1:ncol(betas)) {
			if (is.na(betas[i, k])) next
			distr = generate_p_value(betas[i, k], bi[i, 1], bi[i ,2], bi[i, 3])
			betas[i, k] = as.integer(distr$state)
			betas_p[i, k] = distr$p
	  	}
	}
	return(list(betas = betas, betas_p = betas_p))
}


create_consensus_vector_p <- function(betas, group_names, mu0=0.1, mu1=0.9, sigma=0.25) {
	consensus_state <- c()
	consensus_p <- c()
	# For each cell type
	for (i in 1:length(group_names)) {
		S_state = c()
		S_p = c()
		name = group_names[i]
		group_betas = betas[grepl(name, rownames(betas)), ]

		# For each probe
		for (j in 1:ncol(group_betas)) {
			# Penalty for the case when there are only NAs at a probe
			if (all(is.na(group_betas[, j]))) {
				S_state[j] = round(runif(1))
				S_p[j] = 1
				next
			}
			
			prod_0 = 1
			# Calculating P(B | m = 0)
			# P(B | m = 0) = (B_1 | m = 0) * ... * (B_N | m = 0)
			for (k in 1:nrow(group_betas)) {
				# Add a penalty based on number of NAs present
				if (is.na(group_betas[k, j])) {
					next
				}

				# P(B_i | m = 0) ~ N(.1, .1)
				prod_0 = prod_0 * (1 - pnorm(group_betas[k, j], mu0, sigma))
			}

			prod_1 = 1
			# Calculating P(B | m = 1)
			# P(B | m = 1) = (B_1 | m = 1) * ... * (B_N | m = 1)
			for (k in 1:nrow(group_betas)) {
				# Add a penalty based on number of NAs present
				if (is.na(group_betas[k, j])) {
					next
				}

				# P(B_i | m = 1) ~ N(.9, .1)
				prod_1 = prod_1 * pnorm(group_betas[k, j], mu1, sigma)
			}
			# Add a penalty for when there is a large discrepancy in beta values
			#if (sd(group_betas[, j]) > 0.15) {}
			# cg00007420
			# arg max (P(B | m = 0), P(B | m = 1))
			prob = c(prod_0, prod_1)
			ind = which.max(prob)
			S_state[j] = c(0, 1)[ind]
			# p = 1 - P(B | m = m) / (P(B | m = 0) + P(B | m = 1))
			S_p[j] = 1 - prob[ind] / (prod_0 + prod_1)

		}
		consensus_state[[name]] = S_state
		consensus_p[[name]] = S_p
	} 

	return(list(consensus_state=data.frame(consensus_state), consensus_p=data.frame(consensus_p)))
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


generate_p_values_0 <- function(x, mu0, mu1, nsamples=10) {
	BL <- vapply(
		c(mu0, mu1),
		function(mu) {
			dbinom(
				round(x*nsamples), 
				size=nsamples, prob=mu)}, numeric(1))

	ind <- which.max(BL)
	BT <- c('1', '0')[ind]
	BS <- floor(-log10(1-BL[ind] / sum(BL)) * 10)
	list(BT=BT, BL=BL[ind], BS=BS)
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
group_names = c("NK", "Mono", "Tc", "Neu", "Th", "B")
remove_names = c("PCA0612", "WB1148", "A2.mix", "B6.mix", "PCA1547", "WB1147", "B3.mix", "B2.mix", "A4.mix", "B4.mix", "A3.mix", "ST2007", "B5.mix", "B1.mix", "A5.mix", "WB1145", "A6.mix", "A1.mix")

betas_to_consensus_vector <- function(reference, series_matrix_file, betas, group_names, remove_names) {

	betas <- remove_sex_chr(reference, betas)
	save_data(betas, 'betas_no_XY')

	series_matrix <- readRDS(series_matrix_file)

	betas = t(betas)
	rownames(betas) <- series_matrix$title

	betas <- remove_cell_types(betas, remove_names)
	save_data(betas, 'betas_selected_cells')

	consensus_info <- create_consensus_vector_p(betas, group_names)
	consensus_state <- consensus_info$consensus_state
	consensus_p <- consensus_info$consensus_p

	#bi <- generate_bimodal_distribution(betas)
	#save_data(bi, 'bi_distribution')

	#betas <- fit_betas_to_distribution(betas, bi)
	#save_data(betas, 'betas_fit_to_distribution')

	#consensus_vector <- create_consensus_vector(betas, group_names)
	names(consensus_state) <- c('NK.cells', 'Monocytes', 'Tc.cells', 'Neutrophils','Th.cells', 'B.cells')
	rownames(consensus_state) = colnames(betas)
	

	names(consensus_p) <- c('NK.cells', 'Monocytes', 'Tc.cells', 'Neutrophils','Th.cells', 'B.cells')
	rownames(consensus_p) = colnames(betas)
	
	consensus_state <- remove_unapplicable_probes(consensus_state)
	consensus_p <- consensus_p[(rownames(consensus_p) %in% rownames(consensus_state)), ]

	save_data(consensus_state, 'consensus_state')
	save_data(consensus_p, 'consensus_p')

	return(list(consensus_state=consensus_state, consensus_p=consensus_p)
}


