#'Permutation-based significance tests for differentially expressed regions
#'
#'Creates simulated null distribution of average t statistics for regions
#'classified as over- or underexpressed, and obtains p-values for observed
#'over- and underexpressed regions based on this simulated null.
#'
#'
#' @param regions data frame of regions to obtain p-values for: specifically the
#'\code{$states} element of the return of \code{getRegions}.  Must have
#'\code{mean.t} and \code{state} columns.
#' @param num.perms Number of permutations to use to create the null
#'distribution.
#' @param est.params return from \code{getParams} - HMM parameters used when
#'creating \code{regions}
#' @param chromosome Chromosome you are analyzing.  Currently only runs one
#'chromosome at a time.
#' @param verbose If \code{TRUE} statements with status updates are printed at each permutation.
#' @param DF A list containing a \code{DF} DataFrame with the coverage data and a logical Rle with the positions that passed the cutoff. This object is generated using \code{\link{makeDF}}.
#' @param comparison Either \code{twogroup}, \code{multigroup} or \code{expression}. \code{multigroup} will use the F-statistic and \code{expression} tests the intercept-only model.
#' @param group 0/1 vector denoting group labels for the samples used in
#'analysis. Should have the same length as \code{colsubset} if \code{colsubset} is specified.
#' @param chunksize How many rows of the merged table should be processed at a time?
#' @param adjustvars Optional matrix of adjustment variables (e.g. measured
#'confounders, output from SVA, etc.) to use in fitting linear models to each
#'nucleotide.
#' @param colsubset Optional vector of column indices of the input file that
#'denote samples you wish to include in analysis. 
#' @param scalefac A log transformation is used on the count tables, so zero
#'counts present a problem.  What number should we add to the entire matrix
#'before running the models?  Defaults to 32.
#' @param nonzero If TRUE, use the median of only the nonzero counts as the
#'library size adjustment.
#'
#'@return A vector having length equal to the number of rows in \code{regions},
#'giving a p-value for each region of state 3 or 4 in \code{regions}.
#'@author Alyssa Frazee, Leonardo Collado-Torres
#'@export
#'@import IRanges
#'@seealso
#'\code{\link{getLimmaInput.DF}},\code{\link{getTstats}},\code{\link{getRegions}}

get.pvals.DF <- function(regions, num.perms = 1, est.params, chromosome, verbose=FALSE, DF, comparison = c("twogroup", "multigroup", "expression"), group, chunksize=1e+05, adjustvars=NULL, colsubset=NULL, scalefac=32, nonzero=FALSE){
	# ... should indicate other arguments needed for:
	# getLimmaInput, getTstats, getRegions
	# FIX COLSUBSET + other non-required arguments...
	
	## Required pkgs
	require("IRanges")
	
	## Initialize the data to be used
	data <- DF
	
	## Check that the lengths are the same for group and colsubset (if specified)
	if(!is.null(colsubset)) {
		if(length(group) != length(colsubset)) {
			stop("'group' and 'colsubset' must have the same length when 'colsubset' is specified")
		}
		
		## Subset from the beginning, before the permutations. Thus removing the need of using colsubset later on.
		data$DF <- data$DF[, colsubset]
	}
	
	## Initialize vector for the null statistics
	nullstats <- NULL
	
	## Main permutation loop
	for(i in 1:num.perms){
		idx.permute <- sample(1:length(group))
		group.permute <- group[idx.permute]
		
		## Correctly use the adjusted vars according to the permutation if they were specified
		if(!is.null(adjustvars)) {
			adjustvars.permute <- adjustvars[idx.permute, ]
		} else {
			adjustvars.permute <- NULL
		}
		
		
		if(verbose) print(paste("Getting limma input for permutation", i))
		limma.input <- getLimmaInput.DF(DF=DF, comparison=comparison, group=group.permute, chunksize=chunksize, adjustvars=adjustvars.permute, colsubset=NULL, nonzero=nonzero, scalefac=scalefac)
		
		if(verbose) print(paste("Finding t statistics for permutation", i))
		tstats <- getTstats(fit = limma.input$ebobject, trend = TRUE)
		tt <- tstats$tt
		logfchange <- tstats$logfchange
		
		if(verbose) print(paste("Getting regions for permutation", i))
		regions.null <- getRegions(method = "HMM", chromosome = chromosome, pos = limma.input$pos, tstats = tt, stateprobs = est.params$stateprobs, params = est.params$params, includet = TRUE, includefchange = TRUE, fchange = logfchange)
		nullstats <- append(nullstats, regions.null$states$mean.t[regions.null$states$state==3 | regions.null$states$state==4])
	}
	
	## Construct the pvalues vector
	if(verbose) print("Constructing p-values result vector")
	pvals <- rep(NA, dim(regions)[1])
	for(k in which(regions$state==3 | regions$state==4)){
		pvals[k] <- (sum(abs(nullstats) > abs(regions$mean.t[k])) + 1) / (length(nullstats) + 1)
	}

	## Done!
	return(pvals)
}
