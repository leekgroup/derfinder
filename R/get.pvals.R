## calculate permutation p-values for each region classified by the HMM as DE
## arguments:
## regions:  data frame of genomic regions (usually the $states output of getRegions).  must contain "state" and "mean.t".
## dbfile:  name of database containing counts (the one created with makeDb)
## tablename: name of table in dbfile
## num.perms:  how many times should the permutation procedure be repeated?  default 1, but recommendation is to use more (maybe 10?)
## group:  sample labels (0/1)
## est.params:  the return of getParams() on the original t-statistics
## chromosome:  which chromosome are you analyzing? 
## colsubset: optional, used if you aren't using all the samples in dbfile
## adjustvars: if desired, an nxp matrix of adjustment variables (confounders, SVA output, etc.)
## nonzero:  if TRUE, library size adjustment is median of nonzero counts.  If FALSE, library size adjustment is median of all counts.
## scalefac: how much should counts in the table be offset by before taking the log (to avoid taking log(0))? Default 32.

## colmeds: vector with the column medians or NULL if you want to calculate them this time around. To save computing time, you might want to calculate them once with getColmeds().

## return:  a vector of length = nrow(regions), giving a p-value for each region in regions that is of state 3 or 4 (DE).  (regions of state 1 or 2 are assigned NA)




#'Permutation-based significance tests for differentially expressed regions
#'
#'Creates simulated null distribution of average t statistics for regions
#'classified as over- or underexpressed, and obtains p-values for observed
#'over- and underexpressed regions based on this simulated null.
#'
#'
#'@param regions data frame of regions to obtain p-values for: specifically the
#'\code{$states} element of the return of \code{getRegions}.  Must have
#'\code{mean.t} and \code{state} columns.
#'@param dbfile Database file used to obtain the \code{regions} data frame.
#'@param tablename Name of the table in \code{dbfile}.
#'@param num.perms Number of permutations to use to create the null
#'distribution.
#'@param group 0/1 vector denoting group labels for the samples used in
#'analysis
#'@param est.params return from \code{getParams} - HMM parameters used when
#'creating \code{regions}
#'@param chromosome Chromosome you are analyzing.  Currently only runs one
#'chromosome at a time.
#'@param colsubset Colsubset argument used to create \code{regions}
#'@param adjustvars Optional matrix of adjustment variables (e.g. measured
#'confounders, output from SVA, etc.) to use in fitting linear models to each
#'nucleotide.
#'@param nonzero If TRUE, use the median of only the nonzero counts as the
#'library size adjustment.
#'@param scalefac A log transformation is used on the count tables, so zero
#'counts present a problem.  What number should we add to the entire matrix
#'before running the models?  Defaults to 32.
#'@param chunksize How many rows of the database should be processed at a time?
#'@param colmeds If NULL, the column medians are calculated using
#'\code{\link{getColmeds}}. Otherwise, the output of \code{\link{getColmeds}}
#'is expected.
#'@return A vector having length equal to the number of rows in \code{regions},
#'giving a p-value for each region of state 3 or 4 in \code{regions}.
#'@author Alyssa Frazee
#'@export
#'@seealso
#'\code{\link{getLimmaInput}},\code{\link{getTstats}},\code{\link{getRegions}}

get.pvals = function(regions, dbfile, tablename, num.perms = 1, group, est.params, chromosome, colsubset = c(-1), adjustvars=NULL, nonzero=FALSE, scalefac=32, chunksize=1e+05, colmeds=NULL){
	# ... should indicate other arguments needed for:
	# getLimmaInput, getTstats, getRegions
	# FIX COLSUBSET + other non-required arguments...
	
	nullstats = NULL
	
	for(i in 1:num.perms){
		idx.permute <- sample(1:length(group))
		group.permute <- group[idx.permute]
		if(!is.null(colmeds)) {
			colmeds.permute <- colmeds[idx.permute]
		} else {
			colmeds.permute <- NULL
		}
		
		## Correctly use the adjusted vars and colsubset according ot the permutation if they were specified
		if(!is.null(adjustvars)) {
			adjustvars.permute <- adjustvars[idx.permute, ]
		} else {
			adjustvars.permute <- NULL
		}
		if(colsubset[1]==c(-1) & length(colsubset) == 1) {
			## This is the case where the user wants to ignore the pos data (1st column), but we need to pass the permutation info to getLimmaInput()
			colsubset.permute <- (1:length(group) + 1)[idx.permute]
		} else {
			colsubset.permute <- colsubset[idx.permute]	
		}
		
		
		#print("getting limma input...")
		limma.input = getLimmaInput(dbfile = dbfile, tablename = tablename, group = group.permute, colsubset=colsubset.permute, adjustvars=adjustvars.permute, nonzero=nonzero, scalefac=scalefac, chunksize=chunksize, colmeds=colmeds.permute)
		#print("finding t statistics...")
		tstats = getTstats(fit = limma.input$ebobject, trend = TRUE)
		tt = tstats$tt
		logfchange = tstats$logfchange
		#print("getting regions...")
		regions.null = getRegions(method = "HMM", chromosome = chromosome, pos = limma.input$pos, tstats = tt, stateprobs = est.params$stateprobs, params = est.params$params, includet = TRUE, includefchange = TRUE, fchange = logfchange)
		nullstats = append(nullstats, regions.null$states$mean.t[regions.null$states$state==3|regions.null$states$state==4])
	}
	
	pvals = rep(NA,dim(regions)[1])
	for(k in which(regions$state==3|regions$state==4)){
		pvals[k] = (sum(abs(nullstats)>abs(regions$mean.t[k]))+1)/length(nullstats)
	}

	return(pvals)
}
