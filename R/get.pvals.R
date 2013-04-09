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
