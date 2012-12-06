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

## return:  a vector of length = nrow(regions), giving a p-value for each region in regions that is of state 3 or 4 (DE).  (regions of state 1 or 2 are assigned NA)


get.pvals = function(regions, dbfile, tablename, num.perms = 1, group, est.params, chromosome, colsubset = c(-1)){
	# ... should indicate other arguments needed for:
	# getLimmaInput, getTstats, getRegions
	# FIX COLSUBSET + other non-required arguments...
	
	nullstats = NULL
	
	for(i in 1:num.perms){
		group.permute = sample(group)
		print("getting limma input...")
		limma.input = getLimmaInput(dbfile = dbfile, tablename = tablename, group = group.permute, colsubset = colsubset)
		print("finding t statistics...")
		tstats = getTstats(fit = limma.input$ebobject, trend = TRUE)
		tt = tstats$tt
		logfchange = tstats$logfchange
		print("getting regions...")
		regions.null = getRegions(method = "HMM", chromosome = chromosome, pos = limma.input$pos, tstats = tt, stateprobs = est.params$stateprobs, params = est.params$params, includet = TRUE, includefchange = TRUE, fchange = logfchange)
		nullstats = append(nullstats, regions.null$states$mean.t[regions.null$states$state==3|regions.null$states$state==4])
	}
	
	pvals = rep(NA,dim(regions)[1])
	for(k in which(regions$state==3|regions$state==4)){
		pvals[k] = sum(abs(nullstats)>abs(regions$mean.t[k]))/length(nullstats)
	}

	return(pvals)
}
