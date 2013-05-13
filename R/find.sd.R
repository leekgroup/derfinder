## find.sd(): used to find SD of alternative distribution if find.mean.up or down succeeds
## arguments:
## --prev.p: percentile used to find the alternative mean (i.e., the $p return from find.mean.up or down)
## --found.mean: alternative mean found (i.e., the $m return from find.mean.up or down)
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## --up: TRUE if you want the DEup sd, FALSE if you want the DEdown.sd
## return:
## --if this method succeeds, the estimated standard deviation of the alternative distribution (using this method)
## --if this method succeeds, the estimated sd of the alternative distribution using the output from locfdr. warning message is printed.


#'Helper function for getParams
#'
#'Find standard deviation of distribution of over- or underexpressed statistics
#'
#'This function is for experienced users or debugging only - all other users
#'should use \code{getParams}, which calls this function.
#'
#'@param prev.p percentile of null distribution used to find the mean of the
#'distribution of interest - usually the \code{$p} return of
#'\code{find.mean.up} or \code{find.mean.down}.
#'@param found.mean mean of distribution of interest - usually the \code{$m}
#'return of \code{find.mean.up} or \code{find.mean.down}
#'@param null.mean estimated mean of null distribution (usually found with
#'locfdrFit)
#'@param null.sd estimated standard deviation of null distribution (usually
#'found with locfdrFit)
#'@param null.prop estimated proportion of statistics that came from the null
#'distribution
#'@param vals vector of all the observed values from the mixture distribution
#'@param up if TRUE, find standard deviation of overexpressed statistics,
#'otherwise find standard deviation of underexpressed statistics
#'@return The estimated standard deviation of the over- or underexpressed
#'statistics
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getParams}}, \code{\link{find.mean.up}},
#'\code{\link{find.mean.down}}

find.sd <- function(prev.p, found.mean, null.mean, null.sd, null.prop, vals, up=T){
	if(up) new.percentile = (1+prev.p)/2 # we increase the percentile we're going to look at
	if(!up) new.percentile = prev.p/2 #decrease it.
	cutoff.val = qnorm(new.percentile, null.mean, null.sd) # x-value that is that percentile of the theoretical null distribution
	num.nulls = null.prop*length(vals)
	coef1 = ifelse(up,1-new.percentile,new.percentile)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	num.alts = (1-null.prop)*0.5*length(vals) # this is the number of DEup values
	if(up) alt.percentile = 1-(num.alts.above/num.alts)
	if(!up) alt.percentile = num.alts.above/num.alts
	zstat = qnorm(alt.percentile)
	if(is.na(zstat)){
		warning("Numerical standard deviation estimation failed.  Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	sigma = (cutoff.val - found.mean)/zstat
    if(num.alts.above<=0 | sigma<=0) {
		warning("Numerical standard deviation estimation failed. Defaulting to sd of estimated null distribution.")
		gp = getParams.failsafe(null.mean, null.sd)
		if(up) return(gp$DEup.sd)
		if(!up) return(gp$DEdown.sd)
	}
	return(sigma)
}
