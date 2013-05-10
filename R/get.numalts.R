## get.numalts():
## arguments:
## --pctil: a percentile (e.g. 0.01 for first percentile, 0.88 for 88th percentile...)
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## --up: TRUE if you're estimating the DEup parameters, FALSE if estimating DEdown parameters.
## return: the estimated number of alternative values (in vals) above (below, if up=F) the pctil'th percentile of the null distribution
## the above number is in element $num, the element $val contains the pctil'th percentile of the null distribution.
## note that this function probably shouldn't be used directly: only used by find.mean and find.sd (below)


#'Helper function for getParams
#'
#'Find estimated number of alternative statistics above a set percentile of the
#'estimated null distribution
#'
#'This function is for experienced users or debugging only - all other users
#'should use \code{getParams}, which calls this function.
#'
#'@param pctil percentile, in (0,1), of the null distribution for which the
#'number of alternative statistics above (if \code{pctil} is greater than 0.5)
#'or below (if \code{pctil} is less than 0.5) is desired.
#'@param null.mean estimated mean of null distribution (usually found with
#'locfdrFit)
#'@param null.sd estimated standard deviation of null distribution (usually
#'found with locfdrFit)
#'@param null.prop estimated proportion of statistics that came from the null
#'distribution
#'@param vals vector of all the observed values from the mixture distribution
#'@param up if TRUE, get the number of overexpressed statistics above the
#'100\code{pctil}-th percentile of the null distribution, else get the number
#'of underexpressed statistics below the 100\code{pctil}-th percentile of the
#'null distribution
#'@return A list with elements
#'\item{num }{the estimated number of alternative values above/below \code{val} (see \code{val} below)}
#'\item{val }{the 100\code{pctil}-th percentile of the null distribution}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getParams}}

get.numalts <- function(pctil,null.mean,null.sd,null.prop,vals,up = TRUE){
	cutoff.val = qnorm(pctil,null.mean,null.sd)
	num.nulls = null.prop*length(vals)
	# I use "above" to mean "above or below". very loose...
	coef1 = ifelse(up,1-pctil,pctil)
	num.nulls.above = coef1*num.nulls
	num.total.above = ifelse(up,sum(vals>cutoff.val),sum(vals<cutoff.val))
	num.alts.above = num.total.above - num.nulls.above
	if(num.alts.above<=0) stop("negative num.alts.above. rethink algorithm. sad day.") # LOL at the stop message (Leo)
	return(list(num=round(num.alts.above),val=cutoff.val))
}
