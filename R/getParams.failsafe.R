## getParams.failsafe:
## arguments: null.mean and null.sd, mean and sd of estimated null distribution
## return: list with $DEup.mean, $DEdown.mean, $DEup.sd, and $DEdown.sd, giving parameters of alternative distributions as estimated by locfdr
## used in case our numerical estimation method fails.
## just use the 5th (for DE down) and 95th (for DE up) percentiles of the null distribution as alternative means
## and then have the alternative SDs be equal to the null SD.


#'Helper function for getParams
#'
#'When numerical methods \code{find.mean} and \code{find.sd} fail,
#'\code{getParams.failsafe} is used to calculate parameters of the
#'distributions of t statistics originating from over- or underexpressed
#'nucleotides.
#'
#'For experienced users/debugging only. Most users should use \code{getParams}
#'directly.
#'
#'@param null.mean Estimated mean of null distribution (usually from
#'\code{locfdrFit})
#'@param null.sd Estimated standard deviation of null distribution (usually
#'from \code{locfdrFit})
#'@return A list with elements
#'\item{DEup.mean }{estimated mean of overexpressed distribution, calculated as the 95th percentile of the estimated null distribution}
#'\item{DEup.sd }{estimated standard deviation of overexpressed distribution, set to be equal to the estimated standard deviation of the null distribution}
#'\item{DEdown.mean }{estimated mean of underexpressed distribution, calculated as the 5th percentile of the estimated null distribution}
#'\item{DEdown.sd }{estimated standard deviation of underexpressed distribution, set to be equal to the estimated standard deviation of the null distribution}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getParams}}

getParams.failsafe <- function(null.mean, null.sd){
	DEup.mean = qnorm(.95, mean=null.mean, sd=null.sd)
	DEup.sd = DEdown.sd = null.sd
	DEdown.mean = qnorm(.05, mean=null.mean, sd=null.sd)
	return(list(DEup.mean = DEup.mean, DEdown.mean = DEdown.mean, DEup.sd = DEup.sd, DEdown.sd = DEdown.sd))
}
