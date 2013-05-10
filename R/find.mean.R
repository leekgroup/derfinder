## just used in case you don't want to use the above find.mean.up or find.mean.down.  arguments/return are the same.


#'Optional helper function for getParams
#'
#'Finds mean of distribution of either over- or underexpressed statistics.
#'
#'For experienced users/debugging only.  Calls \code{find.mean.up} if up=TRUE,
#'else calls \code{find.mean.down}.  Most users should use \code{getParams}
#'rather than \code{find.mean}.
#'
#'@param init.value number in (0, 0.5) representing a percentile of the null
#'distribution to use as starting value
#'@param null.mean estimated mean of null distribution (usually found with
#'locfdrFit)
#'@param null.sd estimated standard deviation of null distribution (usually
#'found with locfdrFit)
#'@param null.prop estimated proportion of statistics that came from the null
#'distribution
#'@param vals vector of all the observed values from the mixture distribution
#'@param up if TRUE, find mean of overexpressed statistics, otherwise find mean
#'of underexpressed statistics
#'@return If numerical method succeeds, a list with elements
#'\item{m }{estimated mean of the underexpressed distribution, and}
#'\item{p }{the percentile of the null distribution used to find this mean}
#'If numerical method fails, a list with elements
#'\item{m }{mean of underexpressed distribution, as estimated by the 5th percentile of the estimated null distribution}
#'\item{s }{standard deviation of underexpressed distribution, as estimated by the standard deviation of the null distribution}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getParams}}, \code{\link{find.mean.up}},
#'\code{\link{find.mean.down}}

find.mean <- function(init.value, null.mean, null.sd, null.prop, vals, up = TRUE){
	if(up) k = find.mean.up(init.value, null.mean, null.sd, null.prop, vals)
	if(!up) k = find.mean.down(init.value, null.mean, null.sd, null.prop, vals)
	return(k)
}
