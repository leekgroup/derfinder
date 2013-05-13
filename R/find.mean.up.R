## find.mean.up():
## arguments:
## --init.value: intial percentile of null distribution we are going to look at
## --null.mean, null.sd: parameters of null distribution as estimated from locfdr
## --null.prop: proportion of null values as estimated from locfdr
## --vals: full vector of t statistics (or whatever) that you estimated your null/alt distributions from
## return:
## --if our method SUCCEEDS: list with elements $m (the estimated mean of the DE-up distribution) and $p (the percentile of the null distribution used to calculate that mean)
## --if our method FAILS: list with elements $m (estimated mean of the DE-up distribution) and $s (estimated sd of DE-up dist), as estimated with getParams.failsafe().


#'Helper function for getParams
#'
#'Finds mean of the distribution of statistics generated from overexpressed
#'nucleotides
#'
#'This function is for experienced users or debugging only - all other users
#'should use \code{getParams}, which calls this function.
#'
#'@param init.value number in (0.5, 1) representing a percentile of the null
#'distribution to use as starting value
#'@param null.mean estimated mean of null distribution (usually found with
#'locfdrFit)
#'@param null.sd estimated standard deviation of null distribution (usually
#'found with locfdrFit)
#'@param null.prop estimated proportion of statistics that came from the null
#'distribution
#'@param vals vector of all the observed values from the mixture distribution
#'@return If numerical method succeeds, a list with elements
#'
#'If numerical method fails, a list with elements
#'\item{m }{estimated mean of the underexpressed distribution, and}
#'\item{p }{the percentile of the null distribution used to find this mean}
#'If numerical method fails, a list with elements
#'\item{m }{mean of underexpressed distribution, as estimated by the 5th percentile of the estimated null distribution}
#'\item{s }{standard deviation of underexpressed distribution, as estimated by the standard deviation of the null distribution}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getParams}}, \code{\link{find.mean.down}},
#'\code{\link{find.sd}}

find.mean.up <- function(init.value, null.mean, null.sd, null.prop, vals){
	if(init.value<=0.5) stop("finding DE up quantile - init.value should be >0.5")
	x = init.value
	target = round(0.5*0.5*(1-null.prop)*length(vals))  #half the expected # of DEup values.
	history = pseq = numseq = qseq = NULL
	counter = 0
	while(TRUE){
		iter = try(get.numalts(x,null.mean,null.sd,null.prop,vals,up=T),silent=T)
		if(class(iter)=="try-error") {
			warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
			gp = getParams.failsafe(null.mean, null.sd)
			return(list(m=gp$DEup.mean, s=gp$DEup.sd))
		}
		if(iter$num == target) return(list(m=iter$val,p = pseq[counter]))
		if(counter>1){
			if(abs(pseq[counter]-pseq[counter-1])<1e-8) return(list(m=iter$val,p=pseq[counter]))
		}
		counter = counter+1
		numseq[counter] = iter$num
		qseq[counter] = iter$val
		pseq[counter] = x

		if(counter>1){
			if(history[counter-1]==1 & numseq[counter]>=numseq[counter-1]){
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}			
			if(history[counter-1]==-1 & numseq[counter]<=numseq[counter-1]){
				warning("Numerical estimation of DE-up mean failed.  Defaulting to mean = 95th percentile of estimated null distribution, sd = sd of estimated null distribution.")
				gp = getParams.failsafe(null.mean, null.sd)
				return(list(m=gp$DEup.mean, s=gp$DEup.sd))
			}
		}

		if(iter$num > target){
			history[counter] = 1
			if(counter==1){ x = (x+1)/2; next }
		}
		if(iter$num < target){
			history[counter] = -1
			if(counter==1){ x = (x+0.5)/2; next }
		}

		if(history[counter-1] != history[counter]) x = (x+pseq[counter-1])/2
		if(history[counter-1] == history[counter]){
				prev = last(which(history != history[counter]))
				if(length(prev)==0 & iter$num > target){x = (x+1)/2; next}
				if(length(prev)==0 & iter$num < target){x = (x+0.5)/2; next}
				x = (x+pseq[prev])/2
			}
		}
}
