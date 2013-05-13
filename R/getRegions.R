## getRegions()
## arguments:
## --method: can be "HMM" for hidden markov model, "smoothcut" for Rafa's method from last summer (smooth t-stats and choose a cutoff), or "CBS" for circular binary segmentation
## --chromosome: name of the chromosome being analyzed. Only used for output purposes (will be printed in returned data frame)
## --pos: vector telling which bases are being tested.  Returned by getLimmaInput(), and is first column of .db file specified earlier in pipeline (e.g., in makeDb())
## --tstats: vector giving t statistics for each base. Same length as pos. 
## --transprobs: used in "HMM" method - first element determines the diagonal of the transition matrix
## --transprobs (cont): second element determines likelihood of transitioning from expressed to DE or vice versa, or from DE up to DE down or vice versa - should be very small
## --transprobs (cont): the remainder of the transition matrix is then fixed. 
## --stateprobs: used in "HMM" method - length 4 vector giving probability of being in each of 0, null, DE up, or DE down states.  Usually obtained with getParams()
## --params: used in "HMM" method - list containing the length-4 vectors $mean and $sd, giving means and sds of normal distributions for 0, null, DE up, or DE down states.  Usually obtained with getParams()
## --K: used in the "smoothcut" method - width of window used in smoothing t-statistics
## --tcut: used in the "smoothcut" method - t statistic cutoff above which a base will be classified as differentially expressed
## --includet: if TRUE, output includes t-statistics for each base in addition to state calls
## --includefchange: if TRUE, output includes the fold change (not log scale) in expression (between groups) for each base
## --fchange: required if includefchange=TRUE - log2 fold changes for each base. These are returned by getTstats().
## return:
## a list, one element is $states, a data frame containing columns chr, start, end, and state, one row for each region identified.
## state=1 means not expressed, 2 means expressed (no DE), 3 means DE up, 4 means DE down
## if includet and/or includefchange are true, $states will contain columns mean.t and/or mean.fold.change, giving the mean t statistic and/or mean fold change for the region.
## return list also contains $states.norle, which is a 3-column data frame giving chr, pos, and state for each position. Useful in plotting.





#'Generate list of regions, classify each as differentially expressed or not
#'
#'Using one of three methods, divides the genome (or chromosome) into regions
#'by putting each nucleotide into a state and grouping contiguous nucleotides
#'of the same state into "regions."  Regions of states 3 and 4 are
#'"differentially expressed."
#'
#'States are labeled numerically in the output as follows: 1="not expressed,"
#'2="equally expressed," 3="overexpressed," 4="underexpressed."
#'
#'@param method Can be one of "HMM" (Hidden Markov Model), "CBS" (circular
#'binary segmentation), or "smoothcut" (t statistics with high enough absolute
#'values are called differentially expressed).
#'@param chromosome Name of chromosome being analyzed - will be printed in
#'output table.
#'@param pos Vector giving genomic positions of the provided t statistics. Must
#'have length equal to that of \code{tstats}.  \code{pos} is returned by
#'\code{getLimmaInput}.
#'@param tstats Vector giving moderated t statistics, in proper genomic order.
#'@param transprobs Vector denoting transition probabilities between states,
#'for use in the "HMM" method. Should have length 2, with first element
#'denoting the probability of staying in the same state (should be large), and
#'the second element denoting the probability of moving directly from a
#'differentially expressed state to an equally expressed state or vice versa,
#'or from an overexpressed state to an underexpressed state or vice versa
#'(should be very small).  Defaults to c(.999, 1e-12).
#'@param stateprobs Marginal probabilities of being in each of the four hidden
#'states, for use with the "HMM" method.  The \code{stateprobs} element of
#'\code{getParams} generates this.
#'@param params Parameters of the normal distributions representing the four
#'states in the "HMM" method. The \code{params} element of \code{getParams}
#'generates this.
#'@param K Smoothing parameter used in the "smoothcut" method: t statistics are
#'smoothed using running median; how wide should the window be?  Default 25.
#'@param tcut Cutoff used in the "smoothcut" method to classify differential
#'expression: how large in absolute value should a moderated t statistic be in
#'order to be classified as having been generated from a differentially
#'expressed nucleotide?  Default 2.
#'@param includet If TRUE, the table in the output will include the average t
#'statistic for each region.
#'@param includefchange If TRUE, the table in the output will include the
#'average estimated fold change (as estimated from the linear models) for each
#'region.
#'@param fchange Required if \code{includefchange = TRUE}. Estimated log2 fold
#'changes from the linear models - should have length equal to that of
#'\code{tstats}. Usually obtained from the \code{logfchange} element of the
#'output of \code{getTstats}.
#'@return A list with elements
#'\item{states.norle }{data frame with one row per nucleotide, giving its genomic location and predicted hidden state}
#'\item{states }{data frame with one row per region, giving its genomic location, length, predicted hidden state, and (if applicable) average t statistic and/or fold change.}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getTstats}}, \code{\link{getParams}}

getRegions <- function(method, chromosome, pos, tstats, transprobs = c(0.999, 1e-12), stateprobs = NULL, params = NULL, K = 25, tcut = 2, includet=F, includefchange=F, fchange=NULL) {
  	if(method!="HMM" & method!="smoothcut" & method!="CBS") stop("Invalid method. Choices are HMM, smoothcut, or CBS")
  	if(sum(sort(pos)!=pos)>0) stop("pos (and probably t-statistics) improperly sorted")
    if (method == "HMM"){
		require(HiddenMarkov)
		if(is.null(params)) stop("HMM method requires parameter estimation - please provide theoretical values or estimate with getParams")
		stayprob = transprobs[1]
    	EtoDE = transprobs[2]
    	EtoZero = 1-stayprob-2*EtoDE
    	transmat = matrix(c(stayprob,(1-stayprob)/3,(1-stayprob)/3,(1-stayprob)/3,				 							         							EtoZero,stayprob,EtoDE,EtoDE,
					     	EtoZero,EtoDE,stayprob,EtoDE,
					     	EtoZero,EtoDE,EtoDE,stayprob),nrow=4,ncol=4,byrow=T)
    	if(is.null(stateprobs)) stop("HMM method requires initial values for each state - please provide theoretical values or estimate with getParams")
    	if(sum(names(params)==c("mean","sd")) < length(params)) stop("params argument formatted incorrectly; please see getParams")
    
   		hmmodel = dthmm(x=tstats, Pi=transmat, delta=stateprobs, distn="norm",pm=params) 
    	state.predictions = Viterbi(hmmodel)
		} #end HMM
		
	if (method == "smoothcut"){	
		smootht = runmed(tstats,k=K,endrule="constant")
		state.predictions = rep(NA,length(smootht))
		state.predictions[which(smootht==0)] = 1
		state.predictions[which(smootht>0 & smootht<tcut)] = 2
		state.predictions[which(smootht<0 & smootht>-tcut)] = 2
		state.predictions[which(smootht>=tcut)] = 3
		state.predictions[which(smootht<=-tcut)] = 4
		} # end smoothcut
		
	if(method == "CBS"){	
		if(is.null(params)) stop("CBS method requires parameter estimation - please provide theoretical values or estimate with getParams")
		require(DNAcopy)
		cna.object = CNA(tstats,chrom=rep(chromosome,length(tstats)),maploc=pos,data.type="logratio",sampleid="all",presorted=TRUE)
		cna.smoothed = smooth.CNA(cna.object)
		segmentation = segment(cna.smoothed,p.method="hybrid",verbose=0) # slow.

		statecalls = segmentation$output$seg.mean

		zprobs = 1-pnorm(abs((statecalls-params$mean[1])/params$sd[1]))
    	nullprobs = 1-pnorm(abs((statecalls-params$mean[2])/params$sd[2]))
		upprobs = 1-pnorm(abs((statecalls-params$mean[3])/params$sd[3]))
		downprobs = 1-pnorm(abs((statecalls-params$mean[4])/params$sd[4]))
		state.predictions.segments = max.col(cbind(zprobs,nullprobs,upprobs,downprobs))
		statetable = list(lengths=segmentation$output$num.mark, values=state.predictions.segments)
		state.predictions = inverse.rle(statetable)
		} # end CBS
		
    states.norle = data.frame(chr=rep(chromosome,length(pos)),pos,state=state.predictions)
    
    pos.full = c(pos[1]:pos[length(pos)])
    state.preds.full = rep(1,length(pos.full)) 
    state.preds.full[pos-pos[1]+1] = state.predictions
    full.rle = rle(state.preds.full)
    END = pos.full[cumsum(full.rle$lengths)]
    START = END-full.rle$lengths+1
    states = data.frame(chr=rep(chromosome,length(START)),start=START,end=END,state=full.rle$values,length=full.rle$lengths)
    if(includet & !includefchange){
    	full.t = rep(0,max(pos))
    	full.t[pos] = tstats
    	tmean = NULL
    	for(i in 1:dim(states)[1]){
    		tmean[i] <- mean(full.t[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.t=tmean)
    }
    if(!includet & includefchange){
    	if(is.null(fchange)) stop("must provide fold changes (log2 scale) if you want to include them in your output (output fold change is not on log scale)")
    	full.fchange = rep(1,max(pos))
    	full.fchange[pos] = 2^fchange
    	fcmean = NULL
    	for(i in 1:dim(states)[1]){
    		fcmean[i] <- mean(full.fchange[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.fold.change=fcmean)
    }
	if(includet & includefchange){
    	full.t = rep(0,max(pos))
    	full.t[pos] = tstats
    	full.fchange = rep(1,max(pos))
    	full.fchange[pos] = 2^fchange
    	tmean = fcmean = NULL
    	for(i in 1:dim(states)[1]){
       		fcmean[i] <- mean(full.fchange[states$start[i]:states$end[i]])
    		tmean[i] <- mean(full.t[states$start[i]:states$end[i]])
    	}
    	states = data.frame(states,mean.t=tmean,mean.fold.change=fcmean)
	}
    out = list(states.norle=states.norle, states=states)
    return(out)

}
