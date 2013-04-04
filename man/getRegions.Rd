\name{getRegions}
\alias{getRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ generate list of regions, classify each as differentially expressed or not
}
\description{
Using one of three methods, divides the genome (or chromosome) into regions by putting each nucleotide into a state and grouping contiguous nucleotides of the same state into "regions."  Regions of states 3 and 4 are "differentially expressed." 
}
\usage{
getRegions(method, chromosome, pos, tstats, transprobs = c(0.999, 1e-12), stateprobs = NULL, params = NULL, K = 25, tcut = 2, includet = F, includefchange = F, fchange = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
Can be one of "HMM" (Hidden Markov Model), "CBS" (circular binary segmentation), or "smoothcut" (t statistics with high enough absolute values are called differentially expressed).}
\item{chromosome}{
Name of chromosome being analyzed - will be printed in output table.
}
  \item{pos}{
Vector giving genomic positions of the provided t statistics. Must have length equal to that of \code{tstats}.  \code{pos} is returned by \code{getLimmaInput}.
}
  \item{tstats}{
Vector giving moderated t statistics, in proper genomic order. 
}
  \item{transprobs}{
Vector denoting transition probabilities between states, for use in the "HMM" method. Should have length 2, with first element denoting the probability of staying in the same state (should be large), and the second element denoting the probability of moving directly from a differentially expressed state to an equally expressed state or vice versa, or from an overexpressed state to an underexpressed state or vice versa (should be very small).  Defaults to c(.999, 1e-12).  
}
  \item{stateprobs}{
Marginal probabilities of being in each of the four hidden states, for use with the "HMM" method.  The \code{stateprobs} element of \code{getParams} generates this.
}
  \item{params}{
Parameters of the normal distributions representing the four states in the "HMM" method. The \code{params} element of \code{getParams} generates this.
}
  \item{K}{
Smoothing parameter used in the "smoothcut" method: t statistics are smoothed using running median; how wide should the window be?  Default 25.
}
  \item{tcut}{
Cutoff used in the "smoothcut" method to classify differential expression: how large in absolute value should a moderated t statistic be in order to be classified as having been generated from a differentially expressed nucleotide?  Default 2. 
}
  \item{includet}{
If TRUE, the table in the output will include the average t statistic for each region.
}
  \item{includefchange}{
If TRUE, the table in the output will include the average estimated fold change (as estimated from the linear models) for each region.
}
  \item{fchange}{
Required if \code{includefchange = TRUE}. Estimated log2 fold changes from the linear models - should have length equal to that of \code{tstats}. Usually obtained from the \code{logfchange} element of the output of \code{getTstats}.
}
}
\details{
States are labeled numerically in the output as follows: 1="not expressed," 2="equally expressed," 3="overexpressed," 4="underexpressed."
}
\value{
a list with elements
\item{states.norle }{data frame with one row per nucleotide, giving its genomic location and predicted hidden state}
\item{states }{data frame with one row per region, giving its genomic location, length, predicted hidden state, and (if applicable) average t statistic and/or fold change.}
}
\author{
Alyssa Frazee
}

\seealso{
\code{\link{getTstats}}, \code{\link{getParams}}
}
