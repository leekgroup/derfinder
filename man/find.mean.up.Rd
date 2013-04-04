\name{find.mean.up}
\alias{find.mean.up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ helper function for getParams
}
\description{
finds mean of the distribution of statistics generated from overexpressed nucleotides
}
\usage{
find.mean.up(init.value, null.mean, null.sd, null.prop, vals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init.value}{
number in (0.5, 1) representing a percentile of the null distribution to use as starting value
}
  \item{null.mean}{
estimated mean of null distribution (usually found with locfdrFit)
}
  \item{null.sd}{
estimated standard deviation of null distribution (usually found with locfdrFit)
}
  \item{null.prop}{
estimated proportion of statistics that came from the null distribution
}
  \item{vals}{
vector of all the observed values from the mixture distribution
}
}
\details{
This function is for experienced users or debugging only - all other users should use \code{getParams}, which calls this function. 
}
\value{
If numerical method succeeds, a list with elements 
\item{m }{estimated mean of the underexpressed distribution, and}
\item{p }{the percentile of the null distribution used to find this mean}
If numerical method fails, a list with elements
\item{m }{mean of underexpressed distribution, as estimated by the 5th percentile of the estimated null distribution}
\item{s }{standard deviation of underexpressed distribution, as estimated by the standard deviation of the null distribution}
}
\author{
Alyssa Frazee
}


\seealso{
\code{\link{getParams}}, \code{\link{find.mean.up}}, \code{\link{find.sd}}
}
