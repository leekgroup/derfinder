\name{find.mean.down}
\alias{find.mean.down}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ helper function for getParams
}
\description{
finds mean of the distribution of statistics generated from underexpressed nucleotides
}
\usage{
find.mean.down(init.value, null.mean, null.sd, null.prop, vals)
}
\arguments{
  \item{init.value}{
number in (0, 0.5) representing a percentile of the null distribution to use as starting value
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getParams}}, \code{\link{find.mean.up}}, \code{\link{find.sd}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
