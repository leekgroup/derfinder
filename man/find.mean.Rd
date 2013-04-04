\name{find.mean}
\alias{find.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ optional helper function for getParams
}
\description{
finds mean of distribution of either over- or underexpressed statistics.
}
\usage{
find.mean(init.value, null.mean, null.sd, null.prop, vals, up = TRUE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{up}{
if TRUE, find mean of overexpressed statistics, otherwise find mean of underexpressed statistics
}
}
\details{
For experienced users/debugging only.  Calls \code{find.mean.up} if up=TRUE, else calls \code{find.mean.down}.  Most users should use \code{getParams} rather than \code{find.mean}. 
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
\code{\link{getParams}}, \code{\link{find.mean.up}}, \code{\link{find.mean.down}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
