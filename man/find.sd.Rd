\name{find.sd}
\alias{find.sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ helper function for getParams
}
\description{
find standard deviation of distribution of over- or underexpressed statistics
}
\usage{
find.sd(prev.p, found.mean, null.mean, null.sd, null.prop, vals, up = T)
}
\arguments{
  \item{prev.p}{
percentile of null distribution used to find the mean of the distribution of interest - usually the \code{$p} return of \code{find.mean.up} or \code{find.mean.down}.
}
  \item{found.mean}{
mean of distribution of interest - usually the \code{$m} return of \code{find.mean.up} or \code{find.mean.down}
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
if TRUE, find standard deviation of overexpressed statistics, otherwise find standard deviation of underexpressed statistics
}
}
\details{
This function is for experienced users or debugging only - all other users should use \code{getParams}, which calls this function. 
}
\value{
the estimated standard deviation of the over- or underexpressed statistics
}

\author{
Alyssa Frazee
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getParams}}, \code{\link{find.mean.up}}, \code{\link{find.mean.down}}
}
