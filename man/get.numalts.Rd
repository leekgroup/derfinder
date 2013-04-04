\name{get.numalts}
\alias{get.numalts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ helper function for getParams
}
\description{
find estimated number of alternative statistics above a set percentile of the estimated null distribution
}
\usage{
get.numalts(pctil, null.mean, null.sd, null.prop, vals, up = TRUE)
}
\arguments{
  \item{pctil}{
percentile, in (0,1), of the null distribution for which the number of alternative statistics above (if \code{pctil} is greater than 0.5) or below (if \code{pctil} is less than 0.5) is desired.
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
if TRUE, get the number of overexpressed statistics above the 100\code{pctil}-th percentile of the null distribution, else get the number of underexpressed statistics below the 100\code{pctil}-th percentile of the null distribution
}
}
\details{
This function is for experienced users or debugging only - all other users should use \code{getParams}, which calls this function. }
\value{
a list with elements
\item{num }{the estimated number of alternative values above/below \code{val} (see \code{val} below)}
\item{val }{the 100\code{pctil}-th percentile of the null distribution}
}
\author{
Alyssa Frazee
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getParams}}
}
