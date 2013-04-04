\name{getParams.failsafe}
\alias{getParams.failsafe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ helper function for getParams
}
\description{
When numerical methods \code{find.mean} and \code{find.sd} fail, \code{getParams.failsafe} is used to calculate parameters of the distributions of t statistics originating from over- or underexpressed nucleotides.
}
\usage{
getParams.failsafe(null.mean, null.sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{null.mean}{
Estimated mean of null distribution (usually from \code{locfdrFit})
}
  \item{null.sd}{
Estimated standard deviation of null distribution (usually from \code{locfdrFit})
}
}
\details{
For experienced users/debugging only. Most users should use \code{getParams} directly.
}
\value{
a list with elements
\item{DEup.mean }{estimated mean of overexpressed distribution, calculated as the 95th percentile of the estimated null distribution}
\item{DEup.sd }{estimated standard deviation of overexpressed distribution, set to be equal to the estimated standard deviation of the null distribution}
\item{DEdown.mean }{estimated mean of underexpressed distribution, calculated as the 5th percentile of the estimated null distribution}
\item{DEdown.sd }{estimated standard deviation of underexpressed distribution, set to be equal to the estimated standard deviation of the null distribution}
}
\author{
Alyssa Frazee
}

\seealso{
\code{\link{getParams}}
}
