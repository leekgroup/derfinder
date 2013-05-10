\name{getParams.failsafe}
\alias{getParams.failsafe}
\title{Helper function for getParams}
\usage{
  getParams.failsafe(null.mean, null.sd)
}
\arguments{
  \item{null.mean}{Estimated mean of null distribution
  (usually from \code{locfdrFit})}

  \item{null.sd}{Estimated standard deviation of null
  distribution (usually from \code{locfdrFit})}
}
\value{
  A list with elements \item{DEup.mean }{estimated mean of
  overexpressed distribution, calculated as the 95th
  percentile of the estimated null distribution}
  \item{DEup.sd }{estimated standard deviation of
  overexpressed distribution, set to be equal to the
  estimated standard deviation of the null distribution}
  \item{DEdown.mean }{estimated mean of underexpressed
  distribution, calculated as the 5th percentile of the
  estimated null distribution} \item{DEdown.sd }{estimated
  standard deviation of underexpressed distribution, set to
  be equal to the estimated standard deviation of the null
  distribution}
}
\description{
  When numerical methods \code{find.mean} and
  \code{find.sd} fail, \code{getParams.failsafe} is used to
  calculate parameters of the distributions of t statistics
  originating from over- or underexpressed nucleotides.
}
\details{
  For experienced users/debugging only. Most users should
  use \code{getParams} directly.
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getParams}}
}

