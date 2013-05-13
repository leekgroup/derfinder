\name{findclosest.behind}
\alias{findclosest.behind}
\title{Find closest long region upstream}
\usage{
  findclosest.behind(ind, regions)
}
\arguments{
  \item{ind}{index of regions data frame}

  \item{regions}{the $states data frame, as returned by
  \code{getRegions}}
}
\value{
  index of closest large region upstream of a region that
  needs to be merged
}
\description{
  Helper function for \code{mergeRegions}
}
\note{
  Not generally used alone - internal function for
  \code{mergeRegions}
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{mergeRegions}}
}

