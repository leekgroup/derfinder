\name{mergeRegions}
\alias{mergeRegions}
\title{
Merge close regions
}
\description{
Regions of the same state that are fewer than 5bp apart from each other are merged.
}
\usage{
mergeRegions(regions, cutoff = 6)
}
\arguments{
  \item{regions}{
data frame of regions, usually the \code{$state} return of \code{getRegions}
}
  \item{cutoff}{
How far apart do regions need to be in order to *not* be merged?  (Defaults to 6, meaning that regions fewer than 5bp apart *are* merged).
}
}
\value{
data frame (comparable to \code{regions} input) with close regions merged.
}
\author{
Alyssa Frazee
}

\seealso{
\code{\link{getRegions}}
}
