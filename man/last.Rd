\name{last}
\alias{last}
\title{Get the last element}
\usage{
  last(x)
}
\arguments{
  \item{x}{any object with a \code{tail} method}
}
\value{
  same as \code{x[length(x)]}
}
\description{
  The function returns the last element in a given object
  (any object with a \code{tail} method will work).
}
\note{
  Helper function for \code{getParams}
}
\examples{
x = c(1:20)
last(x) == 20
}
\author{
  Alyssa Frazee
}

