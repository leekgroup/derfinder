\name{last}
\alias{last}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ get the last element
%%  ~~function to do ... ~~
}
\description{
return last element in a given object (any object with a \code{tail} method will work).
}
\usage{
last(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
any object with a \code{tail} method
}
}
\value{
same as \code{x[length(x)]}
}
\author{
Alyssa Frazee
}
\note{
Helper function for \code{getParams}
}
\examples{
x = c(1:20)
last(x) #returns 20
}