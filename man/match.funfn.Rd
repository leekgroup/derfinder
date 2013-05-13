\name{match.funfn}
\alias{match.funfn}
\title{Generic extended version of R match.fun}
\usage{
  match.funfn(FUN, descend = TRUE)
}
\arguments{
  \item{FUN}{Function, character name of function or
  formula describing function.}

  \item{descend}{logical; control whether to search past
  non-function objects.}
}
\value{
  Returns a function.
}
\description{
  A generic \code{match.fun}.
}
\details{
  The default method is the same as \code{match.fun} and
  the \code{formula} method is the same as
  \code{as.function.formula}.  This function can be used
  within the body of a function to convert a function
  specification whether its a function, character string or
  formula into an actual function.
}
\note{
  This is exactly the same as \code{match.funfn} in the
  \code{gsubfn} package.  Do not load the \code{gsubfn}
  package to use this function, as \code{gsubfn} loads
  \code{tcltk}, which is not advised.
}
\author{
  G. Grothendieck
}

