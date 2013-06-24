\name{strapply}
\alias{strapply}
\title{Helper function for read.csv.sql}
\usage{
  strapply(X, pattern, FUN = function(x, ...) x,
    backref = NULL, ..., empty = NULL, ignore.case = FALSE,
    perl = TRUE, engine = "R", simplify = FALSE,
    USE.NAMES = FALSE, combine = c)
}
\arguments{
  \item{X}{list or (atomic) vector of character strings to
  be used}

  \item{pattern}{character string containing a regular (or
  character string for \code{fixed=TRUE} to be matched in
  the given character vector.}

  \item{FUN}{a function, formula, character string, list or
  proto object to be applied to each element of \code{X}.
  See discussion in \code{\link{gsubfn}.}}

  \item{backref}{see \code{gsubfn}}

  \item{\dots}{optional arguments to \code{gsubfn}}

  \item{empty}{If there is no match to a string return this
  value.}

  \item{ignore.case}{If TRUE then case is ignored in the
  \code{pattern} argument.}

  \item{perl}{If TRUE then \code{engine="R"} is used with
  perl regular expressions.  It is required to keep this
  argument at TRUE, since \code{tcl} engine capabilities
  have been removed from this function.}

  \item{engine}{Should always be set to \code{"R"}, since
  the \code{tcl} engine is not available in the tornado
  package.}

  \item{simplify}{logical or function. If logical, should
  the result be simplified to a vector or matrix, as in
  \code{sapply} if possible? If function, that function is
  applied to the result with each component of the result
  passed as a separate argument. Typically if the form is
  used it will typically be specified as rbind.}

  \item{USE.NAMES}{logical; if \code{TRUE} and if \code{X}
  is character, use \code{X} as 'names' for the result
  unless it had names already.}

  \item{combine}{combine is a function applied to the
  components of the result of \code{FUN}. The default is
  \code{"c"}. \code{"list"} is another common choice. The
  default may change to be \code{"list"} in the future.}
}
\value{
  A list of character strings
}
\description{
  The same as \code{strapply} in the \code{gsubfn} package,
  but with \code{tcltk} capabilities removed.
}
\details{
  See details in \code{gsubfn} package.
}
\note{
  Does not need to be used directly in tornado;
  \code{makeDb} wraps this entirely.
}
\author{
  G. Grothendieck
}
\references{
  http://cran.r-project.org/web/packages/gsubfn/gsubfn.pdf
}
\seealso{
  \code{\link{makeDb}}
}

