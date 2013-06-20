\name{ostrapply}
\alias{ostrapply}
\title{Helper function for read.csv.sql}
\usage{
  ostrapply(X, pattern, FUN = function(x, ...) x,
    ignore.case = FALSE, ..., empty = NULL,
    simplify = FALSE, USE.NAMES = FALSE, combine = c)
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

  \item{ignore.case}{If \code{TRUE} then case is ignored in
  the \code{pattern} argument.}

  \item{\dots}{optional arguments to \code{gsubfn}.}

  \item{empty}{If there is no match to a string return this
  value.}

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
  A list of character strings.
}
\description{
  Same as the \code{strapply} or \code{ostrapply} functions
  in \code{gsubfn} package, but with \code{tcltk}
  capabilities removed.
}
\details{
  See \code{strapply} in \code{gsubfn}.  In the tornado
  package this function should rarely be called on its own.
  It is an internal process of \code{read.csv.sql}, which
  is an internal process of \code{makeDb}.
}
\author{
  G. Grothendieck
}
\references{
  http://cran.r-project.org/web/packages/gsubfn/gsubfn.pdf
}
\seealso{
  \code{\link{makeDb}},\code{link{read.csv.sql}}
}

