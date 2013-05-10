\name{getColmeds}
\alias{getColmeds}
\title{Calculate the medians of the database.}
\usage{
  getColmeds(dbfile, tablename, colsubset = c(-1),
    nonzero = FALSE)
}
\arguments{
  \item{dbfile}{Name/location (as character string) of
  database (usually ".db") file containing nucleotide by
  sample coverage.}

  \item{tablename}{Name of the table the database contains}

  \item{colsubset}{Optional vector of column indices of the
  input file that denote samples you wish to include in
  analysis. Should NEVER include 1 (genomic position).}

  \item{nonzero}{If TRUE, use the median of only the
  nonzero counts as the library size adjustment.}
}
\value{
  A vector with the column medians.
}
\description{
  \code{getColmeds} calculates the median of the samples of
  interest from the database table. It is used in
  \code{\link{getLimmaInput}} and it's output can be saved
  to avoid calculating it multiple times when running
  \code{\link{get.pvals}}.
}
\details{
  It is assumed that the first column in the database is
  called \code{pos} and contains genomic position.  Note
  that "group" must have the one fewer entries than the
  database denoted by \code{dbfile} has columns (or, if
  colsubset is used, one fewer entries than the length of
  colsubset).
}
\examples{
## add example here when we have a vignette
}
\author{
  Alyssa Frazee. Code split from getLimmaInput by Leonardo
  Collado-Torres.
}
\seealso{
  \code{\link{getLimmaInput}}
}

