\name{getF}
\alias{getF}
\title{F statistic from comparing nested limma models}
\usage{
  getF(fit, fit0, theData)
}
\arguments{
  \item{fit}{MArrayLM object from fitting model B}

  \item{fit0}{MArrayLM object from fitting model A}

  \item{theData}{data matrix used to fit both models.
  (matrix of outcome values)}
}
\value{
  data frame with elements \item{fstat }{the F statistic
  from the test} \item{df1 }{first df for F distribution
  used in testing} \item{df2 }{second df for F distribution
  used in testing} \item{f_pval }{p-value resulting from
  said test.}
}
\description{
  Given two model fits (created with limma's lmFit), with
  model A nested in model B, calculate an F statistic for
  the test for whether model B is significantly better than
  model A.
}
\author{
  Andrew Jaffe
}
\seealso{
  \code{\link{getTstats}}, \code{\link{makeDb}}
}

