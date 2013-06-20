\name{get.pvals.DF}
\alias{get.pvals.DF}
\title{Permutation-based significance tests for differentially expressed regions}
\usage{
  get.pvals.DF(regions, num.perms = 1, est.params,
    chromosome, verbose = FALSE, DF,
    comparison = c("twogroup", "multigroup", "expression"),
    group, chunksize = 1e+05, adjustvars = NULL,
    colsubset = NULL, scalefac = 32, nonzero = FALSE)
}
\arguments{
  \item{regions}{data frame of regions to obtain p-values
  for: specifically the \code{$states} element of the
  return of \code{getRegions}.  Must have \code{mean.t} and
  \code{state} columns.}

  \item{num.perms}{Number of permutations to use to create
  the null distribution.}

  \item{est.params}{return from \code{getParams} - HMM
  parameters used when creating \code{regions}}

  \item{chromosome}{Chromosome you are analyzing.
  Currently only runs one chromosome at a time.}

  \item{verbose}{If \code{TRUE} statements with status
  updates are printed at each permutation.}

  \item{DF}{A list containing a \code{DF} DataFrame with
  the coverage data and a logical Rle with the positions
  that passed the cutoff. This object is generated using
  \code{\link{makeDF}}.}

  \item{comparison}{Either \code{twogroup},
  \code{multigroup} or \code{expression}. \code{multigroup}
  will use the F-statistic and \code{expression} tests the
  intercept-only model.}

  \item{group}{0/1 vector denoting group labels for the
  samples used in analysis. Should have the same length as
  \code{colsubset} if \code{colsubset} is specified.}

  \item{chunksize}{How many rows of the merged table should
  be processed at a time?}

  \item{adjustvars}{Optional matrix of adjustment variables
  (e.g. measured confounders, output from SVA, etc.) to use
  in fitting linear models to each nucleotide.}

  \item{colsubset}{Optional vector of column indices of the
  input file that denote samples you wish to include in
  analysis.}

  \item{scalefac}{A log transformation is used on the count
  tables, so zero counts present a problem.  What number
  should we add to the entire matrix before running the
  models?  Defaults to 32.}

  \item{nonzero}{If TRUE, use the median of only the
  nonzero counts as the library size adjustment.}
}
\value{
  A vector having length equal to the number of rows in
  \code{regions}, giving a p-value for each region of state
  3 or 4 in \code{regions}.
}
\description{
  Creates simulated null distribution of average t
  statistics for regions classified as over- or
  underexpressed, and obtains p-values for observed over-
  and underexpressed regions based on this simulated null.
}
\author{
  Alyssa Frazee, Leonardo Collado-Torres
}
\seealso{
  \code{\link{getLimmaInput.DF}},\code{\link{getTstats}},\code{\link{getRegions}}
}

