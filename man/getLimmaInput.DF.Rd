\name{getLimmaInput.DF}
\alias{getLimmaInput.DF}
\title{Fit a linear model to each nucleotide}
\usage{
  getLimmaInput.DF(DF,
    comparison = c("twogroup", "multigroup", "expression"),
    group, chunksize = 1e+05, adjustvars = NULL,
    colsubset = NULL, scalefac = 32, nonzero = FALSE)
}
\arguments{
  \item{DF}{A list containing a \code{DF} DataFrame with
  the coverage data and a logical Rle with the positions
  that passed the cutoff. This object is generated using
  \code{\link{makeDF}}.}

  \item{comparison}{Either \code{twogroup},
  \code{multigroup} or \code{expression}. \code{multigroup}
  will use the F-statistic and \code{expression} tests the
  intercept-only model.}

  \item{group}{a 0/1 vector grouping the samples (columns)
  in the database.}

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
  A list with elements \item{ebobject }{A list of five
  vectors (\code{coefficients}, \code{stdev.unscaled},
  \code{sigma}, \code{df.residual}, and \code{Amean}),
  mimicking the \code{MArrayLM} class in \code{limma}.
  Here, \code{coefficients} and \code{stdev.unscaled} are
  only returned for beta1, the coefficient for
  \code{group}, as it is assumed this is the only covariate
  of interest.} \item{pos }{A vector of the same length as
  those contained in \code{ebobject}, giving the genomic
  positions of each linear model.}
}
\description{
  Fits the linear model log2(count+scalefac) = beta0 +
  beta1*group + beta2*library.size + [optional confounders]
  to each nucleotide.  From these models, this function
  constructs an object which can be directly passed to
  \code{getTstats} to obtain \code{limma}'s moderated t
  statistics for each nucleotide, which we use as a measure
  of strength of association between group and count
  (expression). Reads coverage from a DataFrame which can
  be generated using \code{makeDF} and relies heavily on
  the \code{limma} package, using \code{lmFit} as the main
  workhorse.
}
\details{
  Note \code{adjustvars} must have the same number of rows
  as \code{group} has entries. Larger values of
  \code{chunksize} require more memory; smaller values of
  \code{chunksize} require more computation time.
}
\examples{
## add example here when we have a vignette
}
\author{
  Alyssa Frazee, Leonardo Collado-Torres
}
\references{
  Smyth G (2004).  "Linear models and empirical Bayes
  methods for assessing differential expression in
  microarray experiments." Statistical Applications in
  Genetics and Molecular Biology 3(1): Article 3.
}
\seealso{
  \code{\link{getTstats}}, \code{\link{makeDF}}
}

