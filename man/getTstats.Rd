\name{getTstats}
\alias{getTstats}
\title{Calculate moderated t statistic for each nucleotide}
\usage{
  getTstats(fit, trend = FALSE)
}
\arguments{
  \item{fit}{list with elements $sigma, $df.residual,
  $Amean (if trend=TRUE), $coefficients, and
  $stdev.unscaled.  Usually the output of getLimmaInput
  (the $ebobject item), a more memory-efficient version of
  lmFit.}

  \item{trend}{as in eBayes, if TRUE, allow an intensity
  trend in the priors of variances (across genes/bps)}
}
\value{
  A list with elements \item{tt}{moderated t statistics for
  each gene/bp} \item{logfchange}{estimated log2 fold
  change on coverage between the groups, as estimated by a
  linear model}
}
\description{
  Modified (slimmed down) version of limma's eBayes, which
  uses a shrinkage approach to moderate t statistics from
  tests for differential expression across the genome.
}
\author{
  Alyssa Frazee
}
\references{
  Smyth G (2004).  “Linear models and empirical Bayes
  methods for assessing differential expression in
  microarray experiments.” Statistical Applications in
  Genetics and Molecular Biology 3(1): Article 3.
}
\seealso{
  \code{\link{getLimmaInput}}
}

