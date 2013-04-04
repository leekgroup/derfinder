\name{getTstats}
\alias{getTstats}
\title{ calculate moderated t statistic for each nucleotide
}
\description{
Modified (slimmed down) version of limma's eBayes, which uses a shrinkage approach to moderate t statistics from tests for differential expression across the genome.
}
\usage{
getTstats(fit, trend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
list with elements $sigma, $df.residual, $Amean (if trend=TRUE), $coefficients, and $stdev.unscaled.  Usually the output of getLimmaInput (the $ebobject item), a more memory-efficient version of lmFit.
}
  \item{trend}{
as in eBayes, if TRUE, allow an intensity trend in the priors of variances (across genes/bps)
}
}
\value{
list with elements 
\item{tt}{moderated t statistics for each gene/bp}
\item{logfchange}{estimated log2 fold change on coverage between the groups, as estimated by a linear model}
}
\references{
Smyth G (2004).  “Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.”  Statistical Applications in Genetics and Molecular Biology 3(1): Article 3.
}
\author{
Alyssa Frazee
}
\seealso{
\code{\link{getLimmaInput}}
}
