\name{get.pvals}
\alias{get.pvals}
\title{
Permutation-based significance tests for differentially expressed regions
}
\description{
Creates simulated null distribution of average t statistics for regions classified as over- or underexpressed, and obtains p-values for observed over- and underexpressed regions based on this simulated null.
}
\usage{
get.pvals(regions, dbfile, tablename, num.perms = 1, group, est.params, chromosome, colsubset = c(-1), adjustvars=NULL, nonzero=FALSE, scalefac=32, chunksize=1e+05, colmeds=NULL)
}
\arguments{
  \item{regions}{
data frame of regions to obtain p-values for: specifically the \code{$states} element of the return of \code{getRegions}.  Must have \code{mean.t} and \code{state} columns. 
}
  \item{dbfile}{
Database file used to obtain the \code{regions} data frame.
}
  \item{tablename}{
Name of the table in \code{dbfile}.
}
  \item{num.perms}{
Number of permutations to use to create the null distribution.  
}
  \item{group}{
0/1 vector denoting group labels for the samples used in analysis
}
  \item{est.params}{
return from \code{getParams} - HMM parameters used when creating \code{regions}
}
  \item{chromosome}{
Chromosome you are analyzing.  Currently only runs one chromosome at a time.
}
\item{colsubset}{
Colsubset argument used to create \code{regions}
}
  \item{adjustvars}{
Optional matrix of adjustment variables (e.g. measured confounders, output from SVA, etc.) to use in fitting linear models to each nucleotide.
}
\item{nonzero}{
If TRUE, use the median of only the nonzero counts as the library size adjustment.
}
\item{scalefac}{
A log transformation is used on the count tables, so zero counts present a problem.  What number should we add to the entire matrix before running the models?  Defaults to 32.
}
  \item{chunksize}{
How many rows of the database should be processed at a time? 
}
\item{colmeds}{
If NULL, the column medians are calculated using \code{\link{getColmeds}}. Otherwise, the output of \code{\link{getColmeds}} is expected.
}
}
\value{
vector having length equal to the number of rows in \code{regions}, giving a p-value for each region of state 3 or 4 in \code{regions}.
}
\author{
Alyssa Frazee
}

\seealso{
\code{\link{getLimmaInput}},\code{\link{getTstats}},\code{\link{getRegions}}
}
