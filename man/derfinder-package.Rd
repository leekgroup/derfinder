\name{derfinder-package}
\alias{derfinder-package}
\alias{derfinder}
\docType{derfinder}
\title{
Tornado
}
\description{
Differential expression analysis of per-nucleotide coverage tables.
}
\details{
\tabular{ll}{
Package: \tab derfinder\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-29\cr
License: \tab What license is it under?\cr
}
Creates SQLite database of per-nucleotide coverage files (\code{makeDb}), fits linear model to each nucleotide to determine differential expression (\code{getLimmaInput},\code{getTstats}), fits Hidden Markov Model using moderated t statistics from linear models as state emissions to get a list of differentially expressed regions (\code{getRegions}), connects found regions with known annotation (\code{getAnnotation}, \code{plotRegion}, \code{plotExon}, \code{plotGene}).
}
\author{
Alyssa Frazee <afrazee@jhsph.edu>
}
\references{
Paper coming soon.
}
\keyword{ package }
