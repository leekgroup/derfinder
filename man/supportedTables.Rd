\name{supportedTables}
\alias{supportedTables}
\title{Print list of supported (downloadable) tables for a given genome}
\usage{
  supportedTables(genome)
}
\arguments{
  \item{genome}{Character string giving the genome of
  interest.  Available genomes can be seen with
  \code{rtracklayer:::ucscGenomes()[,"db"]}.}
}
\value{
  prints a list of available tables for \code{genome}.
}
\description{
  To be used with \code{getAnnotation}: provides a list of
  available tables from UCSC for any given genome.
}
\examples{
supportedTables("mm10")
mouse.exons = getAnnotation("mm10","refGene") #refGene appears in printed output of supportedTables("mm10").
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getAnnotation}}
}

