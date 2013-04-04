\name{supportedTables}
\alias{supportedTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ print list of supported (downloadable) tables for a given genome
}
\description{
To be used with \code{getAnnotation}: provides a list of available tables from UCSC for any given genome. 
}
\usage{
supportedTables(genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
Character string giving the genome of interest.  Available genomes can be seen with \code{rtracklayer:::ucscGenomes()[,"db"]}. 
}
}
\value{
prints a list of available tables for \code{genome}.
}

\author{
Alyssa Frazee
}

\seealso{
\code{\link{getAnnotation}}
}
\examples{
supportedTables("mm10")
mouse.exons = getAnnotation("mm10","refGene") #refGene appears in printed output of supportedTables("mm10").
}
