\name{getAnnotation}
\alias{getAnnotation}
\title{Download exon information for a given genome}
\usage{
  getAnnotation(genome, tablename, genes = TRUE,
    verbose = TRUE)
}
\arguments{
  \item{genome}{Genome (species) for which annotation is
  desired.  A list of supported genomes can be found using
  \code{rtracklayer:::ucscGenomes()[,"db"]}; details on
  each genome can be seen using
  \code{rtracklayer:::ucscGenomes()}.}

  \item{tablename}{UCSC table from which to download exon
  information. Use \code{supportedTables(genome)} to get a
  list of supported tables for \code{genome}.}

  \item{genes}{If TRUE (as it is by default), each exon in
  the resulting data frame is labeled with the gene it
  belongs to. Gene information is not available in every
  table.  If FALSE, each exon in the resulting data frame
  is labeled with the transcript it belongs to.}

  \item{verbose}{If TRUE, updates are printed on screen as
  annotation download progresses.}
}
\value{
  A data frame (one row per exon) giving the exon's gene or
  transcript, location (start/end), and possibly other
  information.
}
\description{
  Using the \code{GenomicFeatures} package and the UCSC
  genome browser, creates a data frame of exons (one exon
  per row) for the specified genome.
}
\note{
  This function interacts with the online UCSC Genome
  Browser, so internet connection is required to use this
  function, and connection speed affects function speed.
}
\examples{
\dontrun{
mouse.exons <- getAnnotation("mm9","refGene")
head(mouse.exons)
}
}
\author{
  Alyssa Frazee
}

