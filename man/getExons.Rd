\name{getExons}
\alias{getExons}
\title{Find closest exon(s) to a genomic region}
\usage{
  getExons(region, annotation, verbose = TRUE)
}
\arguments{
  \item{region}{length-3 vector \code{(chromosome, start
  position, end position)} of the genomic region of
  interest.  Note that \code{chromosome} needs to be in the
  same format as the \code{chr} column of
  \code{annotation}.}

  \item{annotation}{Data frame containing exon information
  (one row per exon) for the appropriate genome.  Must
  contain columns \code{chr}, \code{start}, and \code{end}.
  It is recommended that \code{getAnnotation} be used to
  obtain an annotation data frame.}

  \item{verbose}{If TRUE, prints output messages when
  function finishes.}
}
\value{
  A list with elements \item{region }{the region argument
  provided} \item{closestExons }{the rows of
  \code{annotation} corresponding to the closest exon to
  \code{region}} closest exon to \code{region}
}
\description{
  Given any genomic region (chromosome, start, end), return
  the closest known exon.
}
\examples{
\dontrun{
exons = getAnnotation("hg19","knownGene")
theRegion = c("chr22", 18216902, 18218350)
getExons(theRegion, exons)
foo = getExons(theRegion, exons)
foo
foo$closestExons
}
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getAnnotation}}
}

