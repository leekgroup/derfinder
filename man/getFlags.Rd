\name{getFlags}
\alias{getFlags}
\title{Flag regions with genomic events of interest}
\usage{
  getFlags(regions, exons, chromosome, nonDE = FALSE,
    pctcut = 0.8)
}
\arguments{
  \item{regions}{data frame of regions to consider, usually
  the \code{$states} output of \code{getRegions}}

  \item{exons}{data frame containing annotated exons you
  would like to consider.  Should have columns
  \code{chr},\code{start},\code{end},\code{exon_id}, and
  \code{gene}.  Can be created with
  \code{\link{getAnnotation}}.}

  \item{chromosome}{Chromosome you're considering.
  Currently you can only do one chromosome at a time.}

  \item{nonDE}{Should we give you flags for only the
  equally expressed regions? (Usually we just want flags
  for differentially expressed regions, so in general this
  is FALSE).}

  \item{pctcut}{What percentage of an exon should a region
  overlap in order to call that exon differentially
  expressed?  Default 0.8 (meaning 80\%).}
}
\value{
  List with elements having length equal to the number of
  rows in \code{regions}: \item{flags}{Event indicated by
  corresponding region} \item{flag.info}{Genomic features
  (exons) associated with genomic events}
  \item{percent.exon}{The percent of the exon (in
  flag.info) overlapped by this region}
}
\description{
  Connects results of tornado pipeline to existing
  annotation as a way to determine which regions may point
  to interesting biological events.
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getAnnotation}}
}

