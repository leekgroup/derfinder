\name{makeDF}
\alias{makeDF}
\title{Merge the coverage information for a group of samples}
\usage{
  makeDF(chr, datadir = NULL, sampledirs = NULL,
    samplepatt = NULL, cutoff = 5, chrlen = NULL,
    org = "BSgenome.Hsapiens.UCSC.hg19",
    bamterm = "accepted_hits.bam", output = NULL,
    verbose = TRUE)
}
\arguments{
  \item{chr}{Chromosome to read. Should be in simple
  format: X and not chrX.}

  \item{datadir}{The main directory where each of the
  \code{sampledirs} is a sub-directory of \code{datadir}.}

  \item{sampledirs}{A character vector with the names of
  the sample directories. If \code{datadir} is \code{NULL}
  it is then assumed that \code{sampledirs} specifies the
  full path to each sample.}

  \item{samplepatt}{If specified and \code{sampledirs} is
  set to \code{NULL}, then the directories matching this
  pattern in \code{datadir} (set to \code{.} if it's set to
  \code{NULL}) are used as the sample directories.}

  \item{cutoff}{Per base pair, at least one sample has to
  have coverage greater than \code{cutoff} to be included
  in the result.}

  \item{chrlen}{The chromosome length in base pairs.}

  \item{org}{If \code{chrlen} is set to \code{NULL}, then
  the chromosome length is deduced using the BSgenome
  annotation packages for either Human or Mouse. This is
  just meant as an easy to use option for a handful of
  genomes.}

  \item{bamterm}{Name of the BAM file used in each sample.
  By default it is set to \code{accepted_hits.bam} since
  that is the automatic name generated when aligning with
  TopHat.}

  \item{output}{If \code{NULL} then no output is saved in
  disk. If \code{auto} then an automatic name is
  constructed (chrXDF.Rdata for example). If another
  character is specified, then that name is used for the
  output file.}

  \item{verbose}{If \code{TRUE} basic status updates will
  be printed along the way.}
}
\value{
  A list with two objects. The first one, \code{DF}, is a
  DataFrame object where each column represents a sample.
  The number of rows depends on the number of base pairs
  that passed the cutoff. The second one, \code{pos}, is a
  logical Rle with the positions of the chromosome that
  passed the cutoff.
}
\description{
  For a group of samples this function reads the coverage
  information for a specific chromosome directly from the
  BAM files. It then merges them into a DataFrame and
  removes the bases that do not pass the cutoff.
}
\author{
  Leonardo Collado-Torres
}

