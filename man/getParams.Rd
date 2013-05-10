\name{getParams}
\alias{getParams}
\title{Calculate parameters to use as input for HMM}
\usage{
  getParams(tstats, plots = FALSE, plotfile = NULL,
    verbose = F)
}
\arguments{
  \item{tstats}{Vector containing all moderated t
  statistics obtained using \code{getTstats}.}

  \item{plots}{if TRUE, create diagnostic plots as
  parameters are estimated}

  \item{plotfile}{Optional string giving a location and PDF
  file name to which plots should be written, if
  \code{plots = TRUE}.  If NULL, plots are created in the
  available graphics device.}

  \item{verbose}{If TRUE, periodic messages are printed
  onscreen during estimation.}
}
\value{
  A list with elements \item{params }{list with elements
  \code{mean} and \code{sd}, both 4-item vectors.
  \code{mean} gives the respective means of the "not
  expressed," "equally expressed," "overexpressed," and
  "underexpressed" distributions; \code{sd} gives their
  respective standard deviations.} \item{stateprobs
  }{vector of percentages of the mixture distribution that
  come from the not expressed," "equally expressed,"
  "overexpressed," and "underexpressed" distributions,
  respectively. It is assumed that "overexpressed" and
  "underexpressed" t statistics comprise equal percentages
  of the mixture.}
}
\description{
  Assumes that the moderated t statistics obtained by
  fitting a linear model to each nucleotide come from a
  Gaussian mixture distribution, where the four
  distributions in the mixture represent distributions of t
  statistics from "underexpressed," "overexpressed,"
  "equally expressed," and "not expressed" nucleotides.
  \code{getParams} estimates the parameters of each of the
  sub-distributions, as well as the percentage of the
  mixture distribution each contributes, in order to use
  these parameters to fit a Hidden Markov Model that
  classifies the nucleotides.
}
\details{
  The standard pipeline here is to feed the output from
  \code{getParams} directly into \code{getRegions} using
  the "HMM" option.
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getRegions}}
}

