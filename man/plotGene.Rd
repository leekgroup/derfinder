\name{plotGene}
\alias{plotGene}
\title{Plot pipeline data/results for a given gene}
\usage{
  plotGene(getRegionObject, ind = NULL, genename = NULL,
    tstats, pos, annotation, counts, group, bppad = 50,
    axpad = 50, prettyskips = T, skiplines = T,
    countsheader = F, countssep = "\t", tabname = NULL,
    plotfile = NULL, width = 900, height = 750,
    plottitle = NULL, chromosome, legendloc = "bottomleft",
    scalefac = 0.5, ylim = c(0, 9))
}
\arguments{
  \item{getRegionObject}{The name of an object created with
  \code{\link{getRegions}}.}

  \item{ind}{index in the provided \code{annotation} of the
  exon you wish to plot.}

  \item{genename}{name of the gene (as listed in the
  provided \code{annotation}) you wish to plot}

  \item{tstats}{Vector of t-statistics that was used to
  create \code{getRegionObject}.}

  \item{pos}{vector of genomic positions corresponding to
  \code{tstats}}

  \item{annotation}{data frame containing exon annotation
  to use (see \code{\link{getAnnotation}}). Must contain a
  "geneName" column listing the exon names. (The column
  \code{gene} in a \code{getAnnotation} object can be
  re-named to \code{geneName}).}

  \item{counts}{Raw coverage data used to obtain
  \code{tstats}. This can be provided in one of three
  forms: (1) a string indicating the location/file name of
  a SQLite database containing the counts, usually created
  with \code{makeDb}; (2) a string indicating the
  location/file name of a text file containing coverage --
  this will get loaded into memory!! -- or (3) an
  already-loaded matrix containing the raw data. Note that
  \code{counts} must have the same number of rows as
  \code{pos} and \code{tstats} have elements, and the rows
  must correspond to genomic position \code{pos}.}

  \item{group}{a vector containing the group labels for the
  columns of counts. Only 2 groups are permitted at this
  time. These labels are used in the plot's legend, so
  generally they are character strings (rather than, say,
  0/1).}

  \item{bppad}{the number of bases to plot outside of the
  designated region (default 50).  Essentially, use this to
  "zoom" in (decrease bppad) or out (increase bppad) on the
  plotted region.}

  \item{axpad}{how much wider (in bases) you'd like the
  x-axis to be, compared to the plotted area}

  \item{prettyskips}{If TRUE, plot
  counts/states/t-statistics contiguously, even if there
  are zero entries between them (i.e., even though pos may
  not indicate that contiguous postions are being plotted).
  Note that in general, when plotting just one exon, this
  is not an issue as exons tend to contain contiguous data.
  So, \code{prettyskips} will only affect areas outside the
  exon, i.e., \code{prettyskips} has larger impact if bppad
  is large. Also, note that \code{prettyskips = FALSE} is
  not allowed at this time.}

  \item{skiplines}{if TRUE, add a light vertical line to
  the plot indicating an eliminated "low-coverage"
  nucleotide}

  \item{countsheader}{If TRUE, the \code{counts} matrix
  contains a header row. Not usually the case if counts is
  a database or already-loaded matrix.}

  \item{countssep}{If reading counts from a text file, the
  separator used in that file.}

  \item{tabname}{If counts is a database file, the name of
  the table that was dumped into that database. (See
  \code{tablename} in \code{makeDb})}

  \item{plotfile}{Optional string containing a file you'd
  like to put the plot into (if NULL, plot appears
  interactively). Should have a .jpg extension.}

  \item{width}{Only used when \code{plotfile} is non-null:
  width (in pixels) of resulting jpg.  Defaults to 900.}

  \item{height}{Only used when \code{plotfile} is non-null:
  width (in pixels) of resulting jpg.  Defaults to 750.}

  \item{plottitle}{Optional main title to use on your plot.
  Defaults to chromosome: start-end (referring to plotted
  REGION)}

  \item{chromosome}{The chromosome corresponding to the
  exon being plotted, in the same format as chromosomes are
  listed in the supplied \code{annotation}.}

  \item{legendloc}{Can be one of
  "topright","bottomright","topleft",or "bottomleft"
  indicating where the legend (indicating group label on
  raw count plot) should be located. Defaults to
  "bottomleft."}

  \item{scalefac}{How should we offset counts (so that we
  can log-transform everything, even zero counts)?
  Defaults to 0.5, derfinder uses an offset of 32.}

  \item{ylim}{Equivalent to ylim argument in plot: what
  should the lower and upper limits be on the y-axis of the
  top panel of plot be?  Defaults to c(0,9).}
}
\value{
  Plots (either on screen or in the supplied .jpg file) the
  following: top panel = genomic position vs. raw coverage
  data, middle panel = genomic position vs. moderated t
  statistic from linear model at that position, bottom
  panel = genomic position vs. predicted state for that
  position, with annotated exons overlaid. States are as
  follows: gray = not expressed, black = equally expressed,
  red = overexpressed (in whatever group = 1 represented in
  \code{getRegions}), green = underexpressed.
}
\description{
  Creates a 3-paneled plot of a selected gene: panel 1 =
  genomic position vs. raw coverage data, panel 2 = genomic
  position vs. moderated t statistic from linear model at
  that position, panel 3 = genomic position vs. predicted
  state for that position, with annotated exons overlaid.
}
\note{
  Provide exactly one of \code{ind} and \code{genename}.
  Recommendation is to provide \code{geneName}.
}
\author{
  Alyssa Frazee
}
\seealso{
  \code{\link{getRegions}}, \code{\link{makeDb}}
}

