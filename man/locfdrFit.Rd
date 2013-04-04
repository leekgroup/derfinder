\name{locfdrFit}
\alias{locfdrFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ modified version of Efron's locfdr
}
\description{
Compute local false discovery rates, following the definitions and description in references listed below. Exactly the same as \code{locfdr}, but returns extra information.}
\usage{
locfdrFit(zz, bre = 120, df = 7, pct = 0, pct0 = 1/4, nulltype = 1, type = 0, plot = 1, mult, mlests, main = " ", sw = 0, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{
A vector of summary statistics, one for each case under simultaneous consideration. The calculations assume a large number of cases, say length(zz) exceeding 200. Results may be improved by transforming zz so that its elements are theoretically distributed as N(0, 1) under the null hypothesis. See the \code{locfdr} vignette for tips on creating \code{zz}.}
  \item{bre}{
Number of breaks in the discretization of the z-score axis, or a vector of breakpoints fully describing the discretization. If length(zz) is small, such as when the number of cases is less than about 1000, set \code{bre} to a number lower than the default of 120. The tornado package keeps this at its default.}
  \item{df}{
Degrees of freedom for ﬁtting the estimated density f(z).  The tornado package keeps this at its default.
}
  \item{pct}{
Excluded tail proportions of zz’s when ﬁtting f(z). \code{pct=0} includes full range
of zz’s. \code{pct} can also be a 2-vector, describing the ﬁtting range. The tornado package keeps this at its default.
}
  \item{pct0}{
Proportion of the zz distribution used in ﬁtting the null density f0(z) by central
matching. If a 2-vector, e.g. \code{pct0=c(0.25,0.60)}, the range \code{[pct0[1], pct0[2]]} is used. If a scalar, \code{[pct0, 1-pct0]} is used. The tornado package keeps this at its default.
}
  \item{nulltype}{
Type of null hypothesis assumed in estimating f0(z), for use in the fdr calculations. 0 is the theoretical null N(0; 1), 1 is maximum likelihood estimation, 2 is central matching estimation, 3 is a split normal version of 2. The tornado package fixes this at 1.
}
  \item{type}{
Type of ﬁtting used for f; 0 is a natural spline, 1 is a polynomial, in either case with degrees of freedom \code{df} [so total degrees of freedom including the intercept is \code{df+1}.] The tornado package fixes this at 0.
}
  \item{plot}{
Plots desired. 0 gives no plots. 1 gives single plot showing the histogram of zz and ﬁtted densities f and p0  f0. 2 also gives plot of fdr, and the right and left tail area Fdr curves. 3 gives instead the f1 cdf of the estimated fdr curve; plot=4 gives all three plots. The tornado package allows choices 0 and 1; equivalent to \code{plots = F} and \code{plots = T} in \code{getParams}.
}
  \item{mult}{
Optional scalar multiple (or vector of multiples) of the sample size for calculation of the corresponding hypothetical Efdr value(s). This is not used in the tornado package.}
  \item{mlests}{
Optional vector of initial values for (delta0, sigma0) in the maximum likelihood iteration. The tornado package includes these values in an updated run of \code{locfdrFit} if they are suggested via warning in the first run.
}
  \item{main}{
Main heading for the histogram plot when \code{plot}>0.
}
  \item{sw}{
Determines the type of output desired. 2 gives a list consisting of the last 5 values listed under Value below. 3 gives the square matrix of dimension \code{bre}-1 representing the inﬂuence function of log(fdr). Any other value of sw returns a list consisting of the ﬁrst 5 (6 if \code{mult} is supplied) values listed below. The tornado package fixes this at 0.}
  \item{verbose}{
if TRUE, various messages are printed onscreen during \code{getParams}.
}
}
\details{
Generally not used directly in the tornado package, but is a workhorse for \code{getParams}.
}
\value{
See the \code{locfdr} manual for returned values. \code{locfdrFit} returns the following additional elements:
\item{yt }{bin heights}
\item{mlest.lo }{if a warning that \code{mlest} should be used in a re-run of \code{locfdrFit}, the suggested first element of \code{mlest}.}
\item{mlest.hi }{if a warning that \code{mlest} should be used in a re-run of \code{locfdrFit}, the suggested second element of \code{mlest}.}
\item{needsfix }{0 if no warning to fix the run with \code{mlest} parameters, 1 otherwise}
\item{nulldens }{a rough estimate of y-axis values for f0(x)}
\item{fulldens }{a rough estimate of y-axis values for f(x)}
}
\references{
http://cran.r-project.org/web/packages/locfdr/locfdr.pdf
}
\author{
Bradley Efron, slight modifications by Alyssa Frazee
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{getParams}
}
