#### NOTE that I (alyssa) made this for comparing two models - NOT just for comparing to null model - so that the library size adjustment/other confounders could be taken into account when as part of the "null" model here.

#'F statistic from comparing nested limma models
#'
#'Given two model fits (created with limma's lmFit), with model A nested in model B,
#' calculate an F statistic for the test for whether model B is significantly better 
#' than model A.
#'
#'@param fit MArrayLM object from fitting model B
#'@param fit0 MArrayLM object from fitting model A
#'@param theData data matrix used to fit both models. (matrix of outcome values)
#'@return data frame with elements
#'\item{fstat }{the F statistic from the test}
#'\item{df1 }{first df for F distribution used in testing}
#'\item{df2 }{second df for F distribution used in testing}
#'\item{f_pval }{p-value resulting from said test.}
#'@author Andrew Jaffe
#'@export
#'@seealso \code{\link{getTstats}}, \code{\link{makeDb}}
getF = function(fit, fit0, theData) {
  
  rss1 = rowSums((fitted(fit)-theData)^2)
  df1 = ncol(fit$coef)
  rss0 = rowSums((fitted(fit0)-theData)^2)
  df0 = ncol(fit0$coef)
  
  fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
  f_pval = pf(fstat, df1-df0, ncol(theData)-df1, lower.tail=FALSE)
  fout = cbind(fstat,df1-df0, ncol(theData)-df1, f_pval)
  colnames(fout)[2:3] = c("df1","df2")
  fout = data.frame(fout)
  return(fout)
}
