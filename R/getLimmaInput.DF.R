#'Fit a linear model to each nucleotide
#'
#'Fits the linear model log2(count+scalefac) = beta0 + beta1*group +
#'beta2*library.size + [optional confounders] to each nucleotide.  From these
#'models, this function constructs an object which can be directly passed to
#'\code{getTstats} to obtain \code{limma}'s moderated t statistics for each
#'nucleotide, which we use as a measure of strength of association between
#'group and count (expression). Reads coverage from a DataFrame which can be generated using \code{makeDF}
#' and relies heavily on the \code{limma} package, using
#'\code{lmFit} as the main workhorse.
#'
#' Note
#'\code{adjustvars} must have the same number of rows as \code{group} has
#'entries. Larger values of \code{chunksize} require more memory; smaller
#'values of \code{chunksize} require more computation time.
#'
#'@param DF A list containing a \code{DF} DataFrame with the coverage data and a logical Rle with the positions that passed the cutoff. This object is generated using \code{\link{makeDF}}.
#'@param comparison Either \code{twogroup}, \code{multigroup} or \code{expression}. \code{multigroup} will use the F-statistic and \code{expression} tests the intercept-only model.
#'@param group a 0/1 vector grouping the samples (columns) in the database.
#'@param chunksize How many rows of the merged table should be processed at a time?
#'@param adjustvars Optional matrix of adjustment variables (e.g. measured
#'confounders, output from SVA, etc.) to use in fitting linear models to each
#'nucleotide.
#'@param colsubset Optional vector of column indices of the input file that
#'denote samples you wish to include in analysis. 
#'@param scalefac A log transformation is used on the count tables, so zero
#'counts present a problem.  What number should we add to the entire matrix
#'before running the models?  Defaults to 32.
#'@param nonzero If TRUE, use the median of only the nonzero counts as the
#'library size adjustment.
#'
#'@return A list with elements
#'\item{ebobject }{A list of five vectors (\code{coefficients}, \code{stdev.unscaled}, \code{sigma}, \code{df.residual}, and \code{Amean}), mimicking the \code{MArrayLM} class in \code{limma}.  Here, \code{coefficients} and \code{stdev.unscaled} are only returned for beta1, the coefficient for \code{group}, as it is assumed this is the only covariate of interest.}
#'\item{pos }{A vector of the same length as those contained in \code{ebobject}, giving the genomic positions of each linear model.}
#'@author Alyssa Frazee, Leonardo Collado-Torres
#'@export
#'@import IRanges
#'@seealso \code{\link{getTstats}}, \code{\link{makeDF}}
#'@references Smyth G (2004).  "Linear models and empirical Bayes methods for
#'assessing differential expression in microarray experiments." Statistical
#'Applications in Genetics and Molecular Biology 3(1): Article 3.
#'@examples
#'
#'## add example here when we have a vignette
#'
getLimmaInput.DF <- function(DF, comparison = c("twogroup", "multigroup", "expression"), group, chunksize = 1e+05, adjustvars = NULL, colsubset = NULL, scalefac = 32, nonzero = FALSE){
	
	## Load required pkgs
	require("limma")
	require("IRanges")
	#require("multicore") # To be implemented later
	
	## Subset the DataFrame to use only the columns of interest
	if(!is.null(colsubset)) {
		data <- DF$DF[, colsubset]
	} else {
		data <- DF$DF
	}
	## Get the positions
	pos <- which(DF$pos)
	
	## Determine total and loop sizes
	N <- nrow(data)
	lastloop <- trunc(N/chunksize)
	
	## Fix the lastloop in case that the N is a factor of chunksize
	if(N %% chunksize == 0 & lastloop > 0)  lastloop <- lastloop - 1

	## Create model matrix
	numcol <- ncol(data)
		
	## Get the medians of the columns if you are doing differential expression:
	if(comparison!="expression"){
		if(nonzero) {
			colmeds <- sapply(data, function(y) { median(y[y > 0]) })
		} else {
			colmeds <- sapply(data, median)
		}
		
	}
	
	if(!is.null(adjustvars) & comparison=="expression") stop("expression-only test does not require adjustment variables - please set adjustvars = NULL when using comparison = \"expression\"")
  
	if(!is.null(adjustvars)){
		## Construct the model matrix 'x' when adjusted variables are given
		string1 <- ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste0("av", i, " <- adjustvars[,", i, "]")))
			string1 <- paste(string1, paste0("av", i), sep="+")
		}
		eval(parse(text=paste0("x = model.matrix(~ as.factor(group) + colmeds", string1, ")")))
	} else {
		if(comparison == "expression"){ 
			## For expression-test, add a column of 1s to the model matrix
			int <- rep(1, numcol)
			x <- model.matrix(~ 0 + int) #intercept-only model, in model.matrix form.
		}
		if(comparison == "twogroup" | comparison == "multigroup"){
			## For differential expresion
			x <- model.matrix(~ as.factor(group) + colmeds)
		}
	}	
	
	## Define modeling function to apply to each chunk:
	lmFit.apply <- function(i) {
		## Subset the DataFrame to the current chunk and transform to a regular matrix
  		if(i!=lastloop) {
			mymat.DF <- data[(chunksize * i + 1):(chunksize * (i + 1)),] 
		} else  {
			## Lastloop
			mymat.DF <- data[(chunksize * i + 1):N,] 
		}
		## Transform out of the Rle world
		mymat <- as.matrix(as.data.frame(mymat.DF))
		
		## Proceed with regular operations
  		mymat <- log2(mymat + scalefac)
  		Amean <- rowMeans(mymat) 
		
		## Finally use lmFIt
  		fit <- lmFit(mymat, x)
  		return(list(fit=fit, Amean=Amean))
  	}
	
	## Fit a model to each row (chunk) of database:
	lmFit.output <- lapply(0:lastloop, lmFit.apply)
	
	## Gather results from all the models and return them:
	if(comparison=="multigroup"){
		stop("multigroup functionality not implemented yet - check back soon!")
	} else if (comparison=="twogroup"){
		coef <- unlist(lapply(lmFit.output, function(x) { x$fit$coefficients[, 2] }))
		stdev <- unlist(lapply(lmFit.output, function(x) { x$fit$stdev.unscaled[, 2] }))    
		am <- unlist(lapply(lmFit.output, function(x) { x$fit$Amean }))   
		sma <- unlist(lapply(lmFit.output, function(x) { x$fit$sigma }))
		dfres <- unlist(lapply(lmFit.output, function(x) { x$fit$df.residual }))    
	} else if(comparison=="expression"){
		## Subtract log2(scalefac) since we want to test beta_0 = 0, not beta_0 = log2(scalefac + 0), which is what we see with 0 expression under the transformation.
		coef <- unlist(lapply(lmFit.output, function(x) { x$fit$coefficients })) - log2(scalefac)
		stdev <- unlist(lapply(lmFit.output, function(x) { x$fit$stdev.unscaled }))
		am <- unlist(lapply(lmFit.output, function(x) { x$fit$Amean })) - log2(scalefac)
		sma  <- unlist(lapply(lmFit.output, function(x) { x$fit$sigma }))
		dfres <- unlist(lapply(lmFit.output, function(x) { x$fit$df.residual }))    
	}
	
    ## Done!
	return(list(ebobject = list(coefficients = as.numeric(coef), stdev.unscaled = as.numeric(stdev), sigma = sma, df.residual = dfres, Amean = am), pos = pos, comparison = comparison))
	
}
