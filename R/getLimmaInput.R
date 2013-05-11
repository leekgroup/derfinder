# updated to allow for not just testing one coefficient (i.e., now can test for expression [test intercept] or test multiple groups [F statistic]).

# Minor modification: do not read in all the medians from the columns that won't be used later on

## updated 7/23/12 AF
## database file (containing filtered coverage) --> precursors to eBayes/getting moderated t-statistics.
## i.e. in pipeline: this goes after makeDb.R, before getTstats.R

## getLimmaInput()
## arguments:
## --dbfile: string giving location of .db file (usually created with makeDb()) containing the data
## --tablename: name of table in the .db file (usually created with makeDb())
## (note: dbfile and tablename arguments in makeDb() and this function should match)
## --adjustvars: if desired, an nxp matrix of adjustment variables (confounders, SVA output, etc.)
## --group: a length n 0/1 vector grouping the samples (currently only handles 2 groups)
## --colsubset: if desired, column indices of the input file of the samples you wish to include in analysis. Should NOT include 1 (pos).
## --scalefac: how much should counts in the table be offset by before taking the log (to avoid taking log(0))? Default 32.
## --nonzero:  if TRUE, library size adjustment is median of nonzero counts.  If FALSE, library size adjustment is median of all counts.
## --colmeds: vector with the column medians or NULL if you want to calculate them this time around. To save computing time, you might want to calculate them once with getColmeds().
## return:
## a list containing elements $ebobject, which is a list mimicking the output of running lmFit on the whole dataset,
## and $pos, giving the row indices of all the positions on the chromosome passing the filtering criterion (i.e. bases with nontrivial coverage)
## usage recommendataion: run getLimmaInput(), pass $ebobject to getTstats(), save $pos as rda file in case it is needed later





#'Fit a linear model to each nucleotide
#'
#'Fits the linear model log2(count+scalefac) = beta0 + beta1*group +
#'beta2*library.size + [optional confounders] to each nucleotide.  From these
#'models, this function constructs an object which can be directly passed to
#'\code{getTstats} to obtain \code{limma}'s moderated t statistics for each
#'nucleotide, which we use as a measure of strength of association between
#'group and count (expression). Reads coverage file from a SQLite database (see
#'\code{makeDb}) and relies heavily on the \code{limma} package, using
#'\code{lmFit} as the main workhorse.
#'
#'It is assumed that the first column in the database is called \code{pos} and
#'contains genomic position.  Note that "group" must have the one fewer entries
#'than the database denoted by \code{dbfile} has columns (or, if colsubset is
#'used, one fewer entries than the length of colsubset). Also,
#'\code{adjustvars} must have the same number of rows as \code{group} has
#'entries. Larger values of \code{chunksize} require more memory; smaller
#'values of \code{chunksize} require more computation time.
#'
#'@param dbfile Name/location (as character string) of database (usually ".db")
#'file containing nucleotide by sample coverage.
#'@param tablename Name of the table the database contains
#'@param comparison Either \code{twogroup}, \code{multigroup} or \code{expression}. \code{multigroup} will use the F-statistic and \code{expression} tests the intercept-only model.
#'@param group a 0/1 vector grouping the samples (columns) in the database.
#'@param chunksize How many rows of the database should be processed at a time?
#'@param adjustvars Optional matrix of adjustment variables (e.g. measured
#'confounders, output from SVA, etc.) to use in fitting linear models to each
#'nucleotide.
#'@param colsubset Optional vector of column indices of the input file that
#'denote samples you wish to include in analysis. Should NEVER include 1
#'(genomic position).
#'@param scalefac A log transformation is used on the count tables, so zero
#'counts present a problem.  What number should we add to the entire matrix
#'before running the models?  Defaults to 32.
#'@param nonzero If TRUE, use the median of only the nonzero counts as the
#'library size adjustment.
#'@param colmeds If NULL, the column medians are calculated using
#'\code{\link{getColmeds}}. Otherwise, the output of \code{\link{getColmeds}}
#'is expected.
#'@return A list with elements
#'\item{ebobject }{A list of five vectors (\code{coefficients}, \code{stdev.unscaled}, \code{sigma}, \code{df.residual}, and \code{Amean}), mimicking the \code{MArrayLM} class in \code{limma}.  Here, \code{coefficients} and \code{stdev.unscaled} are only returned for beta1, the coefficient for \code{group}, as it is assumed this is the only covariate of interest.}
#'\item{pos }{A vector of the same length as those contained in \code{ebobject}, giving the genomic positions of each linear model.}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getTstats}}, \code{\link{makeDb}}, \code{\link{lmFit}},
#'\code{\link{MArrayLM-class}}
#'@references Smyth G (2004).  “Linear models and empirical Bayes methods for
#'assessing differential expression in microarray experiments.” Statistical
#'Applications in Genetics and Molecular Biology 3(1): Article 3.
#'@examples
#'
#'## add example here when we have a vignette
#'
getLimmaInput <- function(dbfile, tablename, comparison = c("twogroup", "multigroup", "expression"), group = NULL, chunksize = 100000, adjustvars = NULL, colsubset = c(-1), scalefac = 32, nonzero = FALSE, colmeds=NULL){
	require(limma)
	require(multicore)
	require(Genominator)
	
	# get ready to read from database:
	tab = ExpData(dbfile, tablename)
	pos = tab[,1]$pos
	N = length(pos)
	lastloop = trunc(N/chunksize)

	# create model matrix:
	numcol = length(tab[1,])
		
	## Get the medians of the columns of interest if they are not given, and if you are doing differential expression:
	if(is.null(colmeds) & comparison!="expression"){
		colmeds <- getColmeds(dbfile=dbfile, tablename=tablename, colsubset=colsubset, nonzero=nonzero)
	}
	
  if(!is.null(adjustvars) & comparison=="expression") stop("expression-only test does not require adjustment variables - please set adjustvars = NULL when using comparison = \"expression\"")
  
	if(!is.null(adjustvars)){
		string1 = ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
			string1 = paste(string1, paste("av",i,sep=""),sep="+")
		}
		eval(parse(text=paste("x = model.matrix(~group+colmeds",string1,")",sep="")))
	} else {
    if(comparison == "expression"){ 
      numsamps = ifelse(is.null(colsubset), numcol-1, length(colsubset))
      int = rep(1, numsamps)
      x = model.matrix(~0+int) #intercept-only model, in model.matrix form.
    }
		if(comparison == "twogroup" | comparison == "multigroup"){
      x = model.matrix(~as.factor(group)+colmeds)
		}
	}	
	
	# make sure colsubset will work in the next part:
	if(colsubset[1]==c(-1) & length(colsubset) == 1) colsubset <- c(2:numcol)
	
	# define modeling function to apply to each chunk:
	lmFit.apply = function(i){
  		if(i!=lastloop) mymat <- tab[(chunksize * i + 1):(chunksize * (i + 1)), -1][, colsubset - 1] 
  		else mymat <- tab[(chunksize * i + 1):N, -1][, colsubset - 1] 
  		mymat <- log2(mymat + scalefac)
  		Amean <- rowMeans(mymat) 
  		fit <- lmFit(mymat, x)
  		return(list(fit=fit, Amean=Amean))
  		}
	
	# fit a model to each row (chunk) of database:
	if(interactive()) {
		warning("Cannot use mclapply in interactive session; reverting to single-core lapply")
		lmFit.output = lapply(0:lastloop, lmFit.apply)
		}
	if(!interactive()) {
		warning("mclapply functionality not yet implemented: fitting models with one core.")
		lmFit.output = lapply(0:lastloop, lmFit.apply)
	}
	# gather results from all the models and return them:
	coef = stdev = sma = dfres = am = NULL
	for(i in 1:length(lmFit.output)){
  		coef = append(coef,lmFit.output[[i]]$fit$coefficients[,2])
  		stdev = append(stdev, lmFit.output[[i]]$fit$stdev.unscaled[,2])
  		sma = append(sma, lmFit.output[[i]]$fit$sigma)
  		dfres = append(dfres, lmFit.output[[i]]$fit$df.residual)
  		am = append(am, lmFit.output[[i]]$Amean)
	}
	return(list(ebobject = list(coefficients = as.numeric(coef), stdev.unscaled = as.numeric(stdev), sigma = sma, df.residual = dfres, Amean = am), pos = pos, comparison = comparison))
	
}
