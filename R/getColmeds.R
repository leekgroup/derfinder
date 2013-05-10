## This function gets the column medians
## Ideally, you won't need to calculate them more than once for the whole analysis.



#'Calculate the medians of the database.
#'
#'\code{getColmeds} calculates the median of the samples of interest from the
#'database table. It is used in \code{\link{getLimmaInput}} and it's output can
#'be saved to avoid calculating it multiple times when running
#'\code{\link{get.pvals}}.
#'
#'It is assumed that the first column in the database is called \code{pos} and
#'contains genomic position.  Note that "group" must have the one fewer entries
#'than the database denoted by \code{dbfile} has columns (or, if colsubset is
#'used, one fewer entries than the length of colsubset).
#'
#'@param dbfile Name/location (as character string) of database (usually ".db")
#'file containing nucleotide by sample coverage.
#'@param tablename Name of the table the database contains
#'@param colsubset Optional vector of column indices of the input file that
#'denote samples you wish to include in analysis. Should NEVER include 1
#'(genomic position).
#'@param nonzero If TRUE, use the median of only the nonzero counts as the
#'library size adjustment.
#'@return A vector with the column medians.
#'@author Alyssa Frazee. Code split from getLimmaInput by Leonardo
#'Collado-Torres.
#'@seealso \code{\link{getLimmaInput}}
#'@export
#'@examples
#'
#'## add example here when we have a vignette
#'
getColmeds <- function(dbfile, tablename, colsubset=c(-1), nonzero=FALSE) {
	## Load libraries
	require(Genominator)
	
	## Define the database to use
	tab <- ExpData(dbfile, tablename)
	
	## Get the number of columns
	ncol  <- length(tab[1,])
		
	## Default behavior
	colsToCheck <- 2:ncol
	## Adjust the behavior if the user supplies colsubset manually
	if(length(colsubset) > 1) colsToCheck <- colsubset
		
	## Calculate the medians
	colmeds <- NULL
	for(i in colsToCheck){
		if(nonzero){
			eval(parse(text=paste0("med = median(tab[,", i, "]$", names(tab[, i]), "[tab[,", i, "]$", names(tab[, i]), "!=0]", ")")))
		}
		if(!nonzero){
			eval(parse(text=paste0("med = median(tab[,", i, "]$", names(tab[,i]), ")")))
		}
		colmeds <- c(colmeds, med)
	} #get median of each column to use as adjustment variable
	
	## Finish
	return(colmeds)
}
