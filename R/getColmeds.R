## This function gets the column medians
## Ideally, you won't need to calculate them more than once for the whole analysis.

getColmeds <- function(dbfile, tablename, colsubset=c(-1), nonzero=FALSE) {
	## Load libraries
	library(Genominator)
	
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