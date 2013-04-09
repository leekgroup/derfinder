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



getLimmaInput =
function(dbfile, tablename, group, chunksize = 100000, adjustvars = NULL, colsubset = c(-1), scalefac = 32, nonzero = FALSE, colmeds=NULL){
	require(limma)
	require(multicore)
	require(Genominator)
	
	# get ready to read from database:
	tab = ExpData(dbfile, tablename)
	pos = tab[,1]$pos
	N = length(pos)
	lastloop = trunc(N/chunksize)

	# create model matrix:
	ncol = length(tab[1,])
		
	## Get the medians of the columns of interest if they are not given
	if(is.null(colmeds)){
		colmeds <- getColmeds(dbfile=dbfile, tablename=tablename, colsubset=colsubset, nonzero=nonzero)
	}
	
	if(!is.null(adjustvars)){
		string1 = ""
		for(i in 1:dim(adjustvars)[2]){
			eval(parse(text=paste("av",i," <- adjustvars[,",i,"]",sep="")))
			string1 = paste(string1, paste("av",i,sep=""),sep="+")
		}
		eval(parse(text=paste("x = model.matrix(~group+colmeds",string1,")",sep="")))
	} else {
		x = model.matrix(~group+colmeds)
	}	
	
	# make sure colsubset will work in the next part:
	if(colsubset[1]==c(-1) & length(colsubset) == 1) colsubset <- c(2:ncol)
	
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
	return(list(ebobject = list(coefficients = as.numeric(coef), stdev.unscaled = as.numeric(stdev), sigma = sma, df.residual = dfres, Amean = am), pos = pos))
	
}
