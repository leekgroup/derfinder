### written by JTL
### argument descriptions: see getLimmaInput
### return: object of class sva
### this is used to get surrogate variables (in order to adjust for unmeasured confounders)



#'Do SVA on matrix stored in database
#'
#'Performs surrogate variable analysis on a matrix, usually base-pair by sample
#'coverage, so the surrogate variables can be adjusted for in later analysis
#'steps.  SVA is meant to provide a method to adjust for unknown confounders.
#'

#'
#'@param dbfile Name/location (as character string) of database (usually ".db")
#'file containing nucleotide by sample coverage.
#'@param tablename Name of the table the database contains
#'@param group a 0/1 vector grouping the samples (columns) in the database.
#'@param chunksize How many rows of the database should be processed at a time?
#'@param colsubset Optional vector of column indices of the input file that
#'denote samples you wish to include in analysis. Should NEVER include 1
#'(genomic position).
#'@return An object of class sva - see help files for sva.
#'@author Jeff Leek
#'@export

getSvs <- function(dbfile, tablename, group, chunksize = 100000,colsubset = c(-1)){
  require("limma")
  require("multicore")
  require("Genominator")
  require("sva")
  require("genefilter")
	
  tab = ExpData(dbfile, tablename)
  pos = tab[,1]$pos
  print("here")

  N = length(pos)
  lastloop = trunc(N/chunksize)
  print(lastloop)
  
  vv.apply = function(i){
    if(i!=lastloop) mymat <- tab[(chunksize*i+1):(chunksize*(i+1)),colsubset]
    else mymat <- tab[(chunksize*i+1):N,colsubset]
    mymat <- log2(as.matrix(mymat)+0.5)
    rowV <- rowVars(mymat)
    return(list(rowV=rowV))
  }

  vv.output = lapply(0:lastloop,vv.apply)

  # Collect all the variances
  vv = NULL
  for(i in 1:length(vv.output)){
    vv = append(vv,vv.output[[i]]$rowV)  
  }

  # Set up the model matrices
  ind = which(rank(-vv) <= 10000)
  mod = model.matrix(~as.factor(group))
  mod0 = cbind(mod[,1])

  # Run SVA
  svaObj = sva(as.matrix(tab[ind,colsubset]),mod=mod,mod0=mod0)
  
  return(svaObj)
}
