### written by JTL
### argument descriptions: see getLimmaInput
### return: object of class sva
### this is used to get surrogate variables (in order to adjust for unmeasured confounders)

getSvs <- function(dbfile, tablename, group, chunksize = 100000,colsubset = c(-1)){
  require(limma)
  require(multicore)
  require(Genominator)
  require(sva)
  require(genefilter)
	
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
