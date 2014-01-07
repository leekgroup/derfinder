# CODE FOR DERFINDER PAPER
# uses beta version of derfinder package

#########################################
#########################################
##### chunk 1: male vs. female analysis #
#########################################
#########################################

print(Sys.time())
library(devtools)
install_github('derfinder', 'alyssafrazee') # beta version
library(derfinder)

dbfile = "Y-tophat-revised.db"
tablename = "chrY"
textfile = "tophatY-updated"

# create database:
makeDb(dbfile = dbfile, tablename = tablename, textfile = textfile, cutoff = 5)

sex = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1) #based on the samples in textfile
limma.input = getLimmaInput(dbfile = dbfile, tablename = tablename, group = sex, nonzero = TRUE)
pos = limma.input$pos
save(pos, file="posY-rev.rda")

# get the moderated t stats and fold changes:
tstats = getTstats(fit = limma.input$ebobject, trend = TRUE)
tt = tstats$tt
logfchange = tstats$logfchange
save(tt,file="ttY-rev.rda")

# fit the HMM:
find.them = getParams(tt)
regions = getRegions(method = "HMM", chromosome = "Y", pos = pos, tstats = tt, stateprobs = find.them$stateprobs, params = find.them$params, includet = TRUE, includefchange = TRUE, fchange = logfchange)

# merge the regions:
regions.merged.y = mergeRegions(regions$states)
save("regions.merged.y",file="Ychrom-regions-merged-new-rev.rda")

# get the p-values:
##### define p-value function to print status messages & the total number of null statistics, which we needed for one of the reviewer responses.
get.pvals <- function (regions, dbfile, tablename, num.perms = 1, group, est.params, chromosome, colsubset = c(-1)){
  nullstats = NULL
  for (i in 1:num.perms) {
    group.permute = sample(group)
    print(paste("starting iteration",i))
    limma.input = getLimmaInput(dbfile = dbfile, tablename = tablename,group = group.permute, colsubset = colsubset, nonzero = TRUE)
    tstats = getTstats(fit = limma.input$ebobject, trend = TRUE)
    tt = tstats$tt
    logfchange = tstats$logfchange
    regions.null = getRegions(method = "HMM", chromosome = chromosome,
        pos = limma.input$pos, tstats = tt, stateprobs = est.params$stateprobs,
        params = est.params$params, includet = TRUE, includefchange = TRUE,
        fchange = logfchange)
    nullstats = append(nullstats, regions.null$states$mean.t[regions.null$states$state == 3 | regions.null$states$state == 4])
  }
  print("Number of null stats:")
  print(length(nullstats))
  pvals = rep(NA, dim(regions)[1])
  for (k in which(regions$state == 3 | regions$state == 4)) {
    pvals[k] = (sum(abs(nullstats) > abs(regions$mean.t[k]))+1)/(length(nullstats)+1)
  }
  return(pvals)
}


pvals = get.pvals(regions = regions.merged.y, dbfile = dbfile, tablename = tablename, num.perms = 10, group = sex, est.params = find.them, chromosome = "Y")
save(pvals, file="Y-pvals-new-rev.rda")

# get the flags:
exons = getAnnotation("hg19","knownGene")
myflags = getFlags(regions = regions.merged.y, exons, "chrY", pctcut = 0.8)
save(myflags, file="Y-flags-new-rev.rda")









#########################################
#########################################
##### chunk 2: male vs. male analysis ###
#########################################
#########################################

dbfile = "Y-tophat-revised.db"
tablename = "chrY"
textfile = "tophatY-updated"

# create database: (ALREADY MADE IN CHUNK 1)
#makeDb(dbfile = dbfile, tablename = tablename, textfile = textfile, cutoff = 5)

to.use = c(2,3,6,7,10,11,12,13,16)
set.seed(651)
sex = sample(c(rep(1,5),rep(0,4))) #randomly assign groups of males
limma.input = getLimmaInput(dbfile = dbfile, tablename = tablename, group = sex, colsubset = to.use, nonzero = TRUE)
pos = limma.input$pos

# get the moderated t stats and fold changes:
tstats = getTstats(fit = limma.input$ebobject, trend = TRUE)
tt = tstats$tt
logfchange = tstats$logfchange
save(tt,file="ttY-men-rev.rda")

# fit the HMM:
find.them = getParams(tt)
regions = getRegions(method = "HMM", chromosome = "Y", pos = pos, tstats = tt, stateprobs = find.them$stateprobs, params = find.them$params, includet = TRUE, includefchange = TRUE, fchange = logfchange)

# merge the regions:
regions.merged.y = mergeRegions(regions$states)
save("regions.merged.y",file="Ychrom-regions-merged-men-rev.rda")

# get the p-values: (uses same "debugging get.pvals" as above)
pvals = get.pvals(regions = regions.merged.y, dbfile =dbfile,tablename = tablename, num.perms = 10, group = sex, est.params = find.them, chromosome = "Y", colsubset = to.use)
save(pvals, file="Y-pvals-men-rev.rda")



#########################################
#########################################
##### chunk 3: analyze results ##########
#########################################
#########################################

## load in the exon annotation (ENSEMBL GRCh37)
library(biomaRt)
ensembl = useMart("ensembl")
listDatasets(ensembl)[which(listDatasets(ensembl)$dataset == "hsapiens_gene_ensembl"),]
### VERSION = GRCh37.p12
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
ensexons = getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id", 
  "ensembl_exon_id", "chromosome_name", "exon_chrom_start", 
  "exon_chrom_end", "strand"), mart=ensembl)
# ^takes a few minutes
names(ensexons)[3] = "exon_id"
names(ensexons)[1] = "geneName"
names(ensexons)[4] = "chr"
names(ensexons)[5] = "start"
names(ensexons)[6] = "end"
ensYexons = subset(ensexons, chr=="Y")
nrow(ensYexons) #3748 exons, some duplicated.
ensYexons = ensYexons[!duplicated(ensYexons$exon_id),]
save(ensYexons, file="ensYexons.rda")


#################################################
##### chunk 3a: load resuls from DER Finder #####
#################################################


## load results from our pipeline:
load("ensYexons.rda")
load("Ychrom-regions-merged-men-rev.rda") #regions.merged.y
regions.merged.men = regions.merged.y
load("Ychrom-regions-merged-new-rev.rda") #also regions.merged.y
load("Y-pvals-men-rev.rda") #pvals
pvals.men = pvals
load("Y-pvals-new-rev.rda") #also pvals
regions = data.frame(regions.merged.y,pvals,qvals = p.adjust(pvals, method="fdr"))
regions.men = data.frame(regions.merged.men,pvals = pvals.men, qvals = p.adjust(pvals.men, method="fdr"))

### P VALUE HISTOGRAMS - check
hist(pvals,breaks=30,col="gray70",main="Y chromosome p values: our method",xlab="p values") #yes.
hist(pvals.men,breaks=30,col="gray70",main="Y chromosome p values: our method, men only",xlab="p values") 

### get data frame of only DE regions (with p and q-values)
ders = subset(regions,state==3|state==4)
ders.men = subset(regions.men, state==3|state==4)

# re-create myflags with ensembl exons
myflags = getFlags(regions = regions.merged.y, ensYexons, "Y", pctcut = 0.8)
myflags.men = getFlags(regions = regions.merged.men, ensYexons, "Y", pctcut = 0.8)
save(myflags, file="Y-flags-ensembl.rda")
save(myflags.men, file="Y-flags-men-ensembl.rda")

# if flags already created: 
# load("Y-flags-ensembl.rda")
# load("Y-flags-men-ensembl.rda")

ders = data.frame(ders,flag = myflags$flags)
ders.men = data.frame(ders.men,flag=myflags.men$flags)

# line up exons with p-values (i.e., get rid of the list situation)
info.lengths = sapply(myflags$flag.info, length, USE.NAMES=FALSE)
qvals.ex = ex.names = ex.pct = ex.class = list()
for(i in 1:length(myflags$flags)){
  qvals.ex[[i]] = rep(ders$qvals[i],info.lengths[i])
  ex.names[[i]] = myflags$flag.info[[i]]
  ex.pct[[i]] = myflags$percent.exon[[i]]
  ex.class[[i]] = rep(myflags$flags[i],info.lengths[i])
}
qvals.ex = unlist(qvals.ex)
ex.names = unlist(ex.names)
ex.class = unlist(ex.class)
ex.pct = unlist(ex.pct)

myinfo = data.frame(qvals.ex, ex.names, ex.class, ex.pct, stringsAsFactors=FALSE)
write.table(myinfo, file="myinfo_ensembl.txt",row.names=F,quote=F,sep="\t")

# fix places where exons are covered by >1 different region
#### RULES: if different parts of the same exon were overlapped 
#### by different DE regions, that exon was reduced to just one 
#### row, and the region was then flagged as DE if those combined 
#### regions overlapped the exon by 80%.  The q-value of 
#### the largest region making up that exon's covering was taken
#### as the q-value for that exon. 
#### (So if an exon was at least 80% overlapped by DERs 
#### with q<0.05, it was called DE at the 0.05 level.)  
rep_exons = names(table(myinfo$ex.names))[table(myinfo$ex.names)>1]
repsub = subset(myinfo, ex.names %in% rep_exons)
percentsplit = split(repsub$ex.pct, repsub$ex.names)
qsplit = split(repsub$qvals.ex, repsub$ex.names)
sum(names(qsplit) != rep_exons) ## 0; GOOD
newpct = sapply(percentsplit, sum, USE.NAMES=FALSE)
newq = sapply(1:length(qsplit), function(i){
    qsplit[[i]][which.max(percentsplit[[i]])]
}, USE.NAMES=FALSE)
newclass = ifelse(newpct>=0.8, "DE exons(s)", NA)
myinfo.updated = subset(myinfo, !(ex.names %in% rep_exons))
myinfo.updated = rbind(myinfo.updated, data.frame(qvals.ex = newq, 
    ex.names=names(qsplit), ex.class=newclass, ex.pct=newpct))
rownames(myinfo.updated) = NULL

myinfo.sub = subset(myinfo.updated, ex.class!="novel")
nrow(myinfo.sub) #467 exons

checkunique = function(x){
  sum(ensYexons$exon_id==x)
}
checkexons = unlist(lapply(myinfo.updated$ex.names, checkunique))
table(checkexons) #no id issues
foo = which(unlist(lapply(myinfo.updated$ex.names, checkunique))>1)

# same thing for men:
info.lengths.m = sapply(myflags.men$flag.info, length, USE.NAMES=FALSE)
qvals.ex.m = ex.names.m = ex.pct.m = ex.class.m = list()
for(i in 1:length(myflags.men$flags)){
  qvals.ex.m[[i]] = rep(ders.men$qvals[i],info.lengths.m[i])
  ex.names.m[[i]] = myflags.men$flag.info[[i]]
  ex.pct.m[[i]] = myflags.men$percent.exon[[i]]
  ex.class.m[[i]] = rep(myflags.men$flags[i],info.lengths.m[i])
}
qvals.ex.m = unlist(qvals.ex.m)
ex.names.m = unlist(ex.names.m)
ex.class.m = unlist(ex.class.m)
ex.pct.m = unlist(ex.pct.m)
myinfo.m = data.frame(qvals.ex.m, ex.names.m, ex.class.m, ex.pct.m)
min(myinfo.m$qvals.ex.m) #the min q-value here is 0.86, so nothing will show up in the table. 




#################################################
##### chunk 3a: load resuls from Cufflinks ######
#################################################

# load 'isoform.diff' files from Cuffdiff 
load("tx.sex.cuff-UPDATED.rda")
load("tx.null.cuff-UPDATED.rda")
tx.sex.cuff.y.ok = subset(tx.sex.cuff,substr(locus,1,1)=="Y" & status=="OK")
## recalibrate q values (only using Y chromosome)
tx.sex.cuff.y.ok$qvalue = p.adjust(tx.sex.cuff.y.ok$pvalue,method="fdr")
cuff.t = tx.sex.cuff.y.ok$test_stat
# fix some very, very large (10e308) test stats...
cuff.t = ifelse(cuff.t>20,20,cuff.t)
cuff.t = ifelse(cuff.t<(-20),-20,cuff.t)
cuff.fchange = tx.sex.cuff.y.ok$log2fchange
cuff.fchange = ifelse(cuff.fchange>20,20,cuff.fchange)
cuff.fchange = ifelse(cuff.fchange<(-20),-20,cuff.fchange)
cuff.t[which(cuff.fchange==-20)] = 20

tx.null.cuff.y.ok = subset(tx.null.cuff,substr(locus,1,1)=="Y" & status=="OK")
# recalibrate q:
tx.null.cuff.y.ok$qvalue = p.adjust(tx.null.cuff.y.ok$pvalue,method="fdr")

tx.sex.cuff.y = subset(tx.sex.cuff,substr(locus,1,1)=="Y")
tx.null.cuff.y = subset(tx.null.cuff,substr(locus,1,1)=="Y")

### P VALUE HISTOGRAMS ### - just to check
hist(tx.sex.cuff.y.ok$pvalue, col="gray70",breaks=30, xlab="p values",main="Y chromosome p values: Cufflinks")
hist(tx.null.cuff.y.ok$pvalue, col="gray70",breaks=30, xlab="p values",main="Y chromosome null p values: Cufflinks")

# merged.gtf results, containing the exons of each transcript:
# (result of Cuffmerge)
# (this is for the supplemental table)
merged.gtf = read.table("normal-comparisons/merged.gtf",sep="\t") 
merged.gtf$V9 = as.character(merged.gtf$V9)
summary(merged.gtf$V3) #all exon- good.
elts = unlist(strsplit(merged.gtf$V9,split="; "))
isstandard = function(string){
  isgid = substring(string,1,7)=="gene_id"
  istx = substring(string,1,13)=="transcript_id"
  isexnum = substring(string,1,11)=="exon_number"
  isoid = substring(string,1,3)=="oId"
  istssid = substring(string,1,6)=="tss_id"
  return(isgid|istx|isexnum|isoid|istssid)
}
elts = elts[-which(!isstandard(elts))] #this should work...
merged.gtf.new = merged.gtf[,-9]
names(merged.gtf.new) = c("chr","software","feature","start","end","score","strand","frame")

striplabel = function(x){
  splitme = strsplit(x,split=" ")[[1]][2]
  return(splitme)
}

elts = unlist(lapply(elts,striplabel))
eltsdf = matrix(elts,ncol=5,byrow=TRUE)
merged.gtf.new = data.frame(chr=merged.gtf.new$chr,start=merged.gtf.new$start, end=merged.gtf.new$end, feature=merged.gtf.new$feature, geneid=eltsdf[,1], txid=eltsdf[,2], exnum = as.numeric(eltsdf[,3])) 
merged.gtf.y = subset(merged.gtf.new, chr=="Y")
merged.gtf.y = merged.gtf.y[order(merged.gtf.y$start),] # we want to use THIS - just add in qvals & test status
qvalue = status = NULL

for(i in 1:dim(merged.gtf.y)[1]){
  print(i)
  ind = which(tx.sex.cuff.y$tx_id==merged.gtf.y$txid[i])
  status[i] = as.character(tx.sex.cuff.y$status[ind])
  if(status[i]=="OK"){
    ind2 = which(tx.sex.cuff.y.ok$tx_id==merged.gtf.y$txid[i])
    qvalue[i] = tx.sex.cuff.y.ok$qvalue[ind2]
  }
  if(status[i]!="OK") qvalue[i] = 1
} #takes a while (~ 1 hour)

merged.gtf.y = data.frame(merged.gtf.y, qvalue, status) #sweet.
save(merged.gtf.y,file="cuffdiff-sex/merged.gtf.y.rda")

### IF ABOVE ALREADY DONE:
load("cuffdiff-sex/merged.gtf.y.rda")
merged.gtf.y.ok = subset(merged.gtf.y,status=="OK")
unique.transcripts = merged.gtf.y.ok[!duplicated(as.character(merged.gtf.y.ok$txid)),] 
dim(unique.transcripts) #808 rows

### now get the merged file for the men only: (repeat previous process)
merged.gtf.men = read.table("merged-men.gtf",sep="\t") 
merged.gtf.men$V9 = as.character(merged.gtf.men$V9)
summary(merged.gtf.men$V3) #all exon- good.

elts.men = unlist(strsplit(merged.gtf.men$V9,split="; "))
elts.men = elts.men[-which(!isstandard(elts.men))] #this should work...
merged.gtf.men.new = merged.gtf.men[,-9]
names(merged.gtf.men.new) = c("chr","software","feature","start","end","score","strand","frame")
elts.men = unlist(lapply(elts.men,striplabel))
eltsdf.men = matrix(elts.men,ncol=5,byrow=TRUE)
merged.gtf.men.new = data.frame(chr=merged.gtf.men.new$chr,start=merged.gtf.men.new$start, end=merged.gtf.men.new$end, feature=merged.gtf.men.new$feature, geneid=eltsdf.men[,1], txid=eltsdf.men[,2], exnum = as.numeric(eltsdf.men[,3])) 
merged.gtf.men.y = subset(merged.gtf.men.new, chr=="Y")
merged.gtf.men.y = merged.gtf.men.y[order(merged.gtf.men.y$start),] 
qvalue.m = status.m = NULL
dim(merged.gtf.men.y) #1961 rows

for(i in 1:dim(merged.gtf.men.y)[1]){
  print(i)
  ind = which(tx.null.cuff.y$tx_id==merged.gtf.men.y$txid[i])
  status.m[i] = as.character(tx.null.cuff.y$status[ind])
  if(status.m[i]=="OK"){
    ind2 = which(tx.null.cuff.y.ok$tx_id==merged.gtf.men.y$txid[i])
    qvalue.m[i] = tx.null.cuff.y.ok$qvalue[ind2]
  }
  if(status.m[i]!="OK") qvalue.m[i] = 1
} #again, takes a while (maybe an hour)


merged.gtf.men.y = data.frame(merged.gtf.men.y, qvalue=qvalue.m, status=status.m) #sweet.
save(merged.gtf.men.y,file="cuffdiff-males/merged.gtf.men.y.rda")

##### IF ABOVE ALREADY DONE:
load("cuffdiff-males/merged.gtf.men.y.rda")

merged.gtf.men.y.ok = subset(merged.gtf.men.y,status=="OK")
unique.transcripts.men = merged.gtf.men.y.ok[!duplicated(as.character(merged.gtf.men.y.ok$txid)),] 
dim(unique.transcripts.men) #818 rows


## how many of the cufflinks (or other) DE transcripts are 80% covered by DERs from our pipeline? 
doweagree = function(qcut, unique.transcripts, ders, txinfo, pctcut = 0.8){
  #unique.transcripts is a data frame containing txid and q-value (one per transcript)
  #txinfo is the merged.gtf file, containing the exon locations and txids. can have multiple rows per transcript, but only one row per exon
  k=0
  percent.covered = NULL
  dercands = subset(ders,qvals<=qcut)
  for(i in which(unique.transcripts$qvalue<qcut)){
    print("hey")
    print(i)
    k=k+1
    tname = as.character(unique.transcripts$txid[i])
    theexons = subset(txinfo, txid==tname)
    firstderind = findInterval(theexons$start[1], dercands$start)
    lastderind = findInterval(theexons$end[dim(theexons)[1]], dercands$start)
    if(firstderind==0) firstderind = 1
    if(lastderind==0){percent.covered[k] = 0; next}
    these.ders = dercands[firstderind:lastderind,]
    transcript.pos = NULL
    for(j in 1:dim(theexons)[1]){
      transcript.pos = append(transcript.pos, c(theexons$start[j]:theexons$end[j]))
    }
    ders.pos = NULL
    for(j in 1:dim(these.ders)[1]){
      ders.pos = append(ders.pos, c(these.ders$start[j]:these.ders$end[j]))
    }
    percent.covered[k] = sum(transcript.pos %in% ders.pos)/length(transcript.pos)
  }
  
  numcovered = sum(percent.covered>pctcut)
  return(numcovered)
}

## Among our regions, how many are 80% covered by a DE cufflinks transcript?
dotheyagree = function(qcut, ychr.cuff, ders){
  #ychr.cuff is the merged.gtf file, with q-values and all exons, from cufflinks
  #ders is our ders data frame.
  require(IRanges)
  cover80pct = NULL
  k = 0
  ychr.cuff.sig = subset(ychr.cuff,qvalue<qcut)
  for(i in which(ders$qvals<qcut)){
    print(i)
    k = k+1
    cuffind = findInterval(ders$start[i], ychr.cuff.sig$start)
    endind = findInterval(ders$end[i], ychr.cuff.sig$start)
    if(endind==0){
      cover80pct[k] = 0
      next
    }
    
    endInendind = ders$end[i] %in% c((ychr.cuff.sig$start[endind]):(ychr.cuff.sig$end[endind]))
    
    if(cuffind==0){
      x = IRanges(start = ychr.cuff$start[1:endind], end = ychr.cuff$end[1:endind])
      x = reduce(x)
      tosum = width(x)
      if(endInendind){tosum[length(tosum)] = ders$end[i] - start(x)[length(tosum)] + 1}
      numcovered = sum(tosum)
      percent.covered = numcovered/(ders$length[i])
      cover80pct[k] = ifelse(percent.covered>=0.8,1,0)
      next
    } #end cuffind==0 case
    
    startIncuffind = ders$start[i] %in% c((ychr.cuff.sig$start[cuffind]):(ychr.cuff.sig$end[cuffind]))
    
    if(cuffind==endind){
      if(startIncuffind & endInendind){
        cover80pct[k] = 1
        next
      }
      
      if(startIncuffind & !endInendind){
        numcovered = ychr.cuff$end[cuffind] - ders$start[i] + 1
        percent.covered = numcovered/(ders$length[i])
        cover80pct[k] = ifelse(percent.covered>=0.8,1,0)
        next
      }
      
      if(!startIncuffind){
        cover80pct[k] = 0
        next
      }
    } #end equal case
    
    
    if(cuffind!=endind){
      x = IRanges(start = ychr.cuff$start[cuffind:endind], end = ychr.cuff$end[cuffind:endind])
      x = reduce(x)
      tosum = width(x)
      if(endInendind){tosum[length(tosum)] = ders$end[i] - start(x)[length(tosum)] + 1}
      if(startIncuffind){tosum[1] = end(x)[1] - ders$start[i] + 1}
      if(!startIncuffind){tosum[1] = 0}
      numcovered = sum(tosum)
      percent.covered = numcovered/(ders$length[i])
      cover80pct[k] = ifelse(percent.covered>=0.8,1,0)
    } #end unequal case
    
  } #end for i
  return(sum(cover80pct))
} #end function



#### CREATE TABLE comparing our results to cufflinks results

cufftable = matrix(NA, ncol=9, nrow=17)
cufftable = as.data.frame(cufftable)
names(cufftable) = c("q","numregions","numtranscripts","doweagree","dotheyagree","numregions.men","numtranscripts.men","doweagree.men","dotheyagree.men")
cufftable$q = seq(0,0.8,by=0.05)

for(r in 1:17){
  cufftable$numregions[r] = length(which(ders$qvals<cufftable$q[r])) 
  cufftable$numtranscripts[r] = length(which(unique.transcripts$qvalue<cufftable$q[r]))
  if(cufftable$numregions[r]==0 & cufftable$numtranscripts[r]==0){
    cufftable$dotheyagree[r] = NA
    cufftable$doweagree[r] = NA     
  }else if(cufftable$numtranscripts[r]==0){
    cufftable$dotheyagree[r] = 0
    cufftable$doweagree[r] = NA
  }else if(cufftable$numregions[r]==0){
    cufftable$doweagree[r] = 0
    cufftable$dotheyagree[r] = NA       
  }else{
    cufftable$dotheyagree[r] = dotheyagree(cufftable$q[r], ychr.cuff = merged.gtf.y.ok, ders = ders)
    cufftable$doweagree[r] = doweagree(cufftable$q[r], unique.transcripts = unique.transcripts, txinfo = merged.gtf.y.ok, ders = ders)
  }
  
  cufftable$numregions.men[r] = length(which(ders.men$qvals<cufftable$q[r])) 
  cufftable$numtranscripts.men[r] = length(which(unique.transcripts.men$qvalue<cufftable$q[r]))
  if(cufftable$numregions.men[r]==0 & cufftable$numtranscripts.men[r]==0){
    cufftable$dotheyagree.men[r] = NA
    cufftable$doweagree.men[r] = NA     
  }else if(cufftable$numtranscripts.men[r]==0){
    cufftable$dotheyagree.men[r] = 0
    cufftable$doweagree.men[r] = NA
  }else if(cufftable$numregions.men[r]==0){
    cufftable$doweagree.men[r] = 0
    cufftable$dotheyagree.men[r] = NA
  }else{
    cufftable$dotheyagree.men[r] = dotheyagree(cufftable$q[r], ychr.cuff = merged.gtf.y.men.ok, ders = ders.men)
    cufftable$doweagree[r] = doweagree(cufftable$q[r], unique.transcripts = unique.transcripts.men, txinf = merged.gtf.men.y.ok, ders = ders.men)
  }
}

cufftable   ### supplementary table 2
save(cufftable, file="cufftable-correct-rev.rda")

## if above already done:
load("cufftable-correct-rev.rda")
cufftable

r=11 #q=0.5
doweagree(cufftable$q[r], unique.transcripts = unique.transcripts, txinfo = merged.gtf.y.ok, ders = ders, pctcut = 0.25) #99
sum(unique.transcripts$qvalue<=0.5) #758
99/758 #13.1 %

doweagree(cufftable$q[r], unique.transcripts = unique.transcripts, txinfo = merged.gtf.y.ok, ders = ders, pctcut = 0) #135
135/758 #17.8%



#####################################################
##### chunk 3b: get resuls from DESeq & EdgeR #######
#####################################################

###############################################
###### chunk 3b.1 -- make the exon count tables
###############################################

library(Rsamtools)
library(GenomicFeatures)
load("ensYexons.rda")

# take out same exon annotated in different transcripts:
yx = unique(ensYexons[,-2]) 

exons = GRanges(seqnames = Rle("Y"), 
  ranges = IRanges(start = yx$start, end = yx$end), 
  strand = Rle(yx$strand), 
  exon = yx$exon_id, 
  gene = yx$geneName)
strand(exons) = "*"

countmat = NULL

bamFls = rep(NA,15)
sampleNums = c(1,11,23,2,32,33,3,40,42,43,47,53,55,56,58)
for(i in 1:15){
  bamFls[i] = paste("orbFrontalF",sampleNums[i],"/tophat/accepted_hits.bam",sep="")
  aln = readGAlignments(bamFls[i])
  strand(aln) = "*"
  counts = summarizeOverlaps(features = exons, reads = aln, mode = "Union")
  countmat = cbind(countmat, assays(counts)$count)
}

rownames(countmat) = ensYexons$exon_id
colnames(countmat) = paste0("orbFrontalF",sampleNums)
save(countmat, file="ensembl-exoncounts.rda")





###############################################
###### chunk 3b.2 -- run DESeq and EdgeR ######
###############################################
library(edgeR)
library(DESeq)
library(GenomicRanges)

# load per-exon counts:
load("ensembl-exoncounts.rda") #countmat
length(unique(rownames(countmat))) == nrow(countmat)
#TRUE, passes test.

# label samples:
sex = c(1,1,0,0,1,0,1,0,1,1,1,1,0,0,1)

# make an updated version of edgeR's function, so it will return the group abundances:
exactTest.more = function (object, pair = 1:2, dispersion = "auto", rejection.region = "doubletail", big.count = 900, prior.count.total = 0.5) 
{
  if (!is(object, "DGEList")) 
    stop("Currently only supports DGEList objects as the object argument.")
  if (length(pair) != 2) 
    stop("Pair must be of length 2.")
  rejection.region <- match.arg(rejection.region, c("doubletail", "deviance", "smallp"))
  group <- as.factor(object$samples$group)
  levs.group <- levels(group)
  if (is.numeric(pair)) 
    pair <- levs.group[pair]
  else pair <- as.character(pair)
  if (!all(pair %in% levs.group)) 
    stop("At least one element of given pair is not a group.\n Groups are: ", 
         paste(levs.group, collapse = " "), "\n")
  if (is.null(dispersion)) 
    dispersion <- "auto"
  if (is.character(dispersion)) {
    dispersion <- match.arg(dispersion, c("auto", "common", "trended", "tagwise"))
    dispersion <- switch(dispersion, common = object$common.dispersion, trended = object$trended.dispersion, tagwise = object$tagwise.dispersion, auto = getDispersion(object))
    if (is.null(dispersion)) 
      stop("specified dispersion not found in object")
  }
  ldisp <- length(dispersion)
  ntags <- nrow(object$counts)
  if (ldisp != 1 && ldisp != ntags) 
    stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
  if (ldisp == 1) 
    dispersion <- rep(dispersion, ntags)
  group <- as.character(group)
  j <- group %in% pair
  y <- object$counts[, j, drop = FALSE]
  lib.size <- object$samples$lib.size[j]
  norm.factors <- object$samples$norm.factors[j]
  group <- group[j]
  if (is.null(rownames(y))) 
    rownames(y) <- paste("tag", 1:ntags, sep = ".")
  lib.size <- lib.size * norm.factors
  offset <- log(lib.size)
  lib.size.average <- exp(mean(offset))
  abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
  logCPM <- (abundance + log(1e+06))/log(2)
  prior.count <- lib.size
  prior.count <- prior.count.total * prior.count/sum(prior.count)
  j1 <- group == pair[1]
  n1 <- sum(j1)
  if (n1 == 0) 
    stop("No libraries for", pair[1])
  y1 <- y[, j1, drop = FALSE]
  abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], ntags, n1, byrow = TRUE), offset = offset[j1])
  j2 <- group == pair[2]
  n2 <- sum(j2)
  if (n1 == 0) 
    stop("No libraries for", pair[2])
  y2 <- y[, j2, drop = FALSE]
  abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], ntags, n2, byrow = TRUE), offset = offset[j2])
  logFC <- (abundance2 - abundance1)/log(2)
  e <- exp(abundance)
  input.mean <- matrix(e, ntags, n1)
  output.mean <- input.mean * lib.size.average
  input.mean <- t(t(input.mean) * lib.size[j1])
  y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean, dispersion = dispersion)
  input.mean <- matrix(e, ntags, n2)
  output.mean <- input.mean * lib.size.average
  input.mean <- t(t(input.mean) * lib.size[j2])
  y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean, dispersion = dispersion)
  exact.pvals <- switch(rejection.region, doubletail = exactTestDoubleTail(y1, y2, dispersion = dispersion, big.count = big.count), deviance = exactTestByDeviance(y1, y2, dispersion = dispersion, big.count = big.count), smallp = exactTestBySmallP(y1, y2, dispersion = dispersion, big.count = big.count))
  de.out <- data.frame(logFC = logFC, logCPM = logCPM, PValue = exact.pvals, ab1 = exp(abundance1), ab2=exp(abundance2))
  rownames(de.out) <- rownames(object$counts)
  new("DGEExact", list(table = de.out, comparison = pair, genes = object$genes))
}


### edgeR - men v women
y = DGEList(counts = countmat, group = sex)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)
et = exactTest.more(y)
edger.table = topTags(et, n = dim(et$table)[1], adjust.method = "BH")
head(edger.table$table,n=20) 
edger.results = edger.table$table[is.finite(edger.table$table$logCPM),] ## remove exons with 0 counts in both groups
pdf(file="phist-edgeR-sex.pdf")
    hist(edger.results$PValue, col="gray70",xlab="p value",main="EdgeR p values - male vs. female, Y", breaks=30)
dev.off()
save(edger.results,file="edger.results_ensembl.rda") ## in working directory

# edgeR - men v. men
men.table = countmat[,which(sex==1)]
set.seed(651)
y2 = DGEList(counts = men.table, group = sample(c(rep(1,5),rep(0,4))) )
y2 = estimateCommonDisp(y2)
y2 = estimateTagwiseDisp(y2)
et2 = exactTest.more(y2)
edger.table.men = topTags(et2, n = dim(et2$table)[1], adjust.method = "BH")
head(edger.table.men$table,n=10)
edger.results.men = edger.table.men$table[is.finite(edger.table$table$logCPM),]
pdf(file="phist-edgeR-men.pdf")
    hist(edger.results.men$PValue, col="gray70",xlab="p value",main="Edger p values: men only, Y")
dev.off()
save(edger.results.men,file="edger.results.men_ensembl.rda")

# DESeq - men v women
cds = newCountDataSet(countmat, sex)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
des.res = nbinomTest(cds, "1","0")
head(des.res)
deseq.table = des.res[order(des.res$padj),]
pdf(file="phist-DESeq-sex.pdf")
    hist(deseq.table$pval, col="gray70",xlab="p value",main="DESeq p values: male vs. female, Y", breaks=30)
dev.off()
save(deseq.table,file="deseq.table_ensembl.rda")

# DESeq - men v men
set.seed(651)
cds2 = newCountDataSet(men.table, sample(c(rep(1,5),rep(0,4))))
cds2 = estimateSizeFactors(cds2)
cds2 = estimateDispersions(cds2)
des.res.men = nbinomTest(cds2,"1","0")
head(des.res.men)
deseq.table.men = des.res.men[order(des.res.men$padj),]
pdf(file="phist-DESeq-men.pdf")
    hist(deseq.table.men$pval,col="gray70",xlab="p value",main="DESeq p values: men only, Y")
dev.off()
save(deseq.table.men,file="deseq.table.men_ensembl.rda")


###############################################
###### chunk 3b.3 -- compare to DER Finder ####
###############################################

###### make edger/deseq table: (supplementary table 3)
annottable = matrix(NA,ncol=9, nrow=17)
annottable = as.data.frame(annottable)
names(annottable) = c("q","numregions","num.exons.us", "num.exons.edgeR","num.exons.DESeq", "edger.us.agree", "deseq.us.agree", "deseq.edger.agree", "all.agree")
annottable$q = seq(0,0.8,by=0.05)

int2 = function(z) Reduce('intersect',z) # need to take intersection of more than 2 sets -- this is how you do it.

for(r in 1:17){
  annottable$numregions[r] = sum(ders$qvals<annottable$q[r])
  annottable$num.exons.us[r] = length(unique(myinfo.sub$ex.names[myinfo.sub$ex.pct>=0.8 & myinfo.sub$qvals.ex<annottable$q[r]]))
  ourexons = unique(myinfo.sub$ex.names[myinfo.sub$ex.pct>=0.8 & myinfo.sub$qvals.ex<annottable$q[r]])
  annottable$num.exons.edgeR[r] = sum(edger.results$FDR<annottable$q[r])
  edgerexons = rownames(edger.results)[which(edger.results$FDR<annottable$q[r])]
  # need to fix exon names (included gene names in edgeR/deseq tables)
  ulength = length(edgerexons)
  edgerexons = unlist(lapply(edgerexons, function(x) strsplit(x,"-")[[1]][1]))
  ulength2 = length(unique(edgerexons))
  if(ulength!=ulength2) print(paste("OH NO EDGER",r))
  annottable$num.exons.DESeq[r] = length(which((deseq.table$padj<annottable$q[r])))
  deseqexons = deseq.table$id[which((deseq.table$padj<annottable$q[r]))]
  ulength = length(deseqexons)
  deseqexons = unlist(lapply(deseqexons, function(x) strsplit(x,"-")[[1]][1]))
  ulength2 = length(unique(deseqexons))
  if(ulength!=ulength2) print(paste("OH NO DESEQ",r))
  annottable$edger.us.agree[r] = length(intersect(edgerexons,ourexons))
  annottable$deseq.us.agree[r] = length(intersect(deseqexons,ourexons))
  annottable$deseq.edger.agree[r] = length(int2(list(deseqexons, edgerexons)))
  annottable$all.agree[r] = length(int2(list(edgerexons,ourexons,deseqexons)))
}

annottable

length(setdiff(myinfo.sub$ex.names, rownames(edger.results))) ##345

## investigate a few of the exons that EdgeR/DESeq find, but we do not.
r=2 #set q=0.05
ourexons = unique(myinfo.sub$ex.names[myinfo.sub$ex.pct>=0.8 & myinfo.sub$qvals.ex<annottable$q[r]])
edgerexons = rownames(edger.results)[which(edger.results$FDR<annottable$q[r])]
deseqexons = deseq.table$id[which((deseq.table$padj<annottable$q[r]))]
edgeronly = setdiff(edgerexons,ourexons)
length(edgeronly) #47 exons
deseqonly = setdiff(deseqexons,ourexons) 
length(deseqonly) #39 exons
length(union(edgeronly, deseqonly)) ##54 discovered only by edger or deseq


## setup for plotting:
gro = list()
gro$states = regions.merged.y
states.norle.temp = inverse.rle(list(lengths=regions.merged.y$length, values=regions.merged.y$state)) 
load("posY-rev.rda") # should already be loaded, but just in case
load("ttY-rev.rda") # should already be loaded, but just in case
states.norle.temp2 = states.norle.temp[pos-pos[1]+1]
gro$states.norle=data.frame(pos=pos, states=states.norle.temp2)
rm(states.norle.temp, states.norle.temp2);gc();gc();gc()
group = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1) 
group.l = ifelse(group==1,"male","female")
ensYexons.forplot = ensYexons
names(ensYexons.forplot)[1] = "geneName"
names(ensYexons.forplot)[3] = "name"

# plot those exons found only by edgeR/DEseq:
pdf("only_edger_ens.pdf")
for(ex in edgeronly){
  plotExon(gro, exonname=ex, tstats=tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l, scalefac=32,ylim=c(4.5,7.5)) 
}
dev.off()

pdf("only_deseq_ens.pdf")
for(ex in deseqonly){
  plotExon(gro, exonname=ex, tstats=tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l, scalefac=32,ylim=c(4.5,7.5)) 
}
dev.off()

length(intersect(deseqonly, edgeronly)) #32

### look at these for page 14 of DER Finder manuscript

# some miscellaneous statistics (included in results/text of manuscript)
sum(ders$qvals<0.05) #534
sum(ders$qvals<0.05 & ders$state==4) #6
length(which(ders$flag=="novel" & ders$qvals<0.05)) #280
summary(ders$length[which(ders$flag=="novel" & ders$qvals<0.05)]) #range: 1-3814
de.exons = myinfo.sub$ex.names[which(myinfo.sub$qvals.ex<0.05 & myinfo.sub$ex.pct>=0.8)]
length(unique(de.exons)) #411, and "novel" is not included.
genes.represented = ensYexons$geneName[which(ensYexons$exon_id %in% unique(de.exons))]
length(unique(genes.represented)) #33
min(ders.men$qvals) #0.86
min(tx.sex.cuff.y.ok$qvalue) #0.45
min(tx.null.cuff.y.ok$qvalue) #0.63
sum(tx.sex.cuff.y.ok$value_female==0 & tx.sex.cuff.y.ok$value_male!=0) #736




#########################################
#########################################
##### chunk 4: figures ##################
#########################################
#########################################


#### FIGURE 1
load("chr22exons.rda")
xx = exondata22[exondata22$exon_chrom_start>20936000 & exondata22$exon_chrom_end<20946000 & exondata22$exon_chrom_start < 20941000,]
dim(xx) #108 rows
length(unique(xx$ensembl_exon_id)) #40 unique exons by id
length(unique(xx$ensembl_transcript_id)) #15 transcripts

setEPS(width=12, height=6)
postscript("exons.eps")
#colored_exons = NULL
par(mfrow=c(1,2))

### transcript structures
xax = seq(min(xx$exon_chrom_start), max(xx$exon_chrom_end), by=1)
plot(xax, rep(0,length(xax)), ylim=c(0,17), type="n", xlab="Genomic Position", yaxt = "n", ylab="")
title("(a) Annotated Transcripts: Ensembl 61, Chromosome 22")
for(tx in unique(as.character(xx$ensembl_transcript_id))){
    txind = which(unique(xx$ensembl_transcript_id)==tx)
    gtsub = xx[xx$ensembl_transcript_id==tx,]
    gtsub = gtsub[order(gtsub$exon_chrom_start),]
    for(exind in 1:nrow(gtsub)){
      #mycolor = ifelse(as.character(gtsub$ensembl_exon_id[exind]) %in% colored_exons, "white","gray60")
      mycolor = "gray60"
      polygon(x=c(gtsub$exon_chrom_start[exind], 
        gtsub$exon_chrom_start[exind], 
        gtsub$exon_chrom_end[exind], 
        gtsub$exon_chrom_end[exind]), 
        y=c(txind-0.4,txind+0.4,txind+0.4,txind-0.4), 
        col=mycolor)
      if(exind != nrow(gtsub)) lines(c(gtsub$exon_chrom_end[exind],gtsub$exon_chrom_start[exind+1]),c(txind, txind), lty=2, col="gray60")
      #colored_exons = append(as.character(gtsub$ensembl_exon_id[exind]), colored_exons)  
    }
}
#legend("bottomright", pch=15, col="gray60", "unique exon")

### overlapping exons
olap_exons = subset(xx, exon_chrom_start<20941000 & exon_chrom_start>20940500 &
  exon_chrom_end<20942000 & exon_chrom_end>20941500)
olap_exons = olap_exons[which(!duplicated(olap_exons[,c(5,6)])),]
xax = seq(olap_exons$exon_chrom_start[1]-100, max(olap_exons$exon_chrom_end)+100, by=1)
plot(xax, rep(0,length(xax)), ylim=c(0.5,4.5), 
  type="n", xlab="Genomic Position", yaxt = "n", ylab="", 
  xlim=c(20941800, 20941950))
for(i in seq_along(olap_exons[,1])){
  polygon(x=c(olap_exons$exon_chrom_start[i], 
    olap_exons$exon_chrom_start[i],
    olap_exons$exon_chrom_end[i],
    olap_exons$exon_chrom_end[i]), 
    y=c(i-0.4, i+0.4, i+0.4, i-0.4),
    col="gray60")
}
title("(b) Close-up of Exon Annotation Differences")
dev.off()





# percentile vs. %MF plot [FIGURE 4]
regions.menA = regions.men
regions.A = regions
ders.men = subset(regions.men, state==3|state==4)
ders = subset(regions, state==3|state==4)
regions = ders
regions.men = ders.men
regions$exp = rep("sex",dim(regions)[1])
regions.men$exp = rep("men",dim(regions.men)[1])
us.t = c(regions$mean.t,regions.men$mean.t)
us.exp = c(regions$exp, regions.men$exp)
us.state = c(regions$state, regions.men$state)
us.length = c(regions$length, regions.men$length)
usp = regions$pvals
usp[regions$mean.t<0] = -usp[regions$mean.t<0]
unp = regions.men$pvals
unp[regions.men$mean.t<0] = -unp[regions.men$mean.t<0] # pos t means overexpressed in men
us.p = c(usp,unp)
us.res = data.frame(t = us.t, exp = us.exp, p=us.p)
us.res = us.res[us.length>=102,]
us.res = us.res[order(1/us.res$p,decreasing=T),]
# handle ties (randomly sample labels):
p.rle = rle(us.res$p)
p.rle = data.frame(lengths=p.rle$lengths,values=p.rle$values)
for(k in which(p.rle$lengths>1)){
  inds = which(us.res$p==p.rle$values[k])
  us.res$exp[inds] = sample(us.res$exp[inds])
}

us.pct = NULL
for(i in 1:nrow(us.res)){
  us.pct[i] = sum(us.res$exp[1:i]=="sex")/i
}
t.percentile = NULL
for(i in 1:length(us.res$t)){ t.percentile[i] = 1-i/length(us.res$t)}


cuff.t = tx.sex.cuff.y.ok$test_stat
cuff.t = ifelse(cuff.t>20,20,cuff.t)
cuff.t = ifelse(cuff.t<(-20),-20,cuff.t)
cuff.fchange = tx.sex.cuff.y.ok$log2fchange
cuff.fchange = ifelse(cuff.fchange>20,20,cuff.fchange)
cuff.fchange = ifelse(cuff.fchange<(-20),-20,cuff.fchange)

cuff.t.m = tx.null.cuff.y.ok$test_stat
cuff.t.m = ifelse(cuff.t.m>20,20,cuff.t.m)
cuff.t.m = ifelse(cuff.t.m<(-20),-20,cuff.t.m)
cuff.fchange.m = tx.null.cuff.y.ok$log2fchange
cuff.fchange.m = ifelse(cuff.fchange.m>20,20,cuff.fchange.m)
cuff.fchange.m = ifelse(cuff.fchange.m<(-20),-20,cuff.fchange.m)

cuff.t.all = c(cuff.t,cuff.t.m)
cuff.fchange[which(tx.sex.cuff.y.ok$log2fchange<(-20))] = -20
cuff.fc.all = c(cuff.fchange, cuff.fchange.m)
cuff.p.tmp = tx.sex.cuff.y.ok$pvalue
# positive p-values
cuff.p.tmp[which(tx.sex.cuff.y.ok$test_stat<(-20) | (tx.sex.cuff.y.ok$test_stat>0 & tx.sex.cuff.y.ok$test_stat<=20))] = 1-cuff.p.tmp[which(tx.sex.cuff.y.ok$test_stat<(-20) | (tx.sex.cuff.y.ok$test_stat>0 & tx.sex.cuff.y.ok$test_stat<=20))]
# negative p-values
cuff.p.tmp[which(tx.sex.cuff.y.ok$test_stat>20 | (tx.sex.cuff.y.ok$test_stat<0 & tx.sex.cuff.y.ok$test_stat>(-20)))] = cuff.p.tmp[which(tx.sex.cuff.y.ok$test_stat>20 | (tx.sex.cuff.y.ok$test_stat<0 & tx.sex.cuff.y.ok$test_stat>(-20)))] - 1
head(data.frame(tx.sex.cuff.y.ok, cuff.p.tmp),100)

cuff.p.tmp.null = tx.null.cuff.y.ok$pvalue
# positives:
cuff.p.tmp.null[which(tx.null.cuff.y.ok$test_stat<(-20) | (tx.null.cuff.y.ok$test_stat>0 & tx.null.cuff.y.ok$test_stat<=20))] = 1-cuff.p.tmp.null[which(tx.null.cuff.y.ok$test_stat<(-20) | (tx.null.cuff.y.ok$test_stat>0 & tx.null.cuff.y.ok$test_stat<=20))]
# negatives:
cuff.p.tmp.null[which(tx.null.cuff.y.ok$test_stat>20 | (tx.null.cuff.y.ok$test_stat<0 & tx.null.cuff.y.ok$test_stat>(-20)))] = cuff.p.tmp.null[which(tx.null.cuff.y.ok$test_stat>20 | (tx.null.cuff.y.ok$test_stat<0 & tx.null.cuff.y.ok$test_stat>(-20)))]-1
head(data.frame(tx.null.cuff.y.ok, cuff.p.tmp.null), 100)

cuff.p.all = c(cuff.p.tmp, cuff.p.tmp.null)
cuff.exp = c(rep("sex",length(cuff.t)), rep("men",length(cuff.t.m)))
cuff.res = data.frame(t=cuff.t.all, exp = cuff.exp, p = cuff.p.all)
cuff.res = cuff.res[order(cuff.res$p,decreasing=T),]

cuff.p.rle = rle(cuff.res$p)
cuff.p.rle = data.frame(lengths=cuff.p.rle$lengths,values=cuff.p.rle$values)
for(k in which(cuff.p.rle$lengths>1)){
  inds = which(cuff.res$t==cuff.p.rle$values[k])
  cuff.res$exp[inds] = sample(cuff.res$exp[inds])
}
cuff.pct = NULL
for(i in 1:nrow(cuff.res)){
  cuff.pct[i] = sum(cuff.res$exp[1:i]=="sex")/i
}
p.percentile = NULL
for(i in 1:length(cuff.res$p)){ p.percentile[i] = 1-i/length(cuff.res$p)}


edger.results$exp = rep("sex",nrow(edger.results))
edger.results.men$exp = rep("men",nrow(edger.results.men))
edger.fc = c(edger.results$logFC,edger.results.men$logFC)
esp = edger.results$PValue
esp[edger.results$logFC<0] = -esp[edger.results$logFC<0]
enp = edger.results.men$PValue
enp[edger.results.men$logFC<0] = -enp[edger.results.men$logFC<0] # has positive fold changes indicating overexpression in men

edger.p = c(esp, enp)
edger.exp = c(edger.results$exp, edger.results.men$exp)
edger.logCPM = c(edger.results$logCPM, edger.results.men$logCPM)
edger.res = data.frame(fc = edger.fc, exp = edger.exp, p=edger.p)
edger.res = edger.res[is.finite(edger.logCPM),]
edger.res = edger.res[order(1/edger.res$p,decreasing=T),]
# handle ties (randomly sample labels):
p.rle = rle(edger.res$p)
p.rle = data.frame(lengths=p.rle$lengths,values=p.rle$values)
for(k in which(p.rle$lengths>1)){
  inds = which(edger.res$p==p.rle$values[k])
  edger.res$exp[inds] = sample(edger.res$exp[inds])
}
edger.pct = NULL
for(i in 1:nrow(edger.res)){
  edger.pct[i] = sum(edger.res$exp[1:i]=="sex")/i
}

ep.percentile = NULL
for(i in 1:length(edger.res$p)){ ep.percentile[i] = 1-i/length(edger.res$p)}

deseq.table$exp = rep("sex",nrow(deseq.table))
deseq.table.men$exp = rep("men",nrow(deseq.table.men))
deseq.fc = c(-deseq.table$log2FoldChange,-deseq.table.men$log2FoldChange)
dsp = deseq.table$pval
dsp[which(deseq.table$log2FoldChange>0)] = -dsp[which(deseq.table$log2FoldChange>0)] #DESeq has negative fold changes indicating overexpression in men.
dnp = deseq.table.men$pval
dnp[which(deseq.table.men$log2FoldChange>0)] = -dnp[which(deseq.table.men$log2FoldChange>0)]
deseq.p = c(dsp, dnp)

deseq.exp = c(deseq.table$exp, deseq.table.men$exp)
deseq.res = data.frame(fc = deseq.fc, exp = deseq.exp, p=deseq.p)
deseq.res = deseq.res[order(1/deseq.p,decreasing=T),]
deseq.res = deseq.res[!is.nan(deseq.res$fc),]

# handle ties (randomly sample labels):
p.rle.d = rle(deseq.res$p)
p.rle.d = data.frame(lengths=p.rle.d$lengths,values=p.rle.d$values)
for(k in which(p.rle.d$lengths>1)){
  inds = which(deseq.res$p==p.rle.d$values[k])
  deseq.res$exp[inds] = sample(deseq.res$exp[inds])
}
deseq.pct = NULL
for(i in 1:nrow(deseq.res)){
  deseq.pct[i] = sum(deseq.res$exp[1:i]=="sex")/i
}
d.percentile = NULL
for(i in 1:length(deseq.res$p)){ d.percentile[i] = 1-i/length(deseq.res$p)}


## finally, the plot:
#mycols = c("dodgerblue2","orange","green3","darkorchid3")
plot(t.percentile,us.pct,type="l",ylim=c(0.55,1),xlab="p-value percentile",ylab="Percent from Male vs. Female Comparisons",lwd=2,xlim=c(0.5,1), lty=1)
lines(p.percentile,cuff.pct,lwd=2, lty=2)
lines(ep.percentile, edger.pct,lwd=2, lty=3)
lines(d.percentile, deseq.pct,lwd=2, lty=6)
legend("topleft",lty=c(1,2,3,6),c("DER Finder","Cufflinks","EdgeR","DESeq"),lwd=c(4,4,4,4))

plot(1-t.percentile,us.pct,type="l",xlab="p-value percentile",ylab="Percent from Male vs. Female Comparisons",lwd=2,xlim=c(0,0.4),ylim=c(0.55,1),xaxt="n",lty=1)
axis(1,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2))
lines(1-p.percentile,cuff.pct,lty=4,lwd=2)
lines(1-ep.percentile, edger.pct,lty=3,lwd=2)
lines(1-d.percentile, deseq.pct,lty=2,lwd=2)
legend("bottomleft",lty=c(1,4,3,2),c("DER Finder","Cufflinks","EdgeR","DESeq"),lwd=c(2,2,2,2))

### make this as an EPS -- THIS WAS THE ONE USED IN THE PAPER
#### COLOR
setEPS()
postscript("p-rank-plot-zoom.eps")
plot(1-t.percentile,us.pct,type="l",xlab="p-value percentile",ylab="Percent from Male vs. Female Comparisons",lwd=4,col=mycols[1],xlim=c(0,0.4),ylim=c(0.55,1),xaxt="n")
axis(1,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2))
lines(1-p.percentile,cuff.pct,col=mycols[2],lwd=4)
lines(1-ep.percentile, edger.pct,col=mycols[3],lwd=4)
lines(1-d.percentile, deseq.pct,col=mycols[4],lwd=4)
legend("bottomleft",col=mycols,c("DER Finder","Cuffdiff","EdgeR","DESeq"),lwd=c(4,4,4,4))
dev.off()

#### BLACK AND WHITE
setEPS()
postscript("p-rank-plot-zoom-bw.eps")
plot(1-t.percentile,us.pct,type="l",xlab="p-value percentile",ylab="Percent from Male vs. Female Comparisons",lwd=2,xlim=c(0,0.4),ylim=c(0.55,1),xaxt="n",lty=1)
axis(1,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2))
lines(1-p.percentile,cuff.pct,lty=4,lwd=2)
lines(1-ep.percentile, edger.pct,lty=3,lwd=2)
lines(1-d.percentile, deseq.pct,lty=2,lwd=2)
legend("bottomleft",lty=c(1,4,3,2),c("DER Finder","Cufflinks","EdgeR","DESeq"),lwd=c(2,2,2,2))
dev.off()



#(b)  MA plots [FIGURE 3]

# for DER Finder - get mean coverage for each sample for each region.
coverage.means = matrix(nrow = dim(regions.merged.y)[1], ncol = 15)
dbfile = "Y-tophat-revised.db"
tablename = "chrY"
library(Genominator)
coverage.file = ExpData(dbFilename = dbfile,tablename="chrY")

for(i in 1:dim(regions.merged.y)[1]){
  print(i)
  if(regions.merged.y$state[i]==1) {coverage.means[i,] <- 0; next}
  firstind = which(pos==regions.merged.y$start[i])
  lastind = which(pos==regions.merged.y$end[i])
  coverage.means[i,] = colMeans(coverage.file[firstind:lastind,-1])
}

save(coverage.means, file="coverage.means.sex-rev.rda")
# load("coverage.means.sex-rev.rda") # if above has already been done


# same for men:
coverage.means.men = matrix(nrow = dim(regions.merged.men)[1], ncol=15)
# this will have data for all the samples, but we're only using one at a time and we won't use a female one.
for(i in 1:dim(regions.merged.men)[1]){
  print(i)
  if(regions.merged.men$state[i]==1) {coverage.means.men[i,] <- 0; next}
  firstind = which(pos==regions.merged.men$start[i])
  lastind = which(pos==regions.merged.men$end[i])
  coverage.means.men[i,] = colMeans(coverage.file[firstind:lastind,-1])
}
save(coverage.means.men, file="coverage.means.men-rev.rda")
# load("coverage.means.men-rev.rda") # if the above has already been done


png("MAplots.png", width=700, height=900) ## need to convert the png --> eps with conversion tool
par(mfrow=c(4,3),oma=c(0,0,0,0),mar=c(1.5,1.8,1.5,1),mgp=c(0.7,0.1,0))

## us:
male.avg = rowMeans(log2(coverage.means[,-c(3,4,7,8,13,14)]+32))
male.avg.nolog = rowMeans(coverage.means[,-c(3,4,7,8,13,14)])
female.avg.nolog = rowMeans(coverage.means[,c(3,4,7,8,13,14)])
maleA.avg.nolog = rowMeans(coverage.means[,c(1,2,9,12,15)])
maleB.avg.nolog = rowMeans(coverage.means[,c(5,6,10,11)])
all.avg = rowMeans(log2(coverage.means+32))
male.diff = log2(maleA.avg.nolog+32)-log2(maleB.avg.nolog+32)
sex.diff = log2(male.avg.nolog+32)-log2(female.avg.nolog+32)
plot(all.avg, sex.diff,xlab="A",ylab="M",col="#FF000050",main="DER Finder - sex",xlim=c(5,10),ylim=c(-1,5),pch=19, cex.axis=0.7,tck=-0.005)
plot(male.avg,male.diff,col="#0000FF50",xlab="A",ylab="M",main="DER Finder - males",xlim=c(5,10),ylim=c(-1,5),pch=19, cex.axis=0.7,tck=-0.005)
plot(all.avg, sex.diff,xlab="A",ylab="M",col="#FF000050",main="DER Finder - overlaid",xlim=c(5,10),ylim=c(-1,5),pch=19, cex.axis=0.7,tck=-0.005)
points(male.avg,male.diff,col="#0000FF50",pch=19)

cufflinksA = (log2(tx.sex.cuff.y.ok$value_male+0.5)+log2(tx.sex.cuff.y.ok$value_female+0.5))/2
cufflinksM = log2(tx.sex.cuff.y.ok$value_male+0.5)-log2(tx.sex.cuff.y.ok$value_female+0.5)
cufflinksA.null = (log2(tx.null.cuff.y.ok$value_m1+0.5)+log2(tx.null.cuff.y.ok$value_m2+0.5))/2
cufflinksM.null = log2(tx.null.cuff.y.ok$value_m1+0.5)-log2(tx.null.cuff.y.ok$value_m2+0.5)
plot(cufflinksA, cufflinksM, xlab="A", ylab="M", main="Cuffdiff - sex", col="#FF000050", ylim=c(-13,13), xlim=c(-1,8),pch=19, cex.axis=0.7,tck=-0.005)
plot(cufflinksA.null, cufflinksM.null, xlab="A", ylab="M", main="Cuffdiff - males", col="#0000FF50", ylim=c(-13,13), xlim=c(-1,8),pch=19, cex.axis=0.7,tck=-0.005)
plot(cufflinksA, cufflinksM, xlab="A", ylab="M", main="Cuffdiff - overlaid", col="#FF000050", ylim=c(-13,13), xlim=c(-1,8),pch=19, cex.axis=0.7,tck=-0.005)
points(cufflinksA.null, cufflinksM.null, xlab="A", ylab="M", col="#0000FF50",pch=19)

## edgeR
plot(edger.results$logCPM, edger.results$logFC,col="#FF000050",xlab="A",ylab="M",ylim=c(-7,16),xlim=c(0,18),main="EdgeR - sex",pch=19, cex.axis=0.7,tck=-0.005)
plot(edger.results.men$logCPM, edger.results.men$logFC, col="#0000FF50", xlab="A",ylab="M",ylim=c(-7,16), xlim=c(0,18),main="EdgeR - males",pch=19, cex.axis=0.7,tck=-0.005)
plot(edger.results$logCPM, edger.results$logFC,col="#FF000050",xlab="A",ylab="M",ylim=c(-7,16),xlim=c(0,18),main="EdgeR - overlaid",pch=19,  cex.axis=0.7,tck=-0.005)
points(edger.results.men$logCPM, edger.results.men$logFC,col="#0000FF50",pch=19)

## DESeq:
plot((log2(deseq.table$baseMeanA+0.5)+log2(deseq.table$baseMeanB+0.5))/2, log2(deseq.table$baseMeanA+0.5)-log2(deseq.table$baseMeanB+0.5),xlab="A",ylab="M",main="DESeq - sex",col="#FF000050",ylim=c(-5,13), xlim=c(-1,13),pch=19,  cex.axis=0.7,tck=-0.005)
plot((log2(deseq.table.men$baseMeanA+0.5)+log2(deseq.table.men$baseMeanB+0.5))/2, log2(deseq.table.men$baseMeanA+0.5)-log2(deseq.table.men$baseMeanB+0.5),xlab="A",ylab="M",main="DESeq - males",col="#0000FF50",xlim=c(-1,13), ylim=c(-4,13),pch=19,  cex.axis=0.7,tck=-0.005)
plot((log2(deseq.table$baseMeanA+0.5)+log2(deseq.table$baseMeanB+0.5))/2, log2(deseq.table$baseMeanA+0.5)-log2(deseq.table$baseMeanB+0.5),xlab="A",ylab="M",main="DESeq - overlaid",col="#FF000050",ylim=c(-5,13), xlim=c(-1,13),pch=19,  cex.axis=0.7,tck=-0.005)
points((log2(deseq.table.men$baseMeanA+0.5)+log2(deseq.table.men$baseMeanB+0.5))/2, log2(deseq.table.men$baseMeanA+0.5)-log2(deseq.table.men$baseMeanB+0.5),col="#0000FF50",pch=19)
dev.off()




### p-value histograms [supplementary figure 6]
setEPS()
postscript("phists-4.eps")
par(mfrow=c(2,2),mar=c(1.9,1.8,1.5,1),mgp=c(0.7,0.1,0))

hist(pvals,breaks=30,col="gray70",main="",xlab="p values",cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
text(0.5,1030,"DER Finder",cex=1.2,font=2) 
#hist(pvals.men,breaks=30,col="gray70",main="",xlab="p values",ylim=c(0,1000),cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
#text(0.5,900,"DER Finder - males",cex=1.2,font=2) 
hist(tx.sex.cuff.y.ok$pvalue, col="gray70",breaks=30, xlab="p values",main="",ylim=c(0,350),cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
text(0.5,350,"Cuffdiff",cex=1.2,font=2) 
#hist(tx.null.cuff.y.ok$pvalue, col="gray70",breaks=30, xlab="p values",main="",ylim=c(0,300),cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
#text(0.5,300,"Cufflinks - males",cex=1.2,font=2) 
hist(edger.results$PValue[is.finite(edger.results$logCPM)], col="gray70",xlab="p value",main="",breaks=30,cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
text(0.5,155,"EdgeR",cex=1.2,font=2) 
#hist(edger.results.men$PValue[is.finite(edger.results.men$logCPM)], col="gray70",xlab="p value",main="",breaks=30,ylim=c(0,25),cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
#text(0.5,23,"EdgeR - males",cex=1.2,font=2) 
hist(deseq.table$pval,col="gray70",xlab="p value",main="",breaks=30,cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
text(0.5,130,"DESeq",cex=1.2,font=2) 
#hist(deseq.table.men$pval,col="gray70",xlab="p value",main="",breaks=30,cex.axis=0.7,tck=-0.005, cex.lab = 0.9)
#text(0.5,40,"DESeq - males",cex=1.2,font=2) 

dev.off()


## find some exons:
#(make sure to first run the "setup for plotting" section earlier)

myinfo.cmp = subset(myinfo.sub, ex.names %in% rownames(edger.results))
plotExon(gro, exonname = "ENSE00001435537", tstats = tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l, scalefac=32, ylim=c(4.2,8), bppad=500)

myinfo.cmp[108,] #q=0.001, 61.3% of an exon
edger.results[which(rownames(edger.results)=="ENSE00001435537"),] #q=1
des.res[des.res$id=="ENSE00001435537",] #q=1

which(regions.merged.y$length>400 & regions.merged.y$state==3)[100]
plotRegion(gro, ind=13047, tstats=tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="chrY", group=group.l,scalefac=32, ylim=c(4.5, 7.5))

setEPS()
postscript("wrong-exon-II.eps")
plotExon(gro, exonname = "ENSE00001435537", tstats = tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l, scalefac=32, ylim=c(4.2,8), bppad=500)
dev.off()

regions.merged.y[45,]
ders[which(ders$start== 2714305),] #q = 0.006

derinds = as.numeric(rownames(ders))
derinds.interesting = derinds[which(ders$length>250 & ders$flag=="novel DE region")] #these are not indices but are row names for regions.merged.y
plotinds = which(rownames(regions.merged.y) %in% derinds.interesting)
pdf("novel_plots.pdf")
for(i in plotinds){
  plotRegion(gro, ind=i, tstats=tt, pos=pos, annotation=Yexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="chrY", group=group.l,scalefac=32,ylim=c(5,7), legendloc = "topleft")
  
}
dev.off()

# after examining the "novel_plots.pdf" file:
which(regions.merged.y$start==20662506) #11154
plotRegion(gro, ind=11154, tstats=tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="chrY", group=group.l,scalefac=32,ylim=c(5,7), legendloc = "topleft")
# great, plus overlaps ESTs but no genes.
# get the q-value:
which(ders$start==20662506) #1670
ders[1670,] #q=0.001


## combine these two figures into one: [FIGURE 2 in paper]
source("multipanel.R") # just allows you to plot >1 pretty plot in same frame (i.e. par call is removed -- otherwise same as defined plotExon function)

setEPS(width=12, height=6)
postscript("examples-combined.eps")
par(mfcol=c(3,2), mar=c(1,3.2,2,1), mgp=c(1.5,0.5,0), cex.lab=1.5, omi=c(0.4,0,0,0))
### panel (a)
plotExon.nopar(gro, exonname = "ENSE00001435537", tstats = tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l, scalefac=32, ylim=c(4.9,8.3), bppad=500, plottitle = "(a) chrY: 22737611 - 22737773")

### panel (b)
plotRegion.nopar(gro, ind=11154, tstats=tt, pos=pos, annotation=ensYexons.forplot, counts="Y-tophat-revised.db", tabname="chrY", chromosome="Y", group=group.l,scalefac=32,ylim=c(5,7.3), legendloc = "topright", plottitle = "(b) chrY: 20662506 - 20662937")
mtext("genomic position", side=1, outer=TRUE, adj=0.22, padj=1)
mtext("genomic position", side=1, outer=TRUE, adj=0.8, padj=1)

dev.off()




### various misc statistics
sum(!is.nan(deseq.table$foldChange))
sum(is.finite(edger.results$logCPM))
sum(edger.results$FDR<0.05) #113
length(which(deseq.table$padj<0.05)) #115
setdiff(rownames(edger.results[edger.results$FDR<0.05,]), deseq.table$id[which(deseq.table$padj<0.05)])
setdiff(deseq.table$id[which(deseq.table$padj<0.05)], rownames(edger.results[edger.results$FDR<0.05,]))
length(intersect(deseq.table$id[which(deseq.table$padj<0.05)], rownames(edger.results[edger.results$FDR<0.05,]))) #97

sum(edger.results.men$FDR<0.05) #0
length(which((deseq.table.men$padj<0.05))) #0
min(deseq.table.men$padj, na.rm=T) #1
min(edger.results.men$FDR) #1

sum(ders$state==3|ders$state==4 & ders$qvals<0.05)
sum(ders$state==4 & ders$qvals<0.05) #30
sum(ders$state==3 & ders$qvals<0.05) #30
length(which((myinfo.updated$ex.class=="novel")))



## a few miscellaneous statistics: 
edger.sig = subset(edger.results, FDR<0.05)
deseq.sig = subset(deseq.table, padj<0.05)
ders.sig = subset(ders,qvals<0.05)
edger.sig.locs = NULL
for(i in 1:dim(edger.sig)[1]){
  ss = strsplit(rownames(edger.sig)[i],"-")[[1]]
  exname = ss[1]
  gname = ss[2]
  exind = which(Yexons$gene==gname & Yexons$exon_id==exname)
  edger.sig.locs = append(edger.sig.locs, Yexons$start[exind]:Yexons$end[exind])
}
deseq.sig.locs = NULL
for(i in 1:dim(deseq.sig)[1]){
  ss = strsplit(deseq.sig$id[i],"-")[[1]]
  exname = ss[1]
  gname = ss[2]
  exind = which(Yexons$gene==gname & Yexons$exon_id==exname)
  deseq.sig.locs = append(deseq.sig.locs, Yexons$start[exind]:Yexons$end[exind])
}
us.sig.locs = NULL
for(i in 1:dim(ders.sig)[1]){
  us.sig.locs = append(us.sig.locs, ders.sig$start[i]:ders.sig$end[i])
}

sum(us.sig.locs %in% deseq.sig.locs)/length(us.sig.locs) #0.336
sum(us.sig.locs %in% edger.sig.locs)/length(us.sig.locs) #0.337



#########################################
#########################################
##### chunk 5: supplementary figures ####
#########################################
#########################################

#(1) given our exons at each percent level, what percent of those did EdgeR discover?  (upward sloping)

erexnames = rownames(edger.results)
deexnames = deseq.table$id
x = seq(0,1,by=0.05) 
min(myinfo.updated$ex.pct, na.rm=TRUE) #0.01
derf.num = edger.num = deseq.num = NULL
for(i in 2:length(x)){
  myinds = which(myinfo.updated$ex.pct <= x[i] & myinfo.updated$qvals.ex<0.05)
  derf.num[i] = length(myinds)
  myexons = myinfo.updated$ex.names[myinds]
  ERinds = which(erexnames %in% myexons)
  DESinds = which(deexnames %in% myexons)
  edger.num[i] = length(which(edger.results$FDR[ERinds] < 0.05))
  deseq.num[i] = length(which(deseq.table$padj[DESinds] < 0.05))
}
plot(derf.num, edger.num)
points(derf.num, deseq.num, col="blue")
edger.num[1:2] = deseq.num[1:2] = c(0,0)
derf.num[1:2] = c(1,1) #dummy, for plot

setEPS() #[SUPPLEMENTARY FIGURE 1]
postscript("suppfig1.eps")
plot(x, edger.num/derf.num-0.01,type="l",ylim=c(0,1),lwd=2,col="blue",xlab="Maximum Percentage of Exon Covered", ylab="Percentage Discovered by Identify-then-annotate Methods")
points(x, deseq.num/derf.num,col="red",type="l",lwd=2)
legend("bottomright",col=c("blue","red"), c("EdgeR", "DESeq"), lty=c(1,1), lwd=c(2,2))
dev.off()


#(2) given edgeR's exons, what percent did we discover (at each percent coverage) (downward sloping)
erx = erexnames[which(edger.results$FDR < 0.05)]
length(erx)==length(unique(erx)) #TRUE
sum(erx %in% myinfo.updated$ex.names) #103, good...
yvar2 = NULL
for(i in 1:length(x)){
  countvar = NULL
  for(j in 1:length(erx)){
    derf.ind = which(myinfo.updated$ex.names==erx[j])
    if(length(derf.ind)==0){countvar[j] <- 0; next}
    if(length(derf.ind)>1){print("nonunique exon.")}
    countvar[j] = ifelse(myinfo.updated$qvals.ex[derf.ind]<0.05 & myinfo.updated$ex.pct[derf.ind]>=x[i],1,0)
  }
  yvar2[i] = sum(countvar)/length(countvar)
}

dex = deexnames[which(deseq.table$padj < 0.05)]
length(dex)==length(unique(dex)) #TRUE
sum(dex %in% myinfo.updated$ex.names) #107, good...
yvar3 = NULL
for(i in 1:length(x)){
  countvar = NULL
  for(j in 1:length(dex)){
    derf.ind = which(myinfo.updated$ex.names==dex[j])
    if(length(derf.ind)==0){countvar[j] <- 0; next}
    if(length(derf.ind)>1){print("nonunique exon.")}
    countvar[j] = ifelse(myinfo.updated$qvals.ex[derf.ind]<0.05 & myinfo.updated$ex.pct[derf.ind]>=x[i],1,0)
  }
  yvar3[i] = sum(countvar)/length(countvar)
}

setEPS() #[SUPPLEMENTARY FIGURE 2]
postscript("suppfig2.eps")
plot(x,yvar2, type="l",col="blue", lwd=2,ylim=c(0.3,1),xlab="Minimum Percentage of Exon Covered",ylab="Percentage of Identify-then-Annotate Exons Discovered" )
points(x,yvar3,col="red", type="l",lwd=2)
legend("topright",col=c("blue","red"), c("EdgeR", "DESeq"), lty=c(1,1), lwd=c(2,2))
dev.off()



# (3) autocorrelation plot: (code inspired by jeff's autocorrelation plot code)
library(Genominator)
newcounts = ExpData("Y-tophat-revised.db", tablename="chrY")
sex = c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1)
tmp = newcounts[1:1e5,which(sex==1)+1]
tmp = tmp - rowMeans(tmp)

aa = matrix(NA,nrow=9,ncol=51)
for(i in 1:9){
  aa[i,] = acf(tmp[,i])$acf
}

setEPS() #[SUPPLEMENTARY FIGURE 3]
postscript("autocor.eps")
plot(1:51,aa[1,],type="l",lwd=1.5,col="black",ylim=c(0,1),ylab="Correlation",xlab="Distance between features (base-pairs)")
for(i in 2:9){lines(1:51,aa[i,],type="l",lwd=1.5,col="black")}
lines(colMeans(aa)[2]^(0:50),col="red",lwd=3)
legend("bottomleft", col=c("black","red"), lwd=c(1.5, 3), c("observed correlation","AR(1) prediction"))
dev.off()


# t statistics histogram:
mp = getParams(tt)
mp

setEPS() #[SUPPLEMENTARY FIGURE 4]
postscript("ytdist.eps")
hist(tt,col="gray70",breaks=100,xlim=c(-5,10), xlab="moderated t statistics from Y chromosome", main="Test statistic distribution", freq=F, ylim=c(0,0.45))
xax = seq(-5,10,by=0.01)
lines(xax,mp$stateprobs[2]*dnorm(xax,mean=mp$params$mean[2],sd=mp$params$sd[2]),col="black",lwd=3)
lines(xax,mp$stateprobs[3]*dnorm(xax,mean=mp$params$mean[3],sd=mp$params$sd[3]),col="red",lwd=3)
lines(xax,mp$stateprobs[4]*dnorm(xax,mean=mp$params$mean[4],sd=mp$params$sd[4]),col="green",lwd=3)
legend("topright",lwd=c(3,3,3),col=c("black","red","green"),c("equally expressed","overexpressed-men","overexpressed-women"))
dev.off()


## p-value histograms from the simulated data (cuffdiff vs. us)
# us:
pvals.derfinder = pvals
load("pvals100.rda")  ## this involves a whole bunch of other code also... contact me if you want it, or just the rda file
pvals.sim = pvals
pvals = pvals.derfinder
rm(pvals.derfinder)
hist(pvals.sim, col="gray70",breaks=30, xlab="p values", main="DER Finder")

# cuffdiff
iso = read.table("cuffdiff_simulated/isoform_exp.diff",header=TRUE) # I can send you this also, or point you to it on the cluster
tail(iso)
hist(iso$p_value[iso$status=="OK"], col="gray70", breaks=30, xlab="p values", main="Cuffdiff")

setEPS(width=12, height=6) #[SUPPLEMENTARY FIGURE 5]
postscript("pairedhists.eps")
par(mfrow=c(1,2))
hist(pvals.sim, col="gray70",breaks=30, xlab="p values", main="DER Finder")
hist(iso$p_value[iso$status=="OK"], col="gray70", breaks=30, xlab="p values", main="Cuffdiff")
dev.off()

min(iso$q_value)



## getting some of the derfinder results for different read coverages:
## this involves some extra simulation code
#(1) load
load("regions.merged100.rda")
regions.merged100 = regions.merged
load("regions.merged50.rda")
regions.merged50 = regions.merged
load("regions.merged25.rda")
regions.merged25 = regions.merged
pvals.100 = pvals.sim
rm(pvals.sim)
pvals.derfinder = pvals
load("pvals50.rda")
pvals.50 = pvals
load("pvals25.rda")
pvals.25 = pvals
pvals = pvals.derfinder
rm(pvals.derfinder)


# plot:
regions_withp = data.frame(regions.merged100, pvals = pvals.100)
which_de = which((regions_withp$start + trunc(regions_withp$length/2)) %in% c(debps1, debps2))
null_regions = regions_withp[-which_de,]
pdf('null_pvals.pdf')
hist(null_regions)
dev.off()


#(2) get results
regions100 = data.frame(regions.merged100, pvals=pvals.100, qvals=p.adjust(pvals.100,method="fdr"))
regions50 = data.frame(regions.merged50, pvals=pvals.50, qvals=p.adjust(pvals.50,method="fdr"))
regions25 = data.frame(regions.merged25, pvals=pvals.25, qvals=p.adjust(pvals.25,method="fdr"))
length(which(regions100$qvals<0.05)) #435
length(which(regions50$qvals<0.05)) #433
length(which(regions25$qvals<0.05)) #0
min(regions25$qvals, na.rm=TRUE)
ders100 = subset(regions100,state==3|state==4)
ders50 = subset(regions50,state==3|state==4)
ders25 = subset(regions25,state==3|state==4)
sig100 = subset(regions100,qvals<0.05)
sig50 = subset(regions50, qvals<0.05)
sig25 = subset(regions25, qvals<0.05)

# load in the transcripts and DE transcripts:
### helper functions from Rhelp board (https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html)
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = " ", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer","character", "character", "character",     "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

de1 = gffRead("de-tx1.gtf") #the state 3's
de2 = gffRead("de-tx2.gtf") #the state 4's
nonde = gffRead("genes-clean-small.gtf") 
de1$tx <- getAttributeField(de1$attributes, "transcript_id", attrsep="; ")
de2$tx <- getAttributeField(de2$attributes, "transcript_id", attrsep="; ")
nonde$tx = getAttributeField(nonde$attributes, "transcript_id", attrsep="; ")
nonde = subset(nonde, !(tx %in% union(de1$tx,de2$tx)))
de1gr = GRanges(seqnames=Rle(de1$seqname), ranges=IRanges(start=de1$start, end=de1$end), strand=de1$strand)
de2gr = GRanges(seqnames=Rle(de2$seqname), ranges=IRanges(start=de2$start, end=de2$end), strand=de2$strand)
nondegr = GRanges(seqnames=Rle(nonde$seqname), ranges=IRanges(start=nonde$start, end=nonde$end), strand=nonde$strand)

alltx = gffRead("genes-clean-small.gtf")
alltx$length = alltx$end-alltx$start+1
alltx$tx = getAttributeField(alltx$attributes, "transcript_id", attrsep="; ")
txlens = NULL
for(tx in unique(alltx$tx)){
  txlens[which(unique(alltx$tx)==tx)] <- sum(alltx$length[alltx$tx==tx])
}

median(txlens) #1240 (dist is right skewed)
# Experiment 1: I generated 400 76bp reads/tx (200 pairs of reads), assume each tx has ~1240 bps
# --> so 400*76=30400 bps/tx --> 30400/1240 = 24x coverage.
# Experiment 2: 200 76bp reads/tx (100 pairs of reads)
# --> 200*76 = 15200 --> 15200/1240 = 12x coverage
# Experiment 3 = 6x coverage.

qcuts = seq(0.01,1,by=0.01)
regions100$qvals[is.na(regions100$qvals)] <- 1
regions50$qvals[is.na(regions50$qvals)] <- 1
regions25$qvals[is.na(regions25$qvals)] <- 1
sens100 = sens2.100 = spec100 = NULL
sens50 = sens2.50 = spec50 = NULL
sens25 = sens2.25 = spec25 = NULL

for(qcut in qcuts){
  i=which(qcuts==qcut)
  sig100 = subset(regions100,qvals<=qcut)
  sig50 = subset(regions50, qvals<=qcut)
  sig25 = subset(regions25, qvals<=qcut)
  
  neg100 = subset(regions100,qvals>=qcut)
  neg50 = subset(regions50,qvals>=qcut)
  neg25 = subset(regions25,qvals>=qcut)
  
  derf100gr1 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig100$start[sig100$state==3], end=sig100$end[sig100$state==3]))
  derf100gr2 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig100$start[sig100$state==4], end=sig100$end[sig100$state==4]))
  derf100.neg = GRanges(seqnames=Rle(22), ranges=IRanges(start=neg100$start, end=neg100$end))
  
  derf50gr1 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig50$start[sig50$state==3], end=sig50$end[sig50$state==3]))
  derf50gr2 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig50$start[sig50$state==4], end=sig50$end[sig50$state==4]))
  derf50.neg = GRanges(seqnames=Rle(22), ranges=IRanges(start=neg50$start, end=neg50$end))
  
  aa = try(GRanges(seqnames=Rle(22), ranges=IRanges(start=sig25$start[sig25$state==3], end=sig25$end[sig25$state==3])),silent=TRUE)
  ab = try(GRanges(seqnames=Rle(22), ranges=IRanges(start=sig25$start[sig25$state==4], end=sig25$end[sig25$state==4])),silent=TRUE)
  derf25.neg = GRanges(seqnames=Rle(22), ranges=IRanges(start=neg25$start, end=neg25$end))
  
  
  ### 24x dataset:
  j = countOverlaps(de1gr, derf100gr1)
  jj = countOverlaps(de2gr, derf100gr2)
  jjj = c(j,jj) #my supply of sensible variable names is running really low...sorry guys.
  sens100[i] = sum(jjj>0)/length(jjj)
  tmp.df = data.frame(tx=c(de1$tx,de2$tx), ol=jjj)
  pertxnum = NULL
  for(tx in unique(tmp.df$tx)){
    pertxnum[which(unique(tmp.df$tx)==tx)] <- sum(tmp.df$ol[tmp.df$tx==tx])
  }
  sens2.100[i] = sum(pertxnum>0)/length(pertxnum)
  m = countOverlaps(nondegr, derf100.neg)
  spec100[i] = sum(m>0)/length(m)
  
  ### 12x dataset:
  j = countOverlaps(de1gr, derf50gr1)
  jj = countOverlaps(de2gr, derf50gr2)
  jjj = c(j,jj) 
  sens50[i] = sum(jjj>0)/length(jjj)
  tmp.df = data.frame(tx=c(de1$tx,de2$tx), ol=jjj)
  pertxnum = NULL
  for(tx in unique(tmp.df$tx)){
    pertxnum[which(unique(tmp.df$tx)==tx)] <- sum(tmp.df$ol[tmp.df$tx==tx])
  }
  sens2.50[i] = sum(pertxnum>0)/length(pertxnum)
  m = countOverlaps(nondegr, derf50.neg)
  spec50[i] = sum(m>0)/length(m)
  
  ### 6x dataset:
  if(class(aa)!="try-error" & class(ab)!="try-error"){
    derf25gr1 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig25$start[sig25$state==3], end=sig25$end[sig25$state==3]))
    derf25gr2 = GRanges(seqnames=Rle(22), ranges=IRanges(start=sig25$start[sig25$state==4], end=sig25$end[sig25$state==4]))
    j = countOverlaps(de1gr, derf25gr1)
    jj = countOverlaps(de2gr, derf25gr2)
    jjj = c(j,jj) #my supply of sensible variable names is running really low...sorry guys.
    sens25[i] = sum(jjj>0)/length(jjj)
    tmp.df = data.frame(tx=c(de1$tx,de2$tx), ol=jjj)
    pertxnum = NULL
    for(tx in unique(tmp.df$tx)){
      pertxnum[which(unique(tmp.df$tx)==tx)] <- sum(tmp.df$ol[tmp.df$tx==tx])
    }
    sens2.25[i] = sum(pertxnum>0)/length(pertxnum)
    m = countOverlaps(nondegr, derf25.neg)
    spec25[i] = sum(m>0)/length(m)
  }
}

setEPS() #[SUPPLEMENTARY FIGURE 6]
postscript("roc.eps")
plot(c(0,1-spec100), c(0,sens100), col="dodgerblue2", lty=2, type="l", lwd=2, xlim=c(0,0.25),ylim=c(0,1), xlab="false positive rate", ylab="true positive rate", main="DER Finder ROC Curves")
lines(c(0,1-spec100), c(0,sens2.100), col="dodgerblue2", lwd=2)
lines(c(0,1-spec50), c(0,sens50), col="orange", lty=2, lwd=2)
lines(c(0,1-spec50), c(0,sens2.50), col="orange", lwd=2)
lines(c(0,1-spec25), c(0,sens25), col="green3", lty=2, lwd=2)
lines(c(0,1-spec25), c(0,sens2.25), col="green3", lwd=2)
legend("bottomright", lty=c(1,1,1,1,2), col=c("dodgerblue2","orange","green","black","black"), lwd=rep(2,5), c("24x","12x","6x","transcript-level","feature-level"))
dev.off()






