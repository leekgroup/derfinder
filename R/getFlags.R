## getFlags()
## arguments:
## regions:  data frame of genomic regions (usually the $states output of getRegions).  must contain "state," "start", "end"
## exons:  data frame containing exon information.  see getAnnotation().  must contain columns "start", "end", "exon_id", "gene", and "chr"
## chromosome:  which chromosome are you analyzing?  must match with chromosome names in exons$chr
## nonDE:  if TRUE, get p-values for regions of state 2 rather than 3 and 4. (rarely used)
## pctcut:  what percent of an exon does a region need to overlap in order to have that region make a statement about the exon?  default 0.8 (80%)



#'Flag regions with genomic events of interest
#'
#'Connects results of tornado pipeline to existing annotation as a way to
#'determine which regions may point to interesting biological events.
#'
#'
#'@param regions data frame of regions to consider, usually the \code{$states}
#'output of \code{getRegions}
#'@param exons data frame containing annotated exons you would like to
#'consider.  Should have columns
#'\code{chr},\code{start},\code{end},\code{exon_id}, and \code{gene}.  Can be
#'created with \code{\link{getAnnotation}}.
#'@param chromosome Chromosome you're considering.  Currently you can only do
#'one chromosome at a time.
#'@param nonDE Should we give you flags for only the equally expressed regions?
#'(Usually we just want flags for differentially expressed regions, so in
#'general this is FALSE).
#'@param pctcut What percentage of an exon should a region overlap in order to
#'call that exon differentially expressed?  Default 0.8 (meaning 80%).
#'@return List with elements having length equal to the number of rows in
#'\code{regions}:
#'\item{flags}{Event indicated by corresponding region}
#'\item{flag.info}{Genomic features (exons) associated with genomic events}
#'\item{percent.exon}{The percent of the exon (in flag.info) overlapped by this region}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getAnnotation}}

getFlags <- function (regions, exons, chromosome, nonDE = FALSE, pctcut = 0.8){
     exons = subset(exons, chr == chromosome)
     exons = exons[order(exons$start), ]
     if(nonDE){ candidates = subset(regions, state==2) }
     if(!nonDE){ candidates = subset(regions, (state == 3 | state == 4))}
     cand.pos = NULL
     for (i in 1:dim(candidates)[1]) {
         cand.pos = append(cand.pos, c(candidates$start[i]:candidates$end[i]))
     }
     inexon = rep(0, length(cand.pos))
     whichexon = whichgene = rep(NA, length(cand.pos))
     for (i in 1:length(inexon)) {
          ind = findInterval(cand.pos[i], exons$start)
          if (ind == 0) next
          if (cand.pos[i] %in% c(exons$start[ind]:exons$end[ind])) {
              inexon[i] <- 1
              whichexon[i] <- exons$exon_id[ind]
              whichgene[i] <- exons$gene[ind]
          }
      }
      percent.exon = list()
      flags = rep(NA, dim(candidates)[1])
      flag.info = list()
      for (i in 1:dim(candidates)[1]) {
          info = inexon[which(cand.pos == candidates$start[i]):which(cand.pos == candidates$end[i])]
          name.info = whichexon[which(cand.pos == candidates$start[i]):which(cand.pos == candidates$end[i])]
          cand.pos.this = cand.pos[which(cand.pos==candidates$start[i]):which(cand.pos==candidates$end[i])]
          if (length(unique(name.info[!is.na(name.info)])) == 0) {
              percent.exon[[i]] = 0
              if(!nonDE) flags[i] = "novel DE region"
              if(nonDE) flags[i] = "novel region"
              flag.info[[i]] = "novel"
          }
          if (length(unique(name.info[!is.na(name.info)])) == 1) {
              myexon = exons[which(exons$exon_id == unique(name.info[!is.na(name.info)])),]
              if (dim(myexon)[1] > 1) {
                  if (length(unique(myexon$start)) > 1 | length(unique(myexon$end)) > 1) {
                         message(paste(myexon$exon_id[1], "is a problematic exon id:"))
                         stop("exon IDs not unique: please check annotation.")
                  }
                  myexon = myexon[1, ]
              }
              percent.exon[[i]] = sum(info)/(myexon$end - myexon$start + 1)
              flag.info[[i]] = unique(name.info[!is.na(name.info)])
              if (percent.exon[[i]] > pctcut) {
                 if(!nonDE) flags[i] = "DE exon"
                 if(nonDE) flags[i] = "expressed exon"
              } 
           }
           if (length(unique(name.info[!is.na(name.info)])) > 1) {
              indlist = split(1:length(name.info), name.info)
              pctvec = NULL
              this.name.info = NULL
              for (j in 1:length(indlist)) {
                 inds = indlist[[j]] 
                 this.name.info[j] = names(indlist)[j] 
                 myexon = exons[which(exons$exon_id == this.name.info[[j]]),]
                 if (dim(myexon)[1] > 1) {
                    if (length(unique(myexon$start)) > 1 | length(unique(myexon$end)) > 1) {
                        message(paste(myexon$exon_id[1]), "is a problematic exon id:")
                        stop("exon IDs not unique: please check annotation")
                      }
                    myexon = myexon[1, ]
                 }
                 #pctvec[j] = sum(info[inds])/(myexon$end - myexon$start + 1) 
                 exinfo = c((myexon$start):(myexon$end))
                 inregion = ifelse(exinfo %in% cand.pos.this,1,0)
                 pctvec[j] = sum(inregion)/length(inregion)
              }
              percent.exon[[i]] = pctvec
              flag.info[[i]] = this.name.info
              if (sum(pctvec > pctcut) == 1) {
                 if(!nonDE) flags[i] = "DE exon"
                 if(nonDE) flags[i] = "expressed exon"
              }
              if (sum(pctvec > pctcut) > 1) {
                 if(!nonDE) flags[i] = "DE exons"
                 if(nonDE) flags[i] = "expressed exons"
                    
               }
            }
      }
            return(list(flags = flags, flag.info = flag.info, percent.exon = percent.exon))
      }
