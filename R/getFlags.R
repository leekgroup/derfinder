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
#'@param pctcut What percentage of an exon should a region overlap in order to
#'call that exon differentially expressed?  Default 0.8 (meaning 80\%).
#'@return List with elements having length equal to the number of rows in
#'\code{regions}:
#'\item{flags}{Event indicated by corresponding region}
#'\item{flag.info}{Genomic features (exons) associated with genomic events}
#'\item{percent.exon}{The percent of the exon (in flag.info) overlapped by this region}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getAnnotation}}

getFlags <- function (regions, exons, chromosome, pctcut = 0.8){
    ## Appeasing R CMD check
    ## More info at http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    chr <- state <- NULL
    require(GenomicRanges)
    
    exons = subset(exons, chr == chromosome)
    exons = exons[order(exons$start), ]

    regions = subset(regions, chr == chromosome)

    stopifnot(length(unique(exons$chr))==1, 
        length(unique(regions$chr)) == 1,
        unique(exons$chr) == unique(regions$chr))

    exgr = GRanges(seqnames = Rle(exons$chr), 
        ranges = IRanges(start = exons$start, end = exons$end))
    
    candidates = subset(regions, state==3 | state==4)

    regionsgr = GRanges(seqnames=Rle(candidates$chr), 
        ranges = IRanges(start = candidates$start, end = candidates$end))

    overlaps = findOverlaps(regionsgr, exgr)

    ex_by_region = split(subjectHits(overlaps), queryHits(overlaps))

    annotate_name = function(i){
        if(i %in% names(ex_by_region)){
            inds = ex_by_region[[ which(names(ex_by_region)==i) ]]
            return(exons$exon_id[inds])
        }else{
            return(NA)
        }
    }

    annotate_pct = function(i){
        if(i %in% names(ex_by_region)){
            inds = ex_by_region[[ which(names(ex_by_region)==i) ]]
            pcts = rep(NA, length(inds))
            rpos = c(candidates$start[i]:candidates$end[i])
            for(exi in 1:length(inds)){
                expos = c(exons$start[inds[exi]]:exons$end[inds[exi]])
                pcts[exi] = length(intersect(rpos, expos)) / length(expos)
            }
            return(pcts)
        }else{
            return(NA)
        }
    }

    flag.info = lapply(1:nrow(candidates), annotate_name)
    percent.exon = lapply(1:nrow(candidates), annotate_pct)
    flags = rep(NA, nrow(candidates))
    flags[is.na(flag.info)] = "novel"

    has_enough_ol = function(i){
        if(!is.na(percent.exon[[i]][1])){
            return(any(percent.exon[[i]]>pctcut))
        }else{
            return(FALSE)
        }
    }

    potential_de = sapply(1:nrow(candidates), has_enough_ol, USE.NAMES=FALSE)
    flags[potential_de] = "DE exon(s)"

    return(list(flags = flags, flag.info = flag.info, percent.exon = percent.exon))

}



