## arguments - 
## --genome: UCSC genome name (in string format)
## --tablename: the table to download from UCSC.  Can see supported tables from a given genome using supportedTables()
## --genes: if TRUE, exons in the output are labeled by gene, otherwise exons are labeled by transcript
## --verbose: if TRUE, messages are printed as the function runs.
## output - if genes is TRUE, a data frame containing one row per exon with exon information for the given genome
## output - if genes is FALSE, a list with elements $datf (the exon table, same as if genes is TRUE)
## and $transcript.names, a named character vector where the names are internal transcript IDs and the elements are external transcript names

getAnnotation <-
function(genome, tablename, genes = TRUE, verbose = TRUE){
	require(GenomicFeatures)
	require(rtracklayer)
	a <- try(makeTranscriptDbFromUCSC(genome=genome, tablename=tablename))
	if(class(a)=="try-error") stop("Problem accessing requested UCSC annotation - likely there is a problem with genome or tablename arguments. Use ucscGenomes() to see acceptable genomes; use supportedTables(genome) to see acceptable tablenames for your genome.")
	
	if(verbose) show(a)
	if(genes){
		if(verbose) print("Labeling exons by gene...")
		grl <- exonsBy(a,"gene")		
		if(length(grl)==0) stop("cannot list exons by gene using this table. please try another table or specify genes = FALSE to list exons by transcript.")
		datf <- IRanges:::as.data.frame(grl)
		names(datf)[1] <- "gene"
	}
	if(!genes){
		if(verbose) print("Labeling exons by transcript...")
		grl <- exonsBy(a,"tx")
		datf <- IRanges:::as.data.frame(grl)
		names(datf)[1] <- "transcript"
		nametable <- id2name(a,feature.type = "tx")
	}
	names(datf)[2] <- "chr"
	datf <- datf[,1:7]
	if(genes) return(datf)
	return(list(exon.table = datf, transcript.names = nametable))
}
