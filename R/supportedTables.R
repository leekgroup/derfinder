## arguments - genome: UCSC genome name (in string format)
## output - prints the tables that can be downloaded (i.e., supplied to getAnnotation) for the supplied genome.


#'Print list of supported (downloadable) tables for a given genome
#'
#'To be used with \code{getAnnotation}: provides a list of available tables
#'from UCSC for any given genome.
#'
#'
#'@param genome Character string giving the genome of interest.  Available
#'genomes can be seen with \code{rtracklayer:::ucscGenomes()[,"db"]}.
#'@return prints a list of available tables for \code{genome}.
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getAnnotation}}
#'@examples
#'
#' \dontrun{
#' supportedTables("mm10")
#' mouse.exons = getAnnotation("mm10","refGene") #refGene appears in printed output of supportedTables("mm10").
#' }
#'
supportedTables <- function(genome){
	message("Getting supported tables; may take several minutes.")
	require("GenomicFeatures")
	genome.tracks = supportedUCSCFeatureDbTracks(genome)
	intersect(genome.tracks,rownames(supportedUCSCtables()))
}
