## arguments - genome: UCSC genome name (in string format)
## output - prints the tables that can be downloaded (i.e., supplied to getAnnotation) for the supplied genome.
supportedTables <-
function(genome){
	message("getting supported tables. may take several minutes.")
	require(GenomicFeatures)
	genome.tracks = supportedUCSCFeatureDbTracks(genome)
	intersect(genome.tracks,rownames(supportedUCSCtables()))
	}
