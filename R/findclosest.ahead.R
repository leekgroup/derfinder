#'Find closest long region downstream
#'
#'Helper function for \code{mergeRegions}
#'
#'
#'@param ind index of regions data frame
#'@param regions the $states data frame, as returned by \code{getRegions}
#'@return index of closest large region downstream of a region that needs to be
#'merged
#'@note Not generally used alone - internal function for \code{mergeRegions}
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{mergeRegions}}
findclosest.ahead = function(ind,regions){
	if(ind==dim(regions)[1]){
		if(regions$length[ind]>5){return(ind)}
		if(regions$length[ind]<6){
			warning("no regions ahead with length > 5")
			return(ind)
		}
	}
	length.of.next = regions$length[ind+1]
	if(length.of.next>5) return(ind+1)
	if(length.of.next<6) findclosest.ahead(ind+1,regions)
}
