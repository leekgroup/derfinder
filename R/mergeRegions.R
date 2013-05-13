## functions for merging regions that are less than <cutoff> bp apart
## regions = data frame of regions, start, end, state.  (getRegions's "$states" output)
## cutoff:  how close together do regions need to be in order to be merged?

#'Merge close regions
#'
#'Regions of the same state that are fewer than 5bp apart from each other are
#'merged.
#'
#'
#'@param regions data frame of regions, usually the \code{$state} return of
#'\code{getRegions}
#'@param cutoff How far apart do regions need to be in order to *not* be
#'merged?  (Defaults to 6, meaning that regions fewer than 5bp apart *are*
#'merged).
#'@return data frame (comparable to \code{regions} input) with close regions
#'merged.
#'@author Alyssa Frazee
#'@export
#'@seealso \code{\link{getRegions}}

mergeRegions <- function(regions, cutoff=6){
	ind = which(regions$length<cutoff)[1]
	
	while(ind<dim(regions)[1]){
		#print(ind)
		anchor.down.ind = findclosest.behind(ind,regions)
		anchor.up.ind = findclosest.ahead(ind,regions)
	
		if(ind==1){
			ind = which(regions$length<6)[2]
			if(is.na(ind)) break
			next
		}
		
		# if the anchors don't match, leave the small region alone.
		if(regions$state[anchor.down.ind]!=regions$state[anchor.up.ind]){
			smalls = which(regions$length<6)
			ind = smalls[which(smalls>ind)][1]
			if(is.na(ind)) break
			next	
		}
		
		# if the anchors match, merge:
		if(regions$state[anchor.down.ind]==regions$state[anchor.up.ind]){
			regions$end[ind-1] = regions$end[ind+1]
			
			if("mean.t" %in% names(regions)){
			regions$mean.t[ind-1] = (regions$length[ind-1]*regions$mean.t[ind-1]+regions$length[ind]*regions$mean.t[ind]+regions$length[ind+1]*regions$mean.t[ind+1])/(regions$length[ind-1]+regions$length[ind]+regions$length[ind+1])
			}
			if("mean.fold.change" %in% names(regions)){
			regions$mean.fold.change[ind-1] = (regions$length[ind-1]*regions$mean.fold.change[ind-1]+regions$length[ind]*regions$mean.fold.change[ind]+regions$length[ind+1]*regions$mean.fold.change[ind+1])/(regions$length[ind-1]+regions$length[ind]+regions$length[ind+1])
			}
			regions$length[ind-1] = regions$length[ind-1]+regions$length[ind]+regions$length[ind+1]
			regions = regions[-c(ind,ind+1),]
			smalls = which(regions$length<6)
			ind = smalls[which(smalls>=ind)][1]
			#print("made 1 change")
			if(is.na(ind)) break
		} #end merging
} #end while loop
	return(regions)
} #end function

