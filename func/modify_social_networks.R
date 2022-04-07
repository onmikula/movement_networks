# Function:
#	socnet_merge: merges social networks
# Arguments:
#	snets: list of matrices encoding social networks
#	which: numeric, which social networks (from a list) to merge
#	what: a kind of merger, union (default) or intersect
#	times: list of starting and ending times of the network (numeric vectors of length two) or of their durations (numeric scalars)
# Value:
#	adjacenecy matrix of the summary network 

socnet_merge <- function(snets, which=seq_along(snets), what=c("union", "intersect"), times=NULL) {
	snets <- snets[which]
	if (!is.null(time)) {
		times <- times[which]
	} else {
		times <- as.list(rep(1, length(which)))
	}
	if (length(times[[1]]) == 2) {
		times <- lapply(times, diff)
	}
	total <- sum(unlist(times))
	snets <- Map("*", snets, times)
	if (what[1] == "union") {
		nams <- sort(unique(unlist(lapply(snets, rownames))))
	} else if (what[1] == "intersect") {
		nams <- sort(unique(Reduce(intersect, lapply(snets, rownames))))
	}
	merged <- matrix(0, length(nams), length(nams), dimnames=list(nams, nams))
	pairs <- matrix(nams[which(upper.tri(merged), arr.ind=TRUE)], (length(nams) * (length(nams) - 1)) / 2, 2)
	for (i in seq(nrow(pairs))) {
		ijpair <- 0
		for (j in seq_along(snets)) {
			if (all(pairs[i,] %in% rownames(snets[[j]]))) {
				ijpair <- ijpair + snets[[j]][pairs[i,1], pairs[i,2]]
			}
		}
		merged[pairs[i,1], pairs[i,2]] <- merged[pairs[i,2], pairs[i,1]] <- ijpair / total
	}
	return(merged)
}


# Function:
#	socnet_subset: subsets social network
# Arguments:
#	snets: list of adjacency matrices (encoding social networks) or a single such matrix
#	subsets: list of character vectors specifying subsets of individuals to be retained in the corresponding networks (if they are present there); if 'snets' is a single matrix it has to be a single character vector
# Value:
#	adjacenecy matrices of the subset networks

socnet_subset <- function(snets, subsets) {
	single <- !is.list(snets) & is.matrix(snets)
	if (single) {
		snets <- list(snets)
		subsets <- list(subsets)
	}
	for (i in seq_along(snets)) {
		retain <- intersect(rownames(snets[[i]]), subsets[[i]])
		snets[[i]] <- snets[[i]][retain, retain, drop=FALSE]
	}
	if (single) {
		snets <- snets[[1]]
	}
	return(snets)
}

