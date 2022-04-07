# Function:
#	boxnetworks: creates box (location) specific networks
# Arguments:
#	record: data frame with movement records, supposed to contain components 'ID' (individual ID), 'Time' (time of entrance/exit), 'Box' (identifier of the location) and 'Enter' (logical vector, where TRUE means entrance and FALSE exit)
#	normalize: whether to normalize link weights so they express proportions of timespan of the record
# Value:
#	array of size n x n x k (n = no. of individuals, k = no. of boxes), comprising k matrices corresponding to box (=location) specific networks 

boxnetworks <- function(record, normalize=TRUE) {
	record$ID <- factor(record$ID)
	record$Time <- as.numeric(record$Time)
	record$Box <- factor(record$Box)
	record$Enter <- factor(ifelse(record$Enter, "en", "ex"))
	daysplit <- split(record[,c("Time","Box","Enter")], record$ID)
	daysplit <- daysplit[order(names(daysplit))]
	daysplit <- unlist(lapply(daysplit, function(x) split(x[,c("Time","Enter")], x$Box)), recursive=FALSE)
	intervals <- lapply(daysplit, function(x) do.call(cbind, split(x$Time, x$Enter)))
	nr <- lapply(intervals, nrow)
	nr <- unlist(ifelse(sapply(nr, length), nr, 0))
	id <- rep(names(intervals), nr)
	idsplit <- strsplit(id, split="\\.")
	ind <- sapply(idsplit, "[", 1)
	box <- factor(sapply(idsplit, "[", 2))
	mat <- do.call(rbind, intervals)
	exen <- outer(mat[,"ex"], mat[,"en"], "-")
	exex <- outer(mat[,"ex"], mat[,"ex"], "-")
	exex <- exex * (exex > 0)
	enen <- outer(mat[,"en"], mat[,"en"], "-")
	enen <- enen * (enen > 0)
	exen <- exen - (enen + exex)
	exen[exen < 0] <- 0
	boxparts <- split.data.frame(exen, box)
	for (i in seq_along(boxparts)) {
		boxparts[[names(boxparts)[i]]] <- as.matrix(boxparts[[i]][,box==names(boxparts)[i]])		
	}
	indparts <- split(ind, box)
	boxes <- array(0, dim=c(nlevels(record$ID),nlevels(record$ID),nlevels(record$Box)), dimnames=list(levels(record$ID),levels(record$ID),levels(record$Box)))
	for (n in levels(record$Box)) {
		lst <- lapply(lapply(split.data.frame(boxparts[[n]], indparts[[n]]), t), split.data.frame, f=indparts[[n]])
		present <- match(names(lst), levels(record$ID))
		boxes[present, present, n]  <- do.call(rbind, lapply(lst, function(x) sapply(x, sum)))
	}
	if (isTRUE(normalize)) {
		boxes <- boxes / diff(range(record$Time))		
	}
	return(boxes)
}


# Function:
#	sumnetwork: sums box specific networks into a single network
# Arguments:
#	boxes: array with box (=location) specific networks, output of 'boxnetworks'
#	diagonal: whether to fill in diagonal with times the individual spent alone in some of the boxes
# Value:
#	matrix of size n x n (n = no. of individuals), which specifies the social network

sumnetwork <- function(boxes, diagonal=FALSE) {
	net <- apply(boxes, c(1,2), sum)
	empty <- which(apply(net == 0, 1, all))
	if (length(empty) > 0) {
		net <- net[-empty,-empty]		
	}
	if (isFALSE(diagonal)) {
		diag(net) <- 0
	}
	return(net)
}


# Function:
#	make_multilayer_network: a wrapper function which applies 'boxnetworks' & 'sumnetwork' over the list of movement records
# Arguments:
#	records: list of data frames with movement records (with components 'ID', 'Time', 'Box' and 'Enter'), each corresponding to a single time layer
#	normalize: whether to normalize link weights so they express proportions of timespan of the record
#	diagonal: whether to fill in diagonal with times the individual spent alone in some of the boxes
# Value:
#	list of matrices specifying layer-specific social networks

make_multilayer_network <- function(records, normalize=TRUE, diagonal=FALSE) {
	return(lapply(lapply(records, boxnetworks, normalize=normalize), sumnetwork, diagonal=diagonal))
}
