# Function:
#	writeMultilayer: writes a multilayer network into a file in Pajek format accepted by Infomap 
# Arguments:
#	layers: list of matrices describing layer-specific social networks
#	file: file name of the output
#	perlayer: whether to write a separate file for each layer
#	return: whether to return data frame

writeMultilayer <- function(layers, file, perlayer=FALSE, return=FALSE) {
	vertices <- sort(unique(unlist(lapply(layers, rownames))), decreasing=FALSE)
	intralinks <- vector("list", length(layers))
	nonempty <- rep(TRUE, length(layers)) 
	for (i in seq_along(layers)) {
		lay <- which(layers[[i]] != 0, arr.ind=TRUE)
		lay <- rownames(layers[[i]])[as.numeric(lay)]
		if (length(lay) > 0) {
			intralinks[[i]] <- matrix(match(lay, vertices), length(lay) / 2, 2)
			intralinks[[i]] <- cbind(i, intralinks[[i]], layers[[i]][layers[[i]] != 0])
		} else {
			nonempty[i] <- FALSE
		}
	}
	if (isFALSE(perlayer)) {
		intralinks <- do.call(rbind, intralinks)
		if (all(nonempty)) {
			nonempty <- NULL
		} else {
			nonempty <- paste("#layers:", paste(which(nonempty), collapse=", "))
			intralinks[,1] <- match(intralinks[,1], sort(unique(intralinks[,1])))
		}
		present <- sort(unique(as.numeric(intralinks[,2:3])))
		vertices <- vertices[present]
		intralinks[,2:3] <- match(intralinks[,2:3], present)
		intralinks <- apply(intralinks, 1, paste, collapse=" ")
		vmatrix <- cbind(seq_along(vertices), paste("\"", vertices, "\"", sep=""), "1.0")
		content <- c(paste("*Vertices", length(vertices)), apply(vmatrix, 1, paste, collapse=" "), "*Intra", "#layer node node [weight]", nonempty, intralinks)	
		writeLines(content, con=file, sep="\n")
	} else {
		intralinks <- intralinks[nonempty]
		stem <- sub("\\.[[:alpha:]]+$", "", file)
		extension <- sub(stem, "", file, fixed=TRUE)
		no <- formatC(seq_along(layers), format="d", flag=0, width=nchar(length(layers)))
		file <- paste0(paste(stem, paste0("layer", no), sep="_"), extension)	
		for (i in seq_along(intralinks)) {
			present <- sort(unique(as.numeric(intralinks[[i]][,2:3])))
			pvertices <- vertices[present]
			intralinks[[i]][,2:3] <- match(intralinks[[i]][,2:3], present)
			intralinks[[i]] <- apply(intralinks[[i]][,-1], 1, paste, collapse=" ")
			vmatrix <- cbind(seq_along(pvertices), paste("\"", pvertices, "\"", sep=""), "1.0")
			content <- c(paste("*Vertices", length(pvertices)), apply(vmatrix, 1, paste, collapse=" "), "*Edges", intralinks[[i]])	
			writeLines(content, con=file[i], sep="\n")
		}
	}
	if (isTRUE(return)) {
		return(file)
	}
}


# Function:
#	writePajek: writes a single-layer network into a file in Pajek format accepted by Infomap 
# Arguments:
#	mat: matrix describing social network
#	file: file name of the output

writePajek <- function(mat, file) {
	vertices <- apply(cbind(seq(nrow(mat)), paste("\"", rownames(mat), "\"", sep="")), 1, paste, collapse=" ")
	edges <- rbind(cbind(which(lower.tri(mat), arr.ind=TRUE), mat[lower.tri(mat)]), cbind(which(upper.tri(mat), arr.ind=TRUE), mat[upper.tri(mat)])) 
	edges <- apply(edges, 1, paste, collapse=" ")
	result <- c(paste("*Vertices", nrow(mat)), vertices, paste("*Edges", nrow(edges)), edges)
	writeLines(result, file)
}


# Function:
#	writeMultilayerClu: writes a multilayer network partition (= clustering of vertices) into a /clu file accepted by Infomap
# Arguments:
#	P: a matrix indicating cluster-membership of individuals (in rows) in the clusters (=modoules, communities...) across successive time layers (in columns)
#	file: file name of the output
#	perlayer: whether to write a separate file for each layer
#	return: whether to return data frame

writeMultilayerClu <- function(P, file, perlayer=FALSE, return=FALSE) {
	P <- as.matrix(P)
	clu <- which(!is.na(P), arr.ind=TRUE)[,2:1]
	clu <- cbind(clu, P[which(!is.na(P))])
	if (isTRUE(perlayer)) {
		clu <- split(as.data.frame(clu[,-1]), clu[,1])
		stem <- sub("\\.[[:alpha:]]+$", "", file)
		extension <- sub(stem, "", file, fixed=TRUE)
		no <- formatC(seq_along(clu), format="d", flag=0, width=nchar(length(clu)))
		file <- paste0(paste(stem, paste0("layer", no), sep="_"), extension)	
		for (i in seq_along(clu)) {
			clu[[i]] <- clu[[i]][order(rownames(clu[[i]])), 2]
			clu[[i]] <- c("# node cluster", paste(seq_along(clu[[i]]), clu[[i]]))
		}
		Map(writeLines, clu, file)		
	} else {
		clu <- c("# layer node cluster", apply(clu, 1, paste, collapse=" "))
		writeLines(clu, file)		
	}
	if (isTRUE(return)) {
		return(file)
	}
}


# Function:
#	extractCompressionRate: extracts compression rate from .clu output file of Infomap
# Arguments:
#	clu: the name of .clu file
# Value:
#	compression rate value

extractCompressionRate <- function(clu) {
	record <- readLines(clu, n=1)
	from_start <- regexpr("from codelength ", record)
	from_end <- regexpr(" in one level", record)
	from <- as.numeric(substring(record, from_start + attr(from_start,"match.length"), from_end - 1))
	to_start <- regexpr("to codelength ", record)
	to_end <- gregexpr(" in", record)
	to <- as.numeric(substring(record, to_start + attr(to_start,"match.length"), rev(to_end[[1]])[1] - 1))
	comprate <- from / to
	return(comprate)	
}
