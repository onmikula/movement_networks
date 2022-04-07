# Function:
#	allocation: deciphers placement (allocation) of modules in boxes
# Arguments:
#	records: list of data frame with movement records, supposed to contain components 'ID' (individual ID), 'Time' (time of entrance/exit), 'Box' (identifier of the location) and 'Enter' (logical vector, where TRUE means entrance and FALSE exit)
#	modules: matrix indicating module membership, whose rownames must match individual IDs and no. of columns must match the length of 'records' (i.e. no. of time-layers)
# Value:
#	array of dimension boxes x modules x layers 

allocation <- function(records, modules) {
	if (is.data.frame(records)) {
		records <- list(records)
		modules <- as.matrix(modules)[,1,drop=FALSE]
	}
	for (i in seq_along(records)) {
		records[[i]]$Box <- factor(records[[i]]$Box)
	}
    k <- max(modules, na.rm=TRUE)

	mat <- setNames(vector("list", length(records)), names(records))
	for (i in seq_along(records)) {
		if (length(records[[i]]$ID) > 0) {
			ipart <- lapply(split(records[[i]], records[[i]]$ID), function(x) split(x, x$Box))
			imat <- matrix(, length(ipart), nlevels(records[[i]]$Box), dimnames=list(names(ipart),levels(records[[i]]$Box))) 
			for (j in seq_along(ipart)) {
				jpart <- sapply(ipart[[j]], function(x) sum(diff(matrix(x$Time, 2, nrow(x)/2))))
				imat[j,match(names(jpart),colnames(imat))] <- jpart
			}
			colnames(imat) <- levels(records[[i]]$Box)
			imat <- imat[rownames(imat) %in% rownames(modules),,drop=FALSE]
			grp <- modules[rownames(imat),i]
			mat[[i]] <- t(do.call(rbind, by(imat, grp, colSums)))
		}
	}
	boxes <- levels(records[[1]]$Box)[order(as.numeric(levels(records[[1]]$Box)))]
	spat <- array(0, c(nlevels(records[[1]]$Box), k, length(records)), dimnames=list(boxes, seq(k), names(records)))
	for (i in seq_along(records)) {
		if (length(records[[i]]$ID) > 0) {
			spat[,,i] <- mat[[i]][match(rownames(spat[,,i]),rownames(mat[[i]])), match(colnames(spat[,,i]),colnames(mat[[i]]))]
		}
	}
    return(spat)
}



# Function:
#	spatsepar: quantifies a degree of spatial separation between communities (modules) in a given time layer
# Arguments:
#	spat: an array of dimension boxes x modules x layers (produced by 'allocation'), which contains information about possession of boxes by different modules across a range of time layers 
# Value:
#	vector of spatial separation indices

spatsepar <- function(spat) {
	spat <- ifelse(is.na(spat), 0, spat)
	w <- apply(spat, c(1,3), sum, na.rm=TRUE)
	for (i in seq(dim(spat)[3])) {
		wi <- ifelse(w[,i] == 0, 0, 1 / w[,i])
		spat[,,i] <- diag(wi) %*% spat[,,i]
	}
	s <- suppressWarnings(apply(spat, c(1,3), max, na.rm=TRUE))
	s[s == -Inf] <- NA	
	w <- w %*% diag(1/colSums(w))
	sep <- ifelse(colSums(s) > 0, colSums(s * w, na.rm=TRUE), NaN)
	return(sep)
}


# Function:
#	boxpossession: quantifies a degree of spatial separation between communities (modules) in a given time layer
# Arguments:
#	spat: an array of dimension boxes x modules x layers (produced by 'allocation'), which contains information about possession of boxes by different modules across a range of time layers 
#	digits: precision of the output (no. of decimal places)
# Value:
#	array of the same dimension and structure as 'spat', but containing proportional occupancy of boxes by different communities (modules)

boxpossession <- function(spat, digits=6) {
	proportions <- function(x) x / sum(x, na.rm=TRUE)
	possession <- spat
	for (k in seq(dim(possession)[3])) {
		for (i in seq(dim(possession)[1])) {
			possession[i,,k] <- proportions(possession[i,,k])
		}
	}
	possession[is.nan(possession)] <- NA
	possession <- round(possession, digits)
	return(possession)
}


# Function:
#	visualboxes: creates pie plots indicating box possession by different modules across time layers
# Arguments:
#	spat: an array of dimension boxes x modules x layers (produced by 'allocation'), which contains information about possession of boxes by different modules across a range of time layers 
#	palette: sequence of colors indicating module membership (20 colors by default)
#	lwd: width of the circle outline
#	border: color of the circle outline
#	res: resolution for drawing of circle outline
#	gap: gap between circles as a proportion of circle diameter
#	mai: outer margin size in inches
#	device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
#	file: name of pdf file if device == "pdf"
#	width: width of the device in inches

visualboxes <- function(spat, palette, lwd=1, border="black", res=200, gap=0.05, mai=0, device, file, width=7) {
	proportions <- function(x) x / sum(x, na.rm=TRUE)
	if (missing(palette)) {
		palette <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080")
	}
	palette <- palette[seq(dim(spat)[2])]
	gap <- rep_len(gap, 2)
	nr <- dim(spat)[1]
	nc <- dim(spat)[3]
	ng <- dim(spat)[2]
	 
	indices <- expand.grid(1:nr, 1:nc)
	yy <- nr:1 + gap[2] * seq(1, 2 * nr - 1, by=2)[nr:1] - 0.5
	xx <- 1:nc + gap[1] * seq(1, 2 * nc - 1, by=2)[1:nc] - 0.5
	xlim <- c(0, max(xx) + 0.5 + gap[1])
	ylim <- c(0, max(yy) + 0.5 + gap[2])
	height <- width * ylim[2] / xlim[2]
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "modules.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}

	par(mai=rep_len(mai, 4))
	plot(0, 0, type="n", xaxs="i", yaxs="i", axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
	for (i in seq(nrow(indices))) {
		prop <- proportions(spat[indices[i,1],,indices[i,2]])
		prop[is.nan(prop)] <- NA
		if (all(is.na(prop))) {
			col <- "white"
			prop <- 1
		} else {
			col <- palette[seq_along(prop)][!is.na(prop)]
			prop <- na.omit(prop)
		}
		partCircle(prop, col=col, x=xx[indices[i,2]], yy[indices[i,1]], size=1, lwd=lwd, border=border, res=res)
	}
	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
}


# Function:
#	partCircle: draws a pie plot
# Arguments:
#	prop: proportions
#	col: sequence of colors in order corresponding to prop
#	x, y: coordinates of the center
#	size: circle diameter
#	lwd: width of the circle outline
#	border: color of the circle outline
#	res: resolution for drawing of circle outline

partCircle <- function(prop, col, x=0, y=0, size=1, lwd=1, border=NA, res=200) {
	rs <- seq(0, 2 * pi, len=res)
	pts <- cbind(0.5 * cos(rs), 0.5 * sin(rs))
	pts <- size * pts + rep(1, res) %*% t(c(x,y))
	prop <- prop / sum(prop)
	col <- col[seq_along(prop)]
	parts <- round(quantile(seq(nrow(pts)), probs=c(0,cumsum(prop))))
	for (i in seq_along(prop)) {
		polygon(rbind(c(x,y), pts[parts[i]:parts[i+1],], c(x,y)), border="transparent", col=col[i]) 
	}
	if (!is.na(border)) {
		polygon(pts, border=border, col="transparent")
	}
}
