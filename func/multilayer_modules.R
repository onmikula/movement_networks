# Function:
#	timeclusters: extracts modules from the Infomap output
# Arguments:
#	clu: .clu output of Infomap
#	net: .net input of Infomap
#	relabel: logical, whether to label clusters (=modules) according to their appearance in the table
#	ord: logical, whether to sort rows in the order of (1) module no., (2) earliest layer of occurence and (3) the timespan of the record 
# Value:
#	matrix with individuals in rows and time layers in columns, filled with numbers indicating module membership (NA = no record of an individual in that time layer) 

timeclusters <- function(clu, net=NULL, relabel=FALSE, ord=FALSE) {
	x <- readLines(clu)
	h <- which(substr(x, 1, 1) == "#")
	header <- unlist(strsplit(gsub(":", "", gsub("# ", "", rev(x[h])[1])), " "))
	nodes <- do.call(rbind, lapply(sapply(x[-h], strsplit, " "), as.numeric))
	dimnames(nodes) <- list(NULL, header)
	nodes <- nodes[nodes[,"flow"] > 0,]
	mat <- matrix(0, max(nodes[,"node"]), max(nodes[,"layer"]))
	mat[nodes[,"node"] + max(nodes[,"node"]) * (nodes[,"layer"] - 1)] <- nodes[,"cluster"]
	if (!is.null(net)) {
		x <- readLines(net)
		vertices <- x[(grep("^\\*Vertices", x) + 1):(grep("^\\*Intra|^\\*multiplex", x) - 1)]
		vertices <- gsub("\\\"", "", sapply(strsplit(vertices, " "), "[", 2))
		rownames(mat) <- vertices	
		layers <- grep("#\\s*layers\\:", x)
		if (length(layers) > 0) {
			layers <- as.numeric(unlist(strsplit(gsub("#\\s*layers\\:", "", x[layers]), ",\\s*")))
		} else {
			layers <- seq(ncol(mat))
		}
		colnames(mat) <- paste0("L", formatC(layers, format="d", flag=0, width=nchar(max(layers))))
	}
	mat <- mat[rowSums(mat > 0) > 0,]
	mat[mat == 0] <- NA
	if (isTRUE(relabel)) {
		mat[,] <- match(as.numeric(mat), unique(na.omit(as.numeric(mat))))
	}
	if (isTRUE(ord)) {
		group <- setNames(lapply(seq(nrow(mat)), function(i) table(na.omit(mat[i,]))), rownames(mat))
		group <- as.numeric(sapply(group, function(x) names(x)[which.max(x)]))
		start <- apply(!is.na(mat), 1, which.min)
		length <- apply(!is.na(mat), 1, sum)
		ord <- order(group, ncol(mat) - start, ncol(mat) - length, decreasing=FALSE)
		mat <- mat[ord,]
	}
	attr(mat, "layers") <- layers
	return(mat)
}



# Function:
#	visualmodules: visualisation of multilayer modules
# Arguments:
#	mod: matrix or data.frame with entries indicating module membership (the output of 'timeclusters')
#	palette: sequence of colors indicating module membership (20 colors by default)
#	device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
#	file: name of pdf file if device == "pdf"
#	width: width of the device in inches
#	asp: aspect ratio of a cell
#	mai: margin size specified in inches (recycled up to the length of four and adujsted if xaxis or yaxis is TRUE)
#	times: endpoints of layers in units of time
#	xaxis: whether to display x axis (ranks or times of layers)
#	yaxis: whether to display y axis (row names)
#	cex.axis: size of tick labels on x and y axis, when two numbers are given, the second applies to the y axis
#	xlab, ylab: axis labels
#	cex.lab: size of x and y axis labels, when two numbers are given, the second applies to the y axis
#	pos.lab: position of x and y axis labels (~ argument line of mtext), when two numbers are given, the second applies to the y axis
#	nticks: no. of tickmarks on x axis
#	expr: expression adding extra features to the plot (e.g., annotation of rows and columns)

visualmodules <- function(mod, palette, device, file, width=7, asp=1, mai=0.02, times=NULL, xaxis=FALSE, yaxis=FALSE, cex.axis=1, xlab="", ylab="", cex.lab=1, pos.lab=1, nticks=5, expr) {

	findticks <- function(x, k) round(x[round(seq(1,length(x),length=k+2))])

	mod <- as.matrix(mod)
	n <- nrow(mod)
	m <- ncol(mod)
	k <- max(mod, na.rm=TRUE)
	layers <- attr(mod, which="layers")
	if (is.null(layers)) {
		layers <- seq(m)
	}
	if (max(layers) > m) {
		misslayers <- setdiff(seq(max(layers)), layers)
		mod <- cbind(mod, matrix(0, n, length(misslayers)))
		mod <- mod[,order(c(layers, misslayers))]
		m <- ncol(mod)
	}
	if (is.null(times)) {
		times <- seq(m)
	}
	if (missing(palette)) {
		palette <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080")
	}
	palette <- palette[seq(k)]

	cex.axis <- rep_len(cex.axis, 2)
	cex.lab <- rep_len(cex.lab, 2)
	pos.lab <- rep_len(pos.lab, 2)
	pos.axis <- 0.5
	oneline <- 0.2

	strw <- 0.09269206
	strh <- 0.1388889
	mai <- rep_len(mai, 4)
	if (isTRUE(xaxis)) {
		mai[1] <- max(c(mai[1], 0.02 * min(c(width, height)) + strh * cex.axis[1]))
		if (nchar(xlab) > 0) {
			mai[1] <- mai[1] + 1.5 * strh * cex.lab[1] + pos.axis * oneline
			pos.lab[1] <- mai[1] / oneline
		}
		mai[1] <- mai[1] + 2 * strh * cex.axis[1]
		mai[2] <- mai[2] + min(nchar(times)) * strw * cex.axis[1] / 2
		mai[4] <- mai[4] + max(nchar(times)) * strw * cex.axis[1] / 2
	}
	if (isTRUE(yaxis)) {
		mai[2] <- max(c(mai[2], 0.02 * min(c(width, height)) + strw * cex.axis[2]))
		if (nchar(ylab) > 0) {
			mai[2] <- mai[2] + 1.5 * strw * cex.lab[2] + pos.axis * oneline
			pos.lab[2] <- mai[2] / oneline
		}
		mai[2] <- mai[2] + (max(nchar(rownames(mod))) + 1) * strw * cex.axis[2]
		mai[1] <- mai[1] + strh * cex.axis[2] / 2
		mai[3] <- mai[3] + strh * cex.axis[2] / 2
	}

	xlim <- c(0, m)
	height <- asp * width * n / xlim[2]
	width <- width + sum(mai[c(2,4)])
	height <- height + sum(mai[c(1,3)])
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

	par(mai=mai)
	plot(0, 0, xlim=xlim, ylim=c(0,n) + 0.5, type="n", asp=1, xaxs="i", yaxs="i", axes=FALSE, ann=FALSE, bty="n")
	graphics::image(x=times, y=seq(n), z=t(mod[n:1,]), zlim=c(0, k), col=palette, add=TRUE)
	if (isTRUE(xaxis)) {
		at <- findticks(seq(m), k=nticks)
		axis(1, at=at - 0.5, labels=times[at], cex.axis=cex.axis[1], line=pos.axis)
	}
	if (isTRUE(yaxis)) {
		mtext(text=rev(rownames(mod)), side=2, at=seq(n), line=pos.axis, cex=cex.axis[2], las=2)
	}
	if (nchar(xlab) > 0) {
		mtext(side=1, text=xlab, line=pos.lab[1], cex=cex.lab[1])
	}
	if (nchar(ylab) > 0) {
		mtext(side=2, text=ylab, line=pos.lab[2], cex=cex.lab[2])
	}
	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		} else {
			warning("'expr' argument cannot be evaluated as it is not of class 'expression'")
		}
	}
	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
	
}


# Function:
#	propinside: quatifies a proportion of interaction strength allocated inside the modules
# Arguments:
#	A: a matrix or a list of them indicating links in the layer-specific network(s)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
#	G: a vector (possibly factor) or a list of them indicating module membership of vertices from the network(s)
# Value:
#	proportion of weights of links that join vertices belonging to the same modules

propinside <- function(A, G) {
	intramodular <- function(A, G) {
		diag(A) <- 0
		g <- G[match(rownames(A), names(G))]
		g <- match(g, unique(g))
		g <- diag(max(g))[g,]
		d <- g %*% t(g)
		W <- sum(A)
		return(sum(A * d) / W)
	}
	if (inherits(A, "matrix")) {
		A <- list(A)
	}
	if (inherits(G, "matrix") | is.data.frame(G)) {
		G <- lapply(split(t(G), seq(ncol(G))), setNames, nm=rownames(G))
	} else if (!is.list(G)) {
		G <- list(G)
	}
	G <- lapply(G, na.omit)
	return(unlist(Map(intramodular, A, G)))
}
