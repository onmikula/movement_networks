# Function:
#	socnet_plot: plots social network
# Arguments:
#	snet: social network specified as an adjacency matrix or an object of class 'igraph'
#	vertex.color: color of vertices (in the same order as in snet)
#	vertex.shape: shape of vertices (in the same order as in snet)
#	vertex.size: size of vertices
#	edge.width: multiplier of edge widths
#	vertex.label: whether to display vertex labels (default is FALSE)
#	vertex.label.cex: size of vertex labels
#	vertex.label.color: color of vertex labels
#	order: vector of vertex labels, cpecifying their order in circular layout of the graph
#	mai: size of outer margins in inches, recycled up to the length of four
#	device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
#	file: name of pdf file if device == "pdf"
#	width: width of the device in inches
#	height: height of the device in inches
#	... other arguments from igraph.plotting

socnet_plot <- function(snet, vertex.color="skyblue", vertex.shape="circle", vertex.size=10, edge.width=5, vertex.label=FALSE, vertex.label.cex=1, vertex.label.color="black", order=NULL, mai=0.02, device, file, width=7, height=7, ...) {
	mode <- ifelse(all(snet == t(snet)), "undirected", "directed")
	if (inherits(snet, "igraph")) {
		igr <- snet
	} else {
		igr <- igraph::graph_from_adjacency_matrix(snet, mode=mode, weighted=TRUE, diag=FALSE)	
	}
	edge.width <- edge.width * igraph::edge_attr(igr, "weight")
	if (is.null(order)) {
		ord <- seq(igraph::gorder(igr))
	} else {
		ord <- match(order, attr(igraph::V(igr), "name"))
	}
	coords <- igraph::layout_in_circle(igr, order=igraph::V(igr)[ord])
	rownames(coords) <- attr(igraph::V(igr), "names")
	if (isFALSE(vertex.label)) {
		vertex.label <- NA
	} else {
		vertex.label <- attr(igraph::V(igr), "names")
	}

	if (missing(device)) {
		if (.Platform$OS.type == "unix") {
			device <- "quartz"
		} else {
			device <- "x11"
		}
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "network.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}

	par(mai=rep_len(mai, 4))
	plot(igr, layout=coords, vertex.color=vertex.color, vertex.shape=vertex.shape, vertex.size=vertex.size, edge.width=edge.width, vertex.label=vertex.label, vertex.label.cex=vertex.label.cex, vertex.label.color=vertex.label.color, ...)

	if (isTRUE(device == "pdf")) {
		dev.off()
	}

}

