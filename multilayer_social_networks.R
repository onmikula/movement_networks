### TOOLS
library(igraph)
source("func/make_social_networks.R")
source("func/modify_social_networks.R")
source("func/plot_social_networks.R")
source("func/infomaptools.R")
source("func/multilayer_modules.R")
source("func/spatial_structure.R")
findmode <- function(x) as.numeric(names(which.max(table(x))))



# LOAD RADIOTRACKING RECORDS
# loads object 'rfid' - a list of radio-frequency identification records from subsequent time layers
# individual movements within any given time layer are recorded in data frame with four columns:
# ID = individual name, Time = time of box entering/exiting, Box = nest box identifier, Enter = whether the individual enters (TRUE) or exits (FALSE) the box

load("data/rfid.RData")



# SOCIAL NETWORKS
# creates a multilayer network, which is an ordered series of single layer networks that may differ in size and composition
# the networks are represented by their adjacency matrices (whose entries are proportions of time spent together in any time layer)
networks <- make_multilayer_network(rfid, normalize=TRUE, diagonal=FALSE)
save(networks, file="data/networks.RData")

# exports the networks into text file (in so-called pajek format), which is used as an input for multilayer network Infomap
writeMultilayer(networks, "data/pairwise_multilayer_network.net")

# records starting and ending times of the networks
networks_times <- lapply(lapply(rfid, "[[", "Time"), range)
save(networks_times, file="data/networks_times.RData")



### MODULE ESTIMATION (INFOMAP)
# runs Infomap (whose binary is supposed to be found in 'bin' folder of the working directory)
# the relax rate is here set to 0.60
system2(command="./bin/Infomap", args=c("-i", "multiplex", "--multiplex-js-relax-rate", "0.60", "--two-level", "--expanded", "--clu", "data/pairwise_multilayer_network.net", "data"))

# reads in output of Infomap and stores it in the form of matrix with rows corresponding to individuals, columns to time layers and entries indicating module membership (NA if an individual left no record in the particular time layer) 
modules <- timeclusters("data/pairwise_multilayer_network_expanded.clu", net="data/pairwise_multilayer_network.net", relabel=TRUE, ord=TRUE)
write.table(modules, "data/modules.txt", quote=FALSE, sep="\t")



### MODULE STATISTICS
# extracts compression rate as a measure of modularity of network partitioned into map equation (=Infomap) modules
M <- extractCompressionRate("data/pairwise_multilayer_network_expanded.clu")

# calculates proportion of intramodular interactions, i.e., distinctiveness of the modules
P <- propinside(networks, modules)
plot(P, ylim=c(0,1), ylab="Proportion", xlab="Time layer", main="Intramodular interactions", cex.lab=1.25, cex.main=1.25)



### MODULE VISUALIZATION
# reads in the module membership information 
modules <- read.table("data/modules.txt")

# basic plotting
visualmodules(modules)

# export to the file "modules.pdf"
visualmodules(modules, device="pdf", file="plots/modules.pdf")

# annotation of x axis with the ranking of layers and y axis with individual names
visualmodules(modules, xaxis=TRUE, yaxis=TRUE, cex.axis=c(1,0.5), xlab="layers", cex.lab=1.25)

# annotation of x axis with times in units of days and y axis with individual names
days <- as.integer(sapply(lapply(lapply(rfid, "[[", "Time"), range), mean) / (24 * 3600))
visualmodules(modules, xaxis=TRUE, yaxis=TRUE, cex.axis=c(1,0.5), times=days, xlab="days", cex.lab=1.25)



### SPATIAL SEPARATION & BOX POSSESSION
# 'allocation' quantifies time spent by community (module) members in different boxes,
# 'boxpossession' converts the times into proportions
spatial <- allocation(rfid, modules)
possession <- boxpossession(spatial, digits=3)

# quantifies spatiol separation of modules across different time layers
separation <- spatsepar(spatial)

# visualization of 'allocation' (or 'boxpossession') output in the form of a grid of pie plots
which_to_show <- 91:100
visualboxes(spatial[,,which_to_show])
visualboxes(spatial, gap=c(0.1, 1), device="pdf", file="plots/boxpossession.pdf", width=12)



### SUMMARY NETWORKS
# loads social networks characterized by a list of adjacency matrices (one per time-layer)
load("data/networks.RData")

# loads starting and ending times of the networks
# and specifies starting date of the experiment from which the times (in sec) are counted
load("data/networks_times.RData")
start <- as.Date("2013-05-13")

# reads in information about individuals
indiv <- read.delim("data/individuals.txt")

# reads in the module membership information 
modules <- read.table("data/modules.txt")

# a subset of palette created by A. Trubetskoy (https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors)
palette <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe")

# specifies which networks to summarize
which_to_sum <- 96:100

# creates the summary network as a union of specified networks
merged <- socnet_merge(snets=networks, which=which_to_sum, what="union", times=networks_times)

# subsets the summary network (retaining only individuals older than 50 days at the end of the summarized period)
end <- start + max(networks_times[[rev(which_to_sum)[1]]]) / (24 * 60 * 60)
subset <- indiv$ID[(end - as.Date(indiv$Birth)) >= 50]
merged_subs <- socnet_subset(snets=merged, subsets=subset)

# graphical parameters based on the prevailing module membership and sex of the individuals
members <- apply(modules[rownames(merged_subs), which_to_sum], 1, findmode)
color <- palette[members]
shape <- c(F="circle", M="square")[indiv$Sex[match(rownames(merged_subs), indiv$ID)]]
constraint <- names(sort(members))

# plots the summary network
# basic plotting
socnet_plot(merged_subs)
# advanced plotting
socnet_plot(snet=merged_subs, vertex.color=color, vertex.shape=shape, vertex.size=10, edge.width=10, mai=0.02, vertex.label=TRUE, vertex.label.cex=0.4, vertex.label.font=2, order=constraint, device="pdf", file="plots/summary_net_modules.pdf")

