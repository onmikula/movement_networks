# movement_networks
R functions for analysis and visualisation of multilayer social networks from movement-tracking experiments.

In this type of experiments, visits of discrete locations (e.g. nest boxes) are recorded and the principal input is a list of entering and exit times of every individual at every location. Based on these data, the strength of interaction between two individuals in a given time period can be quantified as the total time spent together at the same location. The complete set of pairwise interactions can be formalised as a social network, which is a weighted undirected graph in this case with vertices representing individuals and edges representing their interactions. In the case of time-extended experiment, a sensible approach is to divide the whole time series into discrete time layers and calculate the social network for each of them. An ordered set of these layer-specific networks is the multilayer social network.

Every layer-specific network can contain groups of vertices that are more strongly connected to each other and specific groups can be present in multiple layers. If this is the case, it may be more parsimonious to describe the network as partitioned into discrete clusters of vertices, usually called modules. This is the conceptual basis for the minimum description length clustering formalised via map equation (Rosvall & Bergstrom 2008) and implemented in software Infomap (https://www.mapequation.org). Multilayer adaptation of map equation was introduced by Aslak et al. (2018). The script 'multilayer_social_networks.R' assumes Infomap binary to be placed in the folder 'bin'. Apart from Infomap and base R (R Core Team 2022), the functions depend also on packages igraph (Csardi & Nepusz 2006) and Matrix (Bates & Maechler 2019).

The script is organized into seven sections:
1. Loading raw data. These are data from the house mouse seminatural breeding experiment described in Mikula et al. (2022), namely from its 'domesticus 2013' run. The data are already divided into time layers.
2. Preparing the multilayer social network and writing its structure into the Infomap's input file.
3. Estimation of modules in Infomap.
4. Quantification of module statistics including compression rates as the measure of overall modularity and proportion of intramodular interactions calculated layer-by layer.
5. Visualisation of modules via individual x layer bar plot.
6. Analysis of spatial separation of modules, including calculation of spatial separation index and visualisation of nest box possession by different modules across layers.
7. Handling of layer-specific networks: subsetting, merging into summary networks and visualisation which may include interaction strength, module-membership and possibly additional information (e.g., sex of the individuals).

**When using these functions, please cite:**
- Mikula O, Macholán M, Ďureje Ľ, Hiadlovská Z, Daniszová K, Janotová K, 
Vošlajerová Bímová B (2022) House mouse subspecies do differ in their social structure.

**References:**
- Aslak U, Rosvall M, Lehmann S (2018) Constrained information flows in temporal networks reveal intermittent communities. Phys. Rev. E 97, 062312. (https://doi.org/10.1103/PhysRevE.97.062312)
- Bates D, Maechler M (2019) Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.2-18. (https://CRAN.R-project.org/package=Matrix)
- Csardi G, Nepusz T (2006) The igraph software package for complex network research. InterJournal, Complex Systems 1695, 1–9. (doi: http://igraph.org)
- R Core Team (2022) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
- Rosvall M, Bergstrom CT (2008) Maps of random walks on complex networks reveal community structure. Proc. Natl. Acad. Sci. USA 105, 1118–1123. (doi: https://doi.org/10.1073/pnas.0706851105)
