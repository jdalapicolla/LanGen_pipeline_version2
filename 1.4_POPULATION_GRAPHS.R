####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
###################### STEP 04: POPULATION CONNECTIVITY #######################




### Script prepared by Jeronymo Dalapicolla and Carolina S. Carvalho ###
### Based on DYER, Rodney J. Population graphs and landscape genetics. Annual Review of Ecology, Evolution, and Systematics, v. 46, p. 327-342, 2015. https://doi.org/10.1146/annurev-ecolsys-112414-054150




#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL
#A. THE FILE ".VCF" CLEANED AFTER FILTERING AND WITH DELIMITED GENETIC CLUSTERS, STEP 2.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline





##2. GOALS FOR THIS STEP
#A. ESTIMATE CONNECTIVITY AMONG POPULATIONS USING POPULATIONS GRAPHS.




##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H



##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT
rm(list=ls())





##5. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP
#You can download it in https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
source("functions_LanGen.R")





##6. INSTALL AND LOAD THE PACKAGES
#For r2vcftools do you need install VCFTools in you computer:https://vcftools.github.io/index.html
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
if (!require('r2vcftools'))   remotes::install_github("nspope/r2vcftools");          library('r2vcftools')
if (!require('popgraph'))     devtools::install_github("dyerlab/popgraph");          library('popgraph')
if (!require('gstudio'))      devtools::install_github("dyerlab/gstudio");           library('gstudio')
if (!require('radiator'))     devtools::install_github("thierrygosselin/radiator");  library('radiator')

#From CRAN R:
if (!require('strataG'))      install.packages("strataG");           library('strataG')
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('gghighlight'))  install.packages("gghighlight");       library('gghighlight')
if (!require('ggplot2'))      install.packages("ggplot2");           library('ggplot2')
if (!require('ggpubr'))       install.packages("ggpubr");            library('ggpubr')
if (!require('geosphere'))    install.packages("geosphere");         library('geosphere')
if (!require('sf'))           install.packages("sf");                library('sf')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('igraph'))       install.packages("igraph");            library('igraph')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2')
if (!require('graph4lg'))     install.packages("graph4lg");          library('graph4lg')
if (!require('fields'))       install.packages("fields");            library('fields')

#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, strataG, popgraph, gstudio, radiator, tidyverse, gghighlight, ggplot2, ggpubr, sf, raster, igraph, geosphere, reshape2, graph4lg, fields)




##7. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS
create_dir(c("./Results/Step04/Popgraph"))




##8. CREATE A PATTERN FOR GRAPHIC FOLDERS S TO SAVE THE RESULTS
theme_genetics = theme(axis.text=element_text(size=10, color="black"), #text in ticks axes
                       axis.title=element_text(size=12, face="bold"), #label axes
                       axis.line = element_line(colour = "black", size = 1, linetype = "solid"), #line on axes
                       axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"), #line on ticks
                       axis.ticks.length = unit(.25, "cm"), #ticks length
                       axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), #space between axis and label
                       axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), #space between axis and label
                       strip.text.x = element_text(size = 12, face="bold"), #facets label 
                       panel.grid.major = element_blank(), # remove grids
                       panel.grid.minor = element_blank(), # remove grids
                       panel.background = element_blank(), # remove background
                       panel.border = element_blank()) # remove borders)  





#### ANALYSIS ---- 


#### 1. LOAD FILES -----
#A. Load neutral .vcf file with geographical information and genetic clusters ID, chose in step 2.
#If you have a Genepop file, you can jump to Step #4
snps_neutral = vcfLink("vcf/pilocarpus_filtered_neutral.vcf", overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5277 SNPs.
names(snps_neutral@meta) #verify col names in metafile
head(snps_neutral@meta)





#### 2. ADD INFORMATION ON GEOGRAPHICAL AND POPULATIONS ----
#If you need put the information on locality/population and geographic coordinates in the metafile:
#A. Set Coordinates file path and load it:
coord_file = "Inputs/pilocarpus_coords.csv"
#load the file
coords = read.csv(coord_file)
#verify if it loaded correctly
head(coords)
tail(coords)

#B. Filter coords according to vcf individuals
meta = snps_neutral@meta%>%
 setNames(., c("Sample_Number", "Sample_ID", "Local", "Longitude", "Latitude", "Missing_PC", "ID"))

set_coords = coords %>%
  as.data.frame() %>%
  setNames(., c("Local", "ID", "Longitude", "Latitude", "Aggregate"))
set_coords = semi_join(set_coords, meta, by = "ID")
set_coords = inner_join(meta, set_coords, by = "ID")
head(set_coords)
set_coords = set_coords[,-c(3:5)] %>%
  setNames(., c("Sample_Number", "ID_original", "Missing_PC", "ID", "Local", "Longitude", "Latitude", "Aggregate"))
head(set_coords)


#C. Convert coordinates in UTM if you need. My original data is WGS and the UTM data will be in SAD69 UTM zone 22S. Search in internet the EPSG code for that area. For example https://epsg.io/ #In the case is EPSG:32722.
set_coords = st_as_sf(set_coords, #data frame
                      coords = c("Longitude", "Latitude"), # for point data
                      remove = F, # don't remove these lat/lon cols from df
                      crs = 4326) %>% # add original projection (this is WGS84)
  st_transform(crs =st_crs(32722)) %>% #choose a new EPSG code
  mutate(X_UTM = sf::st_coordinates(.)[,1], Y_UTM = sf::st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  setNames(., c("Sample_Number", "ID_original", "Missing_PC", "ID", "Local", "Longitude", "Latitude", "Aggregate", "Geometry", "X_UTM", "Y_UTM")) %>%
  dplyr::select(., -Geometry)
head(set_coords)


#D. Check the order:
identical(as.character(set_coords$ID_original), as.character(snps_neutral@meta$sample_name))
setdiff(as.character(set_coords$ID_original), as.character(snps_neutral@meta$sample_name))
setdiff(as.character(snps_neutral@meta$sample_name), as.character(set_coords$ID_original))


#E. Add the information in the metafile. You can add more columns if you need, like ecological traits or other clustering classification
head(snps_neutral@meta) #we have sample_num and sample_name
snps_neutral@meta = inner_join(snps_neutral@meta, set_coords, by = c("sample_name" = "ID_original"))
head(snps_neutral@meta)
snps_neutral@meta = snps_neutral@meta[, -c(3:7)]
rownames(snps_neutral@meta) = snps_neutral@meta$sample_name
head(snps_neutral@meta)


#F. Save vcf with correct metafile:
Save(snps_neutral, "vcf/pilocarpus_filtered_neutral.vcf")
#reload the file if you need in Step #1 








#### 3. CONVERTING VCF TO GENEPOP -----
#A. Load using vcfR
snpsR = read.vcfR("vcf/pilocarpus_filtered_neutral.vcf", verbose = T)

#B. Convert vcf into Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#C. Convert Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

#D. Save Gtypes into Genepop. This function do not accept a path to save. You need to move the file to other folder. It will be save in your working directory.
genepopWrite(snps_gtypes, "pilocarpus_filtered_neutral_genepop")







#### 4. LOAD GENEPOP FILES AND EDIT INPUT FOR GRAPHS ----
#A. Load Genepop file from the folder you moved it. Check the names: 
snps_genepop = read_population("./Results/Step04/Popgraph/pilocarpus_filtered_neutral_genepop.txt", type="genepop")
snps_genepop[1:10,1:10]
str(snps_genepop)

#B. Names in Genepop must be identical to metafile or .csv. So, we need to remove the pattern "Default "
snps_genepop2 = snps_genepop %>%
  as.data.frame() %>%
  mutate (sample_name = str_remove(ID, pattern="Default "))
snps_genepop2$sample_name

#C. Select in the metafiles or in a other .csv the information on localities/populations that you wanna test the connectivity
names(snps_neutral@meta) # check columns 
locations = snps_neutral@meta[,c(2,5:9)] # select cols that you want to put in genepop file
head(locations)
#Check if you have populations/locations with less than 3 individuals. We need to remove them! 
table(locations$Local)
table(locations$Aggregate)
#run this line below if you have <3 individuals
locations = locations[locations$Local != names(which(table(locations$Local) < 3)), ]
  

#D. Merge the files
input_popgraph = merge(locations, snps_genepop2, by="sample_name")
#check
input_popgraph[1:10,1:10]
str(input_popgraph)






#### 5. RUN POPULATION GRAPHS USING POPGRAPH ----
#A. Run without spatial information
pop_graph = popgraph(to_mv(input_popgraph), groups=input_popgraph$Local, alpha = 0.05) 
#check
plot(pop_graph)
class(pop_graph)
pop_graph

#Information on pop_graph:
#Node = Vertex = V = strata/populations/localities
#Link = Edge = E = bidirectional because covariance is a symmetric trait 


#B. To plot with spatial information you need coordinates for population/locations. I choose one individual to represent the local but you can calculate a centroid.  
input_popgraph_spatial = as.data.frame(matrix(NA, 21,4)) #21 is the number of my locations
input_popgraph_spatial
loop = 0

#select one individual by location/population
for (i in unique(input_popgraph$Local)) {
  loop = loop+1
  positions = which(input_popgraph$Local == i)
  input_popgraph_spatial[loop,] = input_popgraph[positions[1], c(1,3,4,5)]
}

#rename cols
input_popgraph_spatial = input_popgraph_spatial %>%
  setNames(., c("ID_original", "Local", "Longitude", "Latitude"))
input_popgraph_spatial


#C. Run Population Graphs with spatial information. By location.
pop_graph_spatial = decorate_graph(pop_graph, input_popgraph_spatial, stratum="Local")
plot(pop_graph_spatial)


#D. Plot in the geographic space
ggplot() + 
  geom_edgeset( aes(x=Longitude,y=Latitude), pop_graph_spatial) +
  geom_nodeset( aes(x=Longitude, y=Latitude), pop_graph_spatial, size=4) +
  theme_genetics


#E. Plot in the geographic space with a raster
#define a extension that cover all points for crop. If your study are is wide R can collapse. Better use a QGIS after to visualization
plot(snps_neutral@meta$Longitude, snps_neutral@meta$Latitude)
ext = extent(-50.7, -49.7, -6.5, -5.8)

#Load and crop the raster
landcover = raster::stack("maps/rasters/LC_Carajas_2018.tif")
landcover = crop(landcover, ext)

#Convert raster to dataframe to use in ggplot 
spdf = as(landcover, "SpatialPixelsDataFrame")
df = as.data.frame(spdf)
colnames(df) = c("value", "x", "y")

#plot in ggplot
ggplot() +
  geom_tile(data=df,aes(x=x, y=y, fill=value), alpha=0.8) +
  geom_edgeset( aes(x=Longitude,y=Latitude), pop_graph_spatial, color="white") +
  geom_nodeset( aes(x=Longitude, y=Latitude, size=size), pop_graph_spatial) +
  theme_genetics



#### 6. ESTIMATE METRICS FOR NODES/VERTEX IN POPULATION GRAPHS USING POPGRAPH ----
#A. Create a data frame with all vertex with Lat and long 
df.nodes = data.frame(Pop=V(pop_graph_spatial)$name, Latitude=V(pop_graph_spatial)$Latitude, Longitude=V(pop_graph_spatial)$Longitude)

#B. Closeness:
# a weighted measure of distance from the node to all other nodes.
# measures how many steps is required to access every other vertex from a given vertex.
# evaluates the extent to which a local/Pop is genetically similar to all other locales in the network
# higher closeness values indicate higher genetic distance and less connectivity
df.nodes$closeness = closeness(pop_graph_spatial, mode = c("total")) #for undirectional is a total metric

#C. Betweenness:
# a ranking of how many of the shortest paths through the graph go through a specific node
# sum of the shortest connections
# higher betweenness values indicate a local/Pop has relatively more paths that flow through it from other locales; more connectivity
df.nodes$betweenness = betweenness(pop_graph_spatial)

#D. Degree:
#absolute number of genetic connections a locale has with other locales in the network
#lower values may indicate a locale is more isolated
#higher values of ‘degree’ as indicative of greater gene flow
df.nodes$degree = degree(pop_graph_spatial)

#E. "Eigenvector" Centrality:
#vertices with high eigenvector centralities are those which are connected to many other vertices which are, in turn, connected to many others (and so on).
#higher values as indicative of greater gene flow
#reflects the extent to which each locale is centrally located as a hub of genetic connectivity
df.nodes$eigenCent = evcent(pop_graph_spatial)$vector

#F. "Hub" Centrality
#undirected matrices the adjacency matrix is symmetric and the hub scores are the same as authority scores and eigenvector centralities 
df.nodes$hub = hub.score(pop_graph_spatial)$vector

#G. "Authority" Centrality
#undirected matrices the adjacency matrix is symmetric and the hub scores are the same as authority scores and eigenvector centralities 
df.nodes$authorities = authority.score(pop_graph_spatial)$vector

#H. Strength:
#Summing up the edge weights of the adjacent edges for each vertex.
#Its like a 'degree' metrics that account the weights of links between nodes. Degree is a unweighted version of strength;   
df.nodes$strength = strength(pop_graph_spatial)

#I. Size:
#node size: measure of the within population genetic variance
df.nodes$size = V(pop_graph_spatial)$size

#J. Modularity:
#Module = clustering of nodes with the graph
#Modules:
plot(cluster_walktrap(pop_graph_spatial),pop_graph_spatial)
#Modularity Score:
modularity(pop_graph_spatial, membership(cluster_walktrap(pop_graph_spatial)))
#Dendograms
eb = cluster_edge_betweenness(pop_graph_spatial)
plot(as.dendrogram(eb))
df.nodes$modules = eb$membership

#K. Sample size by locations/Population
df.nodes$Freq = as.data.frame(table(input_popgraph$Local))$Freq

#L. Check and save the results:
head(df.nodes)
str(df.nodes)
write.csv(df.nodes, "./Results/Step04/Popgraph/Metrics_Nodes_Pilocarpus.csv")


#M. Boxplot for all metrics
#select the columns to create a df
df = melt(df.nodes[, -c(2,3,8,9,12,13)], id.vars=c("Pop"))
#save it
pdf("./Results/Step04/Popgraph/Boxplot_NodesMetrics.pdf", onefile =F)
ggplot(df)+
  geom_boxplot(aes(x=variable, y= value), position = position_dodge2(preserve = "single"))+
  scale_x_discrete(labels= c("closeness" = "Closeness", "betweenness" = "Betweenness", "degree" = "Degree", "eigenCent" = "Centrality", "size" = "Size", "strength" = "Strength")) +
  facet_wrap(variable ~., scales="free", ncol =3)+
  theme_genetics+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y= NULL, x = NULL)
dev.off()

#N. Check the outliers in boxplot:
is_outlier = function(x) {  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))}
df.nodes$Pop[which(is_outlier(df.nodes$betweenness))]
df.nodes$Pop[which(is_outlier(df.nodes$closeness))]
df.nodes$Pop[which(is_outlier(df.nodes$degree))]
df.nodes$Pop[which(is_outlier(df.nodes$eigenCent))]
df.nodes$Pop[which(is_outlier(df.nodes$strength))]
df.nodes$Pop[which(is_outlier(df.nodes$size))]





#### 7. ESTIMATE METRICS FOR LINKS/EDGES IN POPULATION GRAPHS USING POPGRAPH ----
#A. Create a data frame with all edges anvertex with Lat and long 
df.edges = data.frame(Edge=as_ids(E(pop_graph_spatial)))

#B. Weigth:
#connection strengths
df.edges$weight = E(pop_graph_spatial)$weight

#C. Edge's Betweenness
df.edges$betweenness = edge.betweenness(pop_graph_spatial)

#D. Check and Save Results
head(df.edges)
write.csv(df.edges, "./Results/Step04/Popgraph/Metrics_Edges_Pilocarpus.csv")

#E. Boxplot for all metrics
#select the columns to create a df
df = melt(df.edges, id.vars=c("Edge"))
#save it
pdf("./Results/Step04/Popgraph/Boxplot_LinkMetrics.pdf", onefile =F)
ggplot(df)+
  geom_boxplot(aes(x=variable, y= value), position = position_dodge2(preserve = "single"))+
  scale_x_discrete(labels= c("betweenness" = "Betweenness", "weight" = "Weigth")) +
  facet_wrap(variable ~., scales="free", ncol =3)+
  theme_genetics+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y= NULL, x = NULL)
dev.off()

#N. Check the outliers in boxplot:
is_outlier = function(x) {  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))}
df.edges$Edge[which(is_outlier(df.edges$betweenness))]
df.edges$Edge[which(is_outlier(df.edges$weight))]





#### 8. EXPORT GRAPHS AND ALL METRICS AS SHAPEFILE ----
#A. Save popgraph as a new object
pop_graph_spatial_full = pop_graph_spatial

#B. Add nodes/vertex metrics. Size is already there!. Need to once at time.
names(df.nodes) 
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Closeness", index = V(pop_graph_spatial_full), df.nodes$closeness)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Betweenness", index = V(pop_graph_spatial_full), df.nodes$betweenness)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Degree", index = V(pop_graph_spatial_full), df.nodes$degree)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Eig_Centrality", index = V(pop_graph_spatial_full), df.nodes$eigenCent)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Hub_Centrality", index = V(pop_graph_spatial_full), df.nodes$hub)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Aut_Centrality", index = V(pop_graph_spatial_full), df.nodes$authorities)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Strength", index = V(pop_graph_spatial_full), df.nodes$strength)
pop_graph_spatial_full = set_vertex_attr(pop_graph_spatial_full, "Modules", index = V(pop_graph_spatial_full), df.nodes$modules)

#C. Add link/edges metrics. Weitgth is already there!. Need to once at time.
pop_graph_spatial_full = set_edge_attr(pop_graph_spatial_full, "Betweenness_edge", index = E(pop_graph_spatial_full), df.edges$betweenness)

#D. Check
pop_graph_spatial_full

#E. Create a a file with geographic information to be the shapefile
crds_patches  = input_popgraph_spatial[,-1]
colnames(crds_patches) = c("ID", "x", "y")
crds_patches

#F. Export
graph_to_shp(graph = pop_graph_spatial_full, crds = crds_patches, mode = "both",
             layer = 'pop_graph_spatial_full_pilocarpus', 
             dir_path = "./maps/shapefiles/",
             metrics = TRUE,
             crds_crs = 4326
)




#### 9. CONDITIONAL GENETIC DISTANCE (cGD) ----
#length of the shortest path through the topology connecting populations and is estimated directly from the adjacency matrix
#length of the edges that connect nodes; The parameter cGD has been used in models of IBD.
#normal edges  = edge length is proportional to that expected under a model of isolation by distance
#extended edges = edges whose length is longer than expected; which are consistent with long-distance dispersal events or range expansion;
#compressed edges = edges whose length is shorter than expected which indicate a reduced permeability of intervening landscape features.

#A. Estimate geographic distance in Km if you used UTM in popgraph
pDist = fields::rdist.earth(cbind( V(pop_graph_spatial)$Longitude, V(pop_graph_spatial)$Latitude), miles=F)
write.csv(pDist, "./Results/Step04/Popgraph/Geographic_Distance_UTM_Pilocarpus.csv")

#B. Estimate geographic distance in Km if you used Lat/Long spheric in popgraph
pilo_coord = as.data.frame(cbind(V(pop_graph_spatial)$name, V(pop_graph_spatial)$Longitude, V(pop_graph_spatial)$Latitude))
pilo_coord[, c(2,3)] = sapply(pilo_coord[, c(2,3)], as.numeric)
str(pilo_coord)
coordinates(pilo_coord) = pilo_coord[,c(2,3)]
projection(pilo_coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")
pDist = geosphere::distm(pilo_coord, fun=distGeo)/1000
write.csv(pDist, "./Results/Step04/Popgraph/Geographic_Distance_WGS_Pilocarpus.csv")

#C. Convert popgraphs in distance matrices to use in Mantel Tests or MLPE. You can use three parameters:
#shortest path among pops/locations:
cGD1 = to_matrix(pop_graph_spatial, mode="shortest path")
cGD1 
write.csv(cGD1, "./Results/Step04/Popgraph/Genetic_Distance_ShortestPath_Pilocarpus.csv")

#binary matrix representing the pairs of connected nodes
cGD2 = to_matrix(pop_graph_spatial, mode="adjacency")
cGD2
write.csv(cGD2, "./Results/Step04/Popgraph/Genetic_Distance_Adjacency_Pilocarpus.csv")

#Similar to the adjacency matrix but using edge weights instead of binary values
cGD3 = to_matrix(pop_graph_spatial, mode="edge weight")
cGD3
write.csv(cGD3, "./Results/Step04/Popgraph/Genetic_Distance_Edge_Weight_Pilocarpus.csv")




#10. MAXIMUM, MINIMUM AND AVARAGE GEOGRAPHIC DISTANCE BETWEEN CONNECTED POPULATIONS/LOCATIONS 
#A. Rename cols and rows of geographic distance
colnames(pDist) = colnames(cGD2) 
rownames(pDist) = rownames(cGD2) 

#B. Convert geographic matrix into data.frame and create a col named link using the two nodes
df_geoDist = as.data.frame(as.table(pDist))
df_geoDist$link = paste(df_geoDist$Var1, df_geoDist$Var2) 
head(df_geoDist)

#C. Convert genetic matrix into data.frame. Use the binary matrix and create a col named link using the two nodes
df_cGD2 = as.data.frame(as.table(cGD2))
df_cGD2$link = paste(df_cGD2$Var1, df_cGD2$Var2) 
head(df_cGD2)

#D. Merge both data.frame by col link
df_full = merge(df_cGD2[,3:4], df_geoDist[,3:4], by="link")
head(df_full)

#E. Rename the columns:
colnames(df_full) = c("Link","Connection","Geo_Dist")
df_full_connected = subset(df_full, Connection > 0)

#F. Check the Histogram of geo distance
hist(df_full_connected$Geo_Dist, xlab="Distance between Connected Locations", main="")

#G. Distances between connected locations:
min(df_full_connected$Geo_Dist)
max(df_full_connected$Geo_Dist) # 68.87 Km
mean(df_full_connected$Geo_Dist)# 18.91 Km

#END;
