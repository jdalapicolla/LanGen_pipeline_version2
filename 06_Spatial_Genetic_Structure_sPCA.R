###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
##############   STEP 06: FINE-SCALE SPATIAL GENETIC STRUCTURE   ##############

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###


#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS TUTORIAL IS ORGANIZED IN DIFFERENT "Steps" AND INSIDE EACH STEP THERE ARE "Actions" DENOTED BY NUMBERS (#1) AND INSIDE EACH ACTION COULD EXIST "Operations" INDICATED BY LETTES (#A.)
#D. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


##2. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline


##3. GOALS FOR THIS STEP:
#A. PERFORM FINE-SCALE SPATIAL GENETIC STRUCTURE BASED ON sPCA MODEL


##4. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H


##5. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##6. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("functions_LanGen.R")


##7. INSTALL AND LOAD THE PACKAGES
#A. install the packages automatically
if("remotes" %in% rownames(installed.packages()) == FALSE){install.packages("remotes")
} else {print (paste0("'remotes' has already been installed in library"))}
if("BiocManager" %in% rownames(installed.packages()) == FALSE){install.packages("BiocManager")
} else {print (paste0("'BiocManager' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("devtools" %in% rownames(installed.packages()) == FALSE){install.packages("devtools")
} else {print (paste0("'devtools' has already been installed in library"))}

if("r2vcftools" %in% rownames(installed.packages()) == FALSE){remotes::install_github("nspope/r2vcftools")
} else {print (paste0("'r2vcftools' has already been installed in library"))}
if("LEA" %in% rownames(installed.packages()) == FALSE){BiocManager::install("LEA")
} else {print (paste0("'LEA' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("geosphere" %in% rownames(installed.packages()) == FALSE){install.packages("geosphere")
} else {print (paste0("'geosphere' has already been installed in library"))}
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("rgdal" %in% rownames(installed.packages()) == FALSE){install.packages("rgdal")
} else {print (paste0("'rgdal' has already been installed in library"))}
if("gridExtra" %in% rownames(installed.packages()) == FALSE){install.packages("gridExtra")
} else {print (paste0("'gridExtra' has already been installed in library"))}
if("dplyr" %in% rownames(installed.packages()) == FALSE){install.packages("dplyr")
} else {print (paste0("'dplyr' has already been installed in library"))}
if("vcfR" %in% rownames(installed.packages()) == FALSE){install.packages("vcfR")
} else {print (paste0("'vcfR' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("spdep" %in% rownames(installed.packages()) == FALSE){BiocManager::install("spdep")
} else {print (paste0("'spdep' has already been installed in library"))}
if("adespatial" %in% rownames(installed.packages()) == FALSE){BiocManager::install("adespatial")
} else {print (paste0("'adespatial' has already been installed in library"))}
if("usdm" %in% rownames(installed.packages()) == FALSE){install.packages("usdm")
} else {print (paste0("'usdm' has already been installed in library"))}
if("splancs" %in% rownames(installed.packages()) == FALSE){install.packages("splancs")
} else {print (paste0("'splancs' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(ggplot2, r2vcftools, geosphere, raster, rgdal, gridExtra, dplyr, LEA,  vcfR, dartR, adegenet, splancs, usdm, adespatial, spdep, usedist, ape)


##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Spatial_Structure/sPCA"))


#------------------------------------------------------------------------------
#                             1. Loading Files 
#------------------------------------------------------------------------------
## Fine-scale population strucutre using sPCA to compare the allelic frequency of individuals to allelic frequencies of their neighbours
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "pilocarpus"


###1.2. LOAD INPUTS
#A. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.

#B. Create a Genind object based on neutral VCFR
vcf_neutral = read.vcfR(paste0("vcf/", project_name,"_filtered_neutral_LEA_DAPC_TESS.vcf"), verbose = FALSE)


#------------------------------------------------------------------------------
#                         2. Subsetting Genotype File 
#------------------------------------------------------------------------------
###2.1. SUBSETTING USING A GENELIGHT OBJECT 
#A. Convert to GENELIGHT
gl = vcfR2genlight(vcf_neutral)

#B. Define pop as POP_ID from TESS
gl@pop = as.factor(snps_neutral@meta$PopID_tess)

#C. Create a subset of snps or individuals
#by individuals/clusters
index.ind = pop(gl) =="3" #selecting by pop
#check if the index worked
table(pop(gl), index.ind)
#by loci number. 500 random SNPs
set.seed(123)
index.loc = sample(nLoc(gl), 500, replace = F)
index.loc

#D. Subsetting
glsub = gl[index.ind, index.loc] # individual and loci
#glsub = gl[index.ind, ] # only individuals
#glsub = gl[,index.loc] # only snps
glsub = gl

#D. Convert to GENEIND
genind = gl2gi(glsub)
genind = gl2gi(gl)
genind
  

#------------------------------------------------------------------------------
#                         3. Geographical Information 
#------------------------------------------------------------------------------
#You can test different connection networks:
#Delaunay triangulation (type=1) = more suited to uniform sampling. Other softwares use this one, for comparasion it is a good one.
#
#Gabriel graph (type=2). The Gabriel graph is a subgraph of the Delaunay triangulation. It can be found in linear time if the Delaunay triangulation is given (Matula & Sokal 1980). The Gabriel graph contains, as subgraphs, the Euclidean minimum spanning tree (type=4), the relative neighborhood graph (type=3), and the nearest neighbor graph (type=6). 
#
#Relative neighbours (type=3)
#
#Minimum spanning tree (type=4)
#
#Neighbourhood by distance (type=5) =  You need to inform minimum (d0) and maximum distance (d1) IN NUMBER OF INDIVIDUALS 0 AND 2 IS A PATTERN. Biologically-defined distance would be better in this case (home-range, dispersion for plants)
#
#K nearests neighbours (type=6) = You need to inform the number of neighbours per point (argument k).Shoul be lower than one-third of the of the number of data points
#
#Inverse distances (type=7) = recommended when clustering of samples is present. This is not a true neighbouring graph: all sites are neighbours, but the spatial weights are directly proportional to the inversed spatial distances. You need to inform the minimum distance between any two distinct points (dmin) and the exponent of the inverse distance matrix (a).
#-----------------------------------------------------------------------------------

###3.1. GEOGRAPHIC DATA AND CONNECTION NETWORK
#A. Geographical information for all samples
coordinates = snps_neutral@meta[,4:5]

#B. Subset coordinates by populations:
for (i in 1:length(unique(snps_neutral@meta$PopID_tess))){
  pop = which(snps_neutral@meta$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}

coord_POP1 = coordinates[pop_TESS_1,]
coord_POP2 = coordinates[pop_TESS_2,]
coord_POP3 = coordinates[pop_TESS_3,]
coord_POP4 = coordinates[pop_TESS_4,]

#C. Create Connection network. Choose Neighbourhood by distance (type=5) for most cases.
#CN1 = chooseCN(coordinates, type=1, plot=T)
#CN2 = chooseCN(coordinates, type=2, plot=T)
#CN3 = chooseCN(coordinates, type=3, plot=T)
#CN4 = chooseCN(coordinates, type=4, plot=T)
CN = chooseCN(coordinates, type=5, d1 = 0, d2 = 2, plot=T)
#CN6 = chooseCN(coordinates, type=6, k = 5, plot=T)
#CN7 = chooseCN(coordinates, type=7, dmin = 3, a = 2, plot=T)

#D. Then convert your CN to a listw using the nb2listw function.
cn = nb2listw(CN)


#------------------------------------------------------------------------------
#                         4. Running sPCA Analysis 
#------------------------------------------------------------------------------
###4.1. RUNNING sPCA:
#Here we will run sPCA following the tutorial provide by Jombart.
#  http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf

#A. Number of PCs
n_pcs = nrow(genind@tab)
  
#B. Create a PCA object
input_scaled = scaleGen (genind, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = n_pcs)

#C. Running sPCA:
mySpca = multispati(pca_input, cn, scannf = F)
#verify number of PC in global and local structure
mySpca
plot(mySpca)
#verify lag vector onto the principal axesdata frame with n rows and (nfposi + nfnega) columns
head(mySpca$ls)


###4.2. TEST LOCAL AND GLOBAL STRUCTURE
# Lollipop plot interpretation: significant if lollipop is far from distribution of permutation values, not significant if lollipop falls along the distribution of permutation values.
# GLOBAL TEST: neighbouring individuals are more similar than expected
# LOCAL TEST: neighbouring individuals are more dissimilar than expected

#A. Global test
myGtest = global.rtest(pca_input$tab, cn, nperm=9999)
myGtest
plot(myGtest)

#B. Local test
myLtest = local.rtest(pca_input$tab, cn, nperm=9999)
myLtest
plot(myLtest)


#------------------------------------------------------------------------------
#                         5. Save the Results 
#------------------------------------------------------------------------------
###5.1. EXTRACT RESULTS TO DO THE INTERPOLATON ON QGIS
#verify lag vector onto the principal axesdata frame with n rows and (nfposi + nfnega) columns
head(mySpca$ls)
head(cbind(mySpca$ls,coordinates))
write.csv(cbind(mySpca$ls,coordinates), paste0("./Results_Spatial_Structure/sPCA/", project_name , "_sPCA_Axis_QGIS.csv"))

###5.2. SAVE RESULTS AS GRAPHS AND TABLES
#A. Tables
results_sPCA = summary(mySpca)
write.csv(results_sPCA, paste0("./Results_Spatial_Structure/sPCA/", project_name , "_sPCA_Summary.csv"))

#B. Spatial and Inertial Components:
pdf("./Results_Spatial_Structure/sPCA/Spatial_Inertia.pdf", onefile =F)
plot.new()
fgraph(mySpca)
dev.off()

#C. Barplot with Global and Local Structure in spectral color pattern
pdf("./Results_Spatial_Structure/sPCA/Global_Local_Structure_Spectral.pdf", onefile =F)
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="black")
dev.off()

#D. Barplot with Global and Local Structure in red color for PC significants
#Global Significant PCs:
n_posi = 2
#For Local PCs test diferent values to put color on graph, for no one PC let ir on 0. Usually you need to plus number of positive PC and Negative to put correctly colors:
n_eig = length(mySpca$eig)-0
#verify colors:
barplot(mySpca$eig, main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(n_posi,n_eig)))

pdf("./Results_Spatial_Structure/sPCA/Global_Local_Structure_Significant.pdf", onefile =F)
barplot(mySpca$eig, main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(n_posi,n_eig)))
dev.off()

#E. Colorplot of mySpca. It represents a cloud of points with colors corresponding to a combination of 1,2 or 3 quantitative variables, assigned to RGB (Red, Green, Blue) channels. For instance, this can be useful to represent up to 3 principal components in space.
#define the number of variables UP TO THREE (Positive + Negatives PCs)
variables = 2
pdf("./Results_Spatial_Structure/sPCA/Colorplot_Significant_PC.pdf", onefile =F)
colorplot(coordinates, as.data.frame(mySpca$ls[,1:variables]), cex=3, col = colors_spca, main="Colorplot of sPCA of Neutral SNPs")
dev.off()

#F. Global Test of mySpca
pdf("./Results_Spatial_Structure/sPCA/Global_Test_Significant.pdf", onefile =F)
plot(myGtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()

#G. Local Test of mySpca
pdf("./Results_Spatial_Structure/sPCA/Local_Test_Significant.pdf", onefile =F)
plot(myLtest, main="Simulation for \n Moran's Eigenvector Maps (MEMs)",  xlab = "Observed R²", ylab = "Frequency")
dev.off()


## Make Interpolated Maps using QGIS 3.4     
#1. Install "Processing" plug-in in QGis
#2. Go to "Add a delimited text file" and open the file with spca axes and coordinates
#3. Save the opened file as a shapefile
#4. Interpolate each axis separately using the Interpolation tool in the Processing Toolbox:
#  Processing Toolbox > Interpolation > Interpolation IDW.
#5. At the Interpolation window: 
# a) select the shapefile with spca axes; 
# b) add the axis in the second box as a vertorial entry layer;
# c) select the "extension";
# d) Adjust the number of columns and rows;
# e) run the interpolation (the raster will appear in the main window);
#6. Repeat the previous steps to each axis
#7. Combine the three interpolated axis raster in a rgb raster using Raster > Micellaneous > Mosaic:
# a) Select the three interpolated axis raster;
# b) Select the option "put each file in a separate band";
# c) Run to build your mosaic.

##END
