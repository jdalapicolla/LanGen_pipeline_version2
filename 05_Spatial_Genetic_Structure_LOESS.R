###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
##############   STEP 05: FINE-SCALE SPATIAL GENETIC STRUCTURE   ##############

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
#A. PERFORM FINE-SCALE SPATIAL GENETIC STRUCTURE BASED ON LOESS MODEL.


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

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(ggplot2, r2vcftools, geosphere, raster, rgdal, gridExtra, dplyr, LEA)


##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Spatial_Structure/LOESS"))


#------------------------------------------------------------------------------
#                        1. Loading Files
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "pilocarpus"


###2.1. LOAD VCF FILES AND DISTANCES: 
#A. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.


#B. Load distance files among all individuals calculate in STEP 3:
#relatedness
REL = read.csv(paste0("Results_Distance/Yang_Reletedness_IND_neutral_", project_name, ".csv"), row.names = 1)
head(REL)
#PCA with broke stick rule
PCA_BS = read.csv(paste0("Results_Distance/PCA_Distance_BSR_IND_neutral_", project_name, ".csv"), row.names = 1)
PCA_BS[1:5, 1:5]
#PCA with 95% of variance explained
PCA_95 = read.csv(paste0("Results_Distance/PCA_Distance_95Var_IND_neutral_", project_name, ".csv"), row.names = 1)
PCA_95[1:5, 1:5]
#Dps
DPS = read.csv(paste0("Results_Distance/Dps_Shared_IND_neutral_", project_name, ".csv"), row.names = 1)
DPS[1:5, 1:5]

#C. Distance by clusters
#relatedness:
REL1 = read.csv(paste0("Results_Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP1.csv"), row.names = 1)
REL2 = read.csv(paste0("Results_Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP2.csv"), row.names = 1)
REL3 = read.csv(paste0("Results_Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP3.csv"), row.names = 1)
REL4 = read.csv(paste0("Results_Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP4.csv"), row.names = 1)


###3. TRANSFORM GENETIC DISTANCES DATAFRAME AS MATRIX AND DIST OBJECTS
#A. Create a list with all Relatedness distance dataframes
Related_files = list(REL1, REL2, REL3, REL4)

#B. Transform Relatedness data frames for populations into matrix
for (k in 1:length(Related_files)){
  rel_temp = Related_files[[k]]
  n_ind = nrow(unique(Related_files[[k]][1])) #number of individuals
  mt = matrix(0, n_ind, n_ind)
  mt[lower.tri(mt, diag=T)] = as.vector(rel_temp$RELATEDNESS_AJK_Yang)
  mt = mt + t(mt) 
  diag(mt) = diag(mt)/2 
  mt
  assign(paste0("mt_REL_", k), mt)
}

#C. Transform Relatedness data frame for all individual into matrix
n_ind = nrow(unique(REL[1])) #number of individuals
n_ind
mt_all = matrix(0, n_ind, n_ind)
mt_all[lower.tri(mt_all, diag=T)] = as.vector(REL$RELATEDNESS_AJK_Yang)
mt_all = mt_all + t(mt_all) 
diag(mt_all) = diag(mt_all)/2 
mt_all
#verify the matrix
class(mt_all)
dim(mt_all)
nrow(mt_all)
ncol(mt_all)


###4. LOAD GEOGRAPHIC INFORMATION:
#A. Create a data frame with the geographical coordenates:
coord_SN = snps_neutral@meta[,c(7,4:5)]
head(coord_SN)

#B. Set long and lat colunms
coordinates(coord_SN) = coord_SN[,c(2,3)]
projection(coord_SN) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

#C. Position of samples by populion by DAPC approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_tess))){
  pop = which(snps_neutral@meta$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}

#D. Subset coord by populations:
coord_POP1 = coord_SN[pop_TESS_1,]
coord_POP1
coord_POP2 = coord_SN[pop_TESS_2,]
coord_POP2
coord_POP3 = coord_SN[pop_TESS_3,]
coord_POP3
coord_POP4 = coord_SN[pop_TESS_4,]
coord_POP4

#D. Create a matrix of geographical distance
#by populations
coord_files = list(coord_POP1, coord_POP2, coord_POP3, coord_POP4)

for (l in 1:length(coord_files)){
  coord = coord_files[[l]]
  sdist_SN  = distm(coord, fun=distGeo)
  sdist_SN  = sdist_SN /1000
  assign(paste0("sdist_POP", l), sdist_SN)
}

#with all samples
sdist_SN  = distm(coord_SN, fun=distGeo)
rownames(sdist_SN) = snps_neutral@sample_id
colnames(sdist_SN) = snps_neutral@sample_id
sdist_SN
sdist_SN  = sdist_SN /1000
#verify the matrix
nrow(sdist_SN)
ncol(sdist_SN)
max(sdist_SN)
dim(sdist_SN)
class(sdist_SN)
sdist_SN


###5. PERFORM THE LOESS. IT WILL TAKE A LONG TIME DEPENDING ON NUMBER OF SAMPLES AND SNPS.
# mat1 is a geographic distance matrix (a full matrix NOT a dist object)
# mat2 is a kinship matrix (a full matrix NOT a dist object)
#All individuals
ot_SN = mDistoLoess(sdist_SN, mt_all, 1000) #numbers indicate the number of replications.

#POP1
ot_SN1 = mDistoLoess(sdist_POP1, mt_REL_1, 1000)

#POP2
ot_SN2 = mDistoLoess(sdist_POP2, mt_REL_2, 1000)

#POP3
ot_SN3 = mDistoLoess(sdist_POP3, mt_REL_3, 1000)

#POP4
ot_SN4 = mDistoLoess(sdist_POP4, mt_REL_4, 1000)


###6. PLOT THE RESULTS
#A. Create a data frame with the results
dd_SN = data.frame(cov = sdist_SN[lower.tri(sdist_SN)], obs = attr(ot_SN$disto, "obs"), lci = attr(ot_SN$disto, "CI")[1,], uci = attr(ot_SN$disto, "CI")[2,])

dd_pop1 = data.frame(cov = sdist_POP1[lower.tri(sdist_POP1)], obs = attr(ot_SN1$disto, "obs"), lci = attr(ot_SN1$disto, "CI")[1,], uci = attr(ot_SN1$disto, "CI")[2,])

dd_pop2 = data.frame(cov = sdist_POP2[lower.tri(sdist_POP2)], obs = attr(ot_SN2$disto, "obs"), lci = attr(ot_SN2$disto, "CI")[1,], uci = attr(ot_SN2$disto, "CI")[2,])

dd_pop3 = data.frame(cov = sdist_POP3[lower.tri(sdist_POP3)], obs = attr(ot_SN3$disto, "obs"), lci = attr(ot_SN3$disto, "CI")[1,], uci = attr(ot_SN3$disto, "CI")[2,])

dd_pop4 = data.frame(cov = sdist_POP4[lower.tri(sdist_POP4)], obs = attr(ot_SN4$disto, "obs"), lci = attr(ot_SN4$disto, "CI")[1,], uci = attr(ot_SN4$disto, "CI")[2,])

#B. Create a plot
LOESS_PLOT = ggplot(dd_pop4, aes(x=cov,y=obs)) + #change for cluster here
  geom_line(size=1.2) +
  geom_ribbon(aes(ymax = uci, ymin = lci), alpha = 0.3)+
  theme_minimal() + 
  geom_hline(yintercept = mean(ot_SN4$disto), lty = 3) + #change for cluster here
  ylab("Mean Relatedness") +
  xlab("Distance (Km)") +
  geom_rug(sides = "b", alpha = 0.02) +
  annotate("text", x = 2.2, y = 0.10, label = "", size=10) +
  theme(axis.title.y = element_text(size=20, color = "black", face = "bold"),
        axis.title.x = element_text(size=20, color = "black", face = "bold")) + 
  theme(axis.text = element_text(face = "bold", color = "black", size = 12))

#C. Save plot as a figure
#verify
plot(LOESS_PLOT)
#save as pdf
pdf(paste0("./Results_Spatial_Structure/LOESS/LOESS_POP4_TESS_", project_name, ".pdf"), family="Helvetica", onefile = F) #change for cluster here
plot.new()
plot(LOESS_PLOT)
dev.off()

##END
