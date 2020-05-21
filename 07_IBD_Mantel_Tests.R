###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
#####################   STEP 07: IBD BY MANTEL TEST   #########################

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
#A. ISOLATION-BY-DISTANCE ANALYSIS USING MANTEL TEST 

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
if("geosphere" %in% rownames(installed.packages()) == FALSE){install.packages("geosphere")
} else {print (paste0("'geosphere' has already been installed in library"))}
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("rgdal" %in% rownames(installed.packages()) == FALSE){install.packages("rgdal")
} else {print (paste0("'rgdal' has already been installed in library"))}
if("ade4" %in% rownames(installed.packages()) == FALSE){install.packages("ade4")
} else {print (paste0("'ade4' has already been installed in library"))}
if("usedist" %in% rownames(installed.packages()) == FALSE){install.packages("usedist")
} else {print (paste0("'usedist' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, usedist, raster, rgdal, geosphere, ade4)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_IBD/Mantel_Tests"))


#------------------------------------------------------------------------------
#                             1. Loading Files 
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "pilocarpus"


###2.1 LOAD VCF FILE WITH GEOGRAPHICAL INFORMATION: 
#A. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.


#------------------------------------------------------------------------------
#                              2. Genetic Data 
#------------------------------------------------------------------------------
###2.1. LOAD GENETIC DISTANCE:
#A. I used PCA_95 because it had similar results from PCA_BS and larger number means more differences

#B. Load distance file among all individuals calculate in STEP 4:
#PCA with 95% of variance explained
PCA_95 = read.csv(paste0("Results_Distance/PCA_Distance_95Var_IND_neutral_", project_name, ".csv"), row.names = 1)
PCA_95[1:5, 1:5]


###2.2. CREATE GENETIC 'DIST' OBJECTS FOR ALL SAMPLES: 
#A. All samples:
genDIST_all = as.dist(PCA_95)


###2.3. CREATE GENETIC 'DIST' OBJECTS BY CLUSTERS:
#A. Position of samples by populion by TESS approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_tess))){
  pop = which(snps_neutral@meta$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}

#B. Subset By Pop 1:
genDIST_pop1= dist_subset(genDIST_all, snps_neutral@sample_id[pop_TESS_1])
#verify dimensions
length(geoDIST_pop1)

#C. Subset By Pop 2:
genDIST_pop2= dist_subset(genDIST_all, snps_neutral@sample_id[pop_TESS_2])
#verify dimensions
length(genDIST_pop2)

#D. Subset By Pop 3:
genDIST_pop3= dist_subset(genDIST_all, snps_neutral@sample_id[pop_TESS_3])
#verify dimensions
length(genDIST_pop3)

#E. Subset By Pop 4:
genDIST_pop4= dist_subset(genDIST_all, snps_neutral@sample_id[pop_TESS_4])
#verify dimensions
length(genDIST_pop4)


#------------------------------------------------------------------------------
#                              3. Geographical Data 
#------------------------------------------------------------------------------
###3.1. LOAD GEOGRAPHIC INFORMATION:
#A. Create a data frame with the geographical coordenates:
coord_SN = snps_neutral@meta[,c(7,4:5)]
head(coord_SN)

#B. Set long and lat colunms
coordinates(coord_SN) = coord_SN[,c(2,3)]
projection(coord_SN) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

#C. Subset coord by populations:
coord_POP1 = coord_SN[pop_TESS_1,]
coord_POP1
coord_POP2 = coord_SN[pop_TESS_2,]
coord_POP2
coord_POP3 = coord_SN[pop_TESS_3,]
coord_POP3
coord_POP4 = coord_SN[pop_TESS_4,]
coord_POP4

###3.2. CREATE MATRIXES FOR GEOGRAPHICAL DISTANCE:
#by populations
coord_files = list(coord_POP1, coord_POP2, coord_POP3, coord_POP4)

for (l in 1:length(coord_files)){
  coord = coord_files[[l]]
  sdist_SN  = distm(coord, fun=distGeo)
  sdist_SN  = sdist_SN /1000
  assign(paste0("sdist_POP", l), sdist_SN)
}

# all samples
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


###3.3. CREATE 'DIST' OBJECTS FOR ALL SAMPLES AND BY CLUSTERS FROM MATRIXES: 
#A. All samples:
geoDIST_all = as.dist(sdist_SN)

#B. By Pop:
geoDIST_pop1 = as.dist(sdist_POP1)
geoDIST_pop2 = as.dist(sdist_POP2)
geoDIST_pop3 = as.dist(sdist_POP3)
geoDIST_pop4 = as.dist(sdist_POP4)


###3.4. VERIFY DIMENSION FOR THE 'DIST' OBJECTS:
length(geoDIST_all)
length(genDIST_all)

length(geoDIST_pop1)
length(genDIST_pop1)

length(geoDIST_pop2)
length(genDIST_pop2)

length(geoDIST_pop3)
length(genDIST_pop3)

length(geoDIST_pop4)
length(genDIST_pop4)


#------------------------------------------------------------------------------
#                              4. Mantel Tests 
#------------------------------------------------------------------------------
###4.1. PERFORM MANTEL TESTS
#A. By Species
mantel_all = mantel.rtest(genDIST_all, geoDIST_all, 10000)

#B. By populations
mantel_pop1 = mantel.rtest(genDIST_pop1, geoDIST_pop1, 10000)
mantel_pop2 = mantel.rtest(genDIST_pop2, geoDIST_pop2, 10000)
mantel_pop3 = mantel.rtest(genDIST_pop3, geoDIST_pop3, 10000)
mantel_pop4 = mantel.rtest(genDIST_pop4, geoDIST_pop4, 10000)

#C. Compile results and save as data.frame
results_mantel = matrix(NA, 5, 2)
colnames(results_mantel) = c("r", "p-value")
rownames(results_mantel) = c("All Samples", "POP1", "POP2", "POP3", "POP4")
results_mantel

results_mantel[1,] = cbind (mantel_all$obs, mantel_all$pvalue)
results_mantel[2,] = cbind (mantel_pop1$obs, mantel_pop1$pvalue)
results_mantel[3,] = cbind (mantel_pop2$obs, mantel_pop2$pvalue)
results_mantel[4,] = cbind (mantel_pop3$obs, mantel_pop3$pvalue)
results_mantel[5,] = cbind (mantel_pop4$obs, mantel_pop4$pvalue)

write.csv(results_mantel, paste0("./Results_IBD/Mantel_Tests/Mantel_Results_", project_name, ".csv"))

#D. Save Plot as PDF for all individuals
pdf("./Results_IBD/Mantel_Tests/Mantel_Plot_All.pdf", onefile=FALSE)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0, max(sdist_SN)), ylim=c(min(genDIST_all), max(genDIST_all)));

points(as.vector(geoDIST_all), as.vector(genDIST_all),  col = 'black', bg = 'white', cex = .9, pch=21)
abline(lm(as.vector(genDIST_all) ~ as.vector(geoDIST_all)), lwd=3, col='black')

axis(1, at=seq(0, max(sdist_SN), by=10), cex.axis=1);
axis(2, at=seq(round(min(genDIST_all),0), round(max(genDIST_all),0), by=10), cex.axis=1, las=1);
mtext(side=1, text='Geographic Distances (Km)', line=2.8, cex=1)
mtext(side=2, text='PC Distance (95%)', line=3.3 , cex=1)
dev.off()


#E. Save Plot as PDF for genetic clusters
pdf("./Results_IBD/Mantel_Tests/Mantel_Plot_POP.pdf", onefile=FALSE)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0, max(sdist_SN)), ylim=c(min(genDIST_all), max(genDIST_all)));

#By population. Keep tha same colors than other analyses:
col_tess= c("yellow", 'darkslategrey', 'red', 'blue')
legend_tess = c("POP1", "POP2", "POP3", "POP4")

points(as.vector(geoDIST_pop4), as.vector(genDIST_pop4),  col = 'black', bg = 'blue', cex = .9, pch=21)
points(as.vector(geoDIST_pop3), as.vector(genDIST_pop3),  col = 'black', bg = 'red', cex = .9, pch=21)
points(as.vector(geoDIST_pop1), as.vector(genDIST_pop1),  col = 'black', bg = 'yellow', cex = .9, pch=21)
points(as.vector(geoDIST_pop2), as.vector(genDIST_pop2),  col = 'black', bg = 'darkslategrey', cex = .9, pch=21)

abline(lm(as.vector(genDIST_pop1) ~ as.vector(geoDIST_pop1)), lwd=3, col='yellow')
abline(lm(as.vector(genDIST_pop2) ~ as.vector(geoDIST_pop2)), lwd=3, col='darkslategrey')
abline(lm(as.vector(genDIST_pop3) ~ as.vector(geoDIST_pop3)), lwd=3, col='red')
abline(lm(as.vector(genDIST_pop4) ~ as.vector(geoDIST_pop4)), lwd=3, col='blue')

axis(1, at=seq(0, max(sdist_SN), by=10), cex.axis=1);
axis(2, at=seq(round(min(genDIST_all),0), round(max(genDIST_all),0), by=10), cex.axis=1, las=1);
mtext(side=1, text='Geographic Distances (Km)', line=2.8, cex=1)
mtext(side=2, text='PC Distance (95%)', line=3.3 , cex=1)

legend("bottomright", legend = legend_tess, cex=1, bty="n",  col = "black", pt.bg = col_tess , pch= 21)
dev.off()


#------------------------------------------------------------------------------
#                       5. Permutation Tests by Population 
#------------------------------------------------------------------------------
###5.1. MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX - PERMUTATION
#A. Define populations positions to subset in a list
populations = list(snps_neutral@sample_id[pop_TESS_1], snps_neutral@sample_id[pop_TESS_2], snps_neutral@sample_id[pop_TESS_3], snps_neutral@sample_id[pop_TESS_4])

#B. Define number of populations
pop_counter = length(populations)

#C. Define a object to save the results
mantel_perm = matrix(NA,pop_counter,2)

#D. Perform the permutation
for(i in 1:pop_counter) {
  samples = !(snps_neutral@sample_id %in% populations[[i]])
  genDIST_temp = dist_subset(genDIST_all, snps_neutral@sample_id[samples])
  geoDIST_temp = dist_subset(geoDIST_all, snps_neutral@sample_id[samples])
  mt = mantel.rtest(genDIST_temp, geoDIST_temp, 10000)
  mantel_perm[i, ] = cbind (mt$obs, mt$pvalue)
}

mantel_perm
colnames(mantel_perm) = c("r", "p-value")
rownames(mantel_perm) = c("POP1", "POP2", "POP3", "POP4")
mantel_perm

#E. Save the results:
write.csv(mantel_perm, paste0("./Results_IBD/Mantel_Tests/Mantel_Results_Permutation_", project_name, ".csv"))

##END
