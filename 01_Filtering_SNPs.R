###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
########################## STEP 01: FILTERING SNPS ############################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###


#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##1. REFERENCES AND MORE INFORMATION:
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS TUTORIAL IS ORGANIZED IN DIFFERENT "Steps" AND INSIDE EACH STEP THERE ARE "Actions" DENOTED BY NUMBERS (#1) AND INSIDE EACH ACTION COULD EXIST "Operations" INDICATED BY LETTES (#A.)
#D. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


##2. INPUTS FOR THIS TUTORIAL:
#A. ALIGNED AND CLEANED ".vcf" OR ".gvcf" FILE AFTER POPULATIONS STEP IN STACKS OR FROM ANOTHER SOFTWARE.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IT IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#C. THE FILE ".CSV" WITH THE GEOGRAPHICAL INFORMATION OF GENETIC SAMPLES IN DECIMALS DEGREES PREFERENCIALLY AND ONE COLUMN WITH LOCALITIES IDs. DIFFERENT INDIVIDUALS MUST HAVE DIFFERENT COORDINATES. IF TWO OR MORE SAMPLES WERE COLLECTED IN THE SAME LOCALITY, CHANGE THE VALUE OF 4TH OR 5TH DECIMAL NUMBER TO SET A DIFFERENCE OF 1 OR 5 METERS. SOME ANALYSES COULD RETURN ERRORS IF THE DIFFERENCES IN METERS BETWEEN INDIVIDUALS IS 0. 


##3. GOALS FOR THIS STEP:
#A. FILTER RAW SNPS DATASET INTO NEUTRAL SNPS DATASET. 
#B. CREATE A DATASET THAT WILL BE USED IN RDA AND LFMM ANALYSES.


##4. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H
#B. CREATE A FOLDER NAMED "vcf" INSIDE YOUR WORKING DIRECTORY AND COPY THE .vcf OR .gvcf FILE THERE.
#C. CREATE A FOLDER NAMED "coords" INSIDE YOUR WORKING DIRECTORY AND COPY THERE, THE FILE .csv WITH GEOGRAPHICAL INFORMATION OF THE SAMPLES


##5. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##6. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("functions_LanGen.R")


##7. INSTALL AND LOAD THE PACKAGES
#A. Install packages automatically
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
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")
} else {print (paste0("'tidyverse' has already been installed in library"))}
if("vcfR" %in% rownames(installed.packages()) == FALSE){install.packages("vcfR")
} else {print (paste0("'vcfR' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("poppr" %in% rownames(installed.packages()) == FALSE){install.packages("poppr")
} else {print (paste0("'poppr' has already been installed in library"))}
if("pcadapt" %in% rownames(installed.packages()) == FALSE){install.packages("pcadapt")
} else {print (paste0("'pcadapt' has already been installed in library"))}
if("qvalue" %in% rownames(installed.packages()) == FALSE){BiocManager::install("qvalue")
} else {print (paste0("'qvalue' has already been installed in library"))}
if("psych" %in% rownames(installed.packages()) == FALSE){install.packages("psych")
} else {print (paste0("'psych' has already been installed in library"))}
if("tess3r" %in% rownames(installed.packages()) == FALSE){devtools::install_github("bcm-uga/TESS3_encho_sen")
} else {print (paste0("'tess3r' has already been installed in library"))}

#B. Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, tidyverse, vcfR, dartR, poppr, adegenet, LEA, pcadapt, qvalue, psych, tess3r)


##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./adapt_var_mapping/Filtering", "./Results_Metafiles", "./Results_Filters/PCA", "./Results_Filters/PCadapt_FST", "./Results_Filters/sNMF_FST", "./Results_Filters/TESS_FST"))


#------------------------------------------------------------------------------
#                       1. Load Files and Verify Data Quality 
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT AND THE PATH FOR VCF/GVCF AND COORDS FILES:
#A. vcf/gvcf file path:
vcf_file = "vcf/pilocarpusredo.vcf"

#B. Project name:
project_name = "pilocarpus"

#C. Coordenates file path and load it:
coord_file = "coords/pilocarpus_coords.csv"
#load the file
coords = read.csv(coord_file)
#verify if it loaded correctly
head(coords)
tail(coords)

#D. Define the columns positions for localities ID (local), sample ID (sample_ID), longitude, and latitude. Add more columns if you need.
local = 1
sample_ID = 2
longitude = 3
latitude = 4


###1.2. LOAD THE VCF FILE
#A. VCF
snps_raw <-  vcfLink(vcf_file, overwriteID=T)
snps_raw

#B. Remove loci that are not SNPs
snps <- Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2), indels="remove")

#C. Basic stats sample size (SZ)
raw_SZ = capture.output(VCFsummary(snps_raw)) 
filtered_SZ = capture.output(VCFsummary(snps)) 
raw_SZ #291 individuals and 36672 SNPs.
filtered_SZ #291 individuals and 35214 SNPs.
snps@sample_id
snps@site_id
snps@meta
Chrom(snps) ##Chromosome, possitions, and IDs


###1.3. REMOVE INDIVIDUALS IF IT IS NECESSARY
#A. Remove one individual
#snps@meta$sample_name
#UNIND <- snps@sample_id[snps@sample_id != "9202_rep_sorted"]
#snps_unind <- Subset(snps, samples=UNIND)
#snps_unind@meta
#VCFsummary(snps_unind) ## 290 individuals and 35214 SNPs.

#B. Remove more than one individual. In the case Granito locality, the rep individual (9202_rep), and individuals (9541, 9561):
snps@meta$sample_name
remove_ind = c("9202_rep_sorted", "9541_sorted", "9561_sorted", "9521_sorted","9522_sorted", "9523_sorted", "9526_sorted","9530_sorted","9531_sorted", "9532_sorted","9533_sorted","9537_sorted","9538_sorted","9540_sorted")
UNIND <- snps@sample_id[!snps@sample_id %in% remove_ind]
snps_unind <- Subset(snps, samples=UNIND)
#verify sample size
length(snps@meta$sample_name)-length(remove_ind)
length(snps_unind@meta$sample_name)
VCFsummary(snps_unind) ## 277 individuals and 35214 SNPs.

#C. If you do not need to remove individuals:
#snps_unind = snps
#VCFsummary(snps_unind) ## 277 individuals and 35214 SNPs.


###1.4. VERIFY THE QUALITY OF DATA
#A. See GenotypeMatrix
genotypes <- GenotypeMatrix(snps_unind)
genotypes[1:10, 1:10] ## -1 is missing; otherwise, gives the number of derived alleles in the genotype -- ie. a 0 and 2 are both homozygotes

#B. Look at depth, quality, HWE, HE, allele frequencies, and Pi
site.depth <- Query(snps_unind, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 63.40 / Median = 47.30
hist(site.depth$MEAN_DEPTH, breaks=30) 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <100]) ### >20 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <250]) ### <200

quality <- Query(snps_unind, type="site-quality") 
summary(quality$QUAL) #Mean = 46.47 / Median = 48.63
hist(quality$QUAL) ## >30

PI <- Query(snps_unind, type="site-pi")
mean(PI$PI) ## Mean nucleotide divergency per-site ###0.241
hist(PI$PI)

HWE <- Query(snps_unind, type="hardy")
summary(HWE$P_HWE) #Mean = 0.129 / Median = 0.0012
hist(HWE$P_HWE)
hist(HWE$P_HWE[HWE$P_HWE<0.0001])

HE <- Query(snps_unind, type="het")
hist(HE$O.HOM) ## O.HOM, E.HOM, N_SITES, F
hist(HE$E.HOM) 
hist(HE$N_SITES) 
hist(HE$F) ## Inbreeding coefficient

freq <- Query(snps_unind, type="freq2")
hist(freq$X.FREQ.)
hist(freq$X.FREQ.1)
head(freq)

#C. Missing per locus
Missing <- apply(GenotypeMatrix(snps_unind), 2, function(x) sum(x < 0)/length(x)*100)
summary(Missing) ## Max missing = 54.9%
hist(Missing)

#D. Missing per individual
Missing_ind <- apply(GenotypeMatrix(snps_unind),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 50.3%
hist(Missing_ind)


###1.5. ADD GEOGRAPHICAL INFORMATION AND MISSING DATA BY INDIVIDUALS TO A VCF FILE:
#A. Save missind data information as data frame
Missingind.df = as.data.frame(Missing_ind)
#remove string "_sorted" from samples ID if you need it
row.names(Missingind.df) = str_remove(row.names(Missingind.df), pattern="_sorted")
#add samples ID to the data frame
Missingind.df$ID = row.names(Missingind.df)
Missingind.df

#B. Remove geographical information of individuals that you excluded in step #3 
missing_coord = coords[coords[,sample_ID] %in% Missingind.df$ID,]
head(missing_coord)
tail(missing_coord)
length(missing_coord[,1])

#C. When coords are in UTM you need to transform to lon/lat in decimals. Use this code:
#C1. remove NAs if they are present.
#missing_coord = missing_coord[!is.na(missing_coord[,longitude]),]

#C2. Add a datum
#My UTM data are in SAD69 UTM zone 22S. Search in internet the EPSG code for that area. For example https://epsg.io/ 
#In the case is EPSG:32722.
#coord_gen_pts = missing_coord[,c(longitude, latitude)] #select lat and long
#coordinates(coord_gen_pts) <- coord_gen_pts
#projection(coord_gen_pts) = crs("+init=epsg:32722") #define utm datum only for coord cols

#C3. Convert:
#coords_gen_dec = spTransform(coord_gen_pts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"))
#create a data.frame with converted coords:
#x = missing_coord[, c(local, sample_ID)]
#y = as.data.frame(coords_gen_dec)
#y = y[,c(longitude, latitude)]
#coords_gen = cbind(x, y)
#colnames(coords_gen) = c("local_ID", "samples", "long", "lat")
#coords_dec = coords_gen
#head(coords_dec)
#tail(coords_dec)
#missing_coord = coords_dec

#D. Sort the geographical information file according to the sample order of the missing data file:
missing_coord_sorted = missing_coord[order(match(missing_coord[,sample_ID], Missingind.df$ID)),]
head(missing_coord_sorted)

#E. Check if there are some difference between files. If they are corrected, function will return "character (0)"
setdiff(Missingind.df$ID, missing_coord_sorted[,sample_ID])

#F. Check if identical samples were selected in the same order. If they are corrected, function will return "TRUE"
identical(as.character(Missingind.df$ID), as.character(missing_coord_sorted[,sample_ID]))

#G. Create a data frame merging missing and geographical information
missing_local = cbind(missing_coord_sorted[,c(local, longitude, latitude)], Missingind.df)
head(missing_local)
tail(missing_local)

#H. Add coords and missing information to vcf file. Be careful to the order!
#verify order. If they are right, function will return "TRUE"
#if you have "_sorted" add in the samples names:
identical(as.character(missing_local$ID), as.character(str_remove(snps_unind@meta$sample_name, pattern="_sorted")))
#if you do not have "_sorted" add in the samples names:
identical(as.character(missing_local$ID), as.character(snps_unind@meta$sample_name))
#add the information. You can add more columns if you need, like ecological traits or other grouping classification 
snps_unind@meta$local_ID = missing_local[,1] #local ID
snps_unind@meta$longitude = missing_local[,2] #longitude
snps_unind@meta$latitude = missing_local[,3] #latitude
snps_unind@meta$missing_PC = missing_local[,4] #% of missing data
snps_unind@meta$ind_ID = missing_local[,5] #sample ID without "_sorted"
#verify the file
head(snps_unind@meta) #verify
tail(snps_unind@meta) #verify

#I. Save metafile with missing data and geographical information:
write_csv(snps_unind@meta, paste0("./Results_Metafiles/missing_data_coords_", project_name, ".csv"))


###1.6. REMOVE INDIVIDUALS WITH LARGE AMOUNT OF MISSING GENOTYPES
#A. Test individuals with large amounts of missing (set to max 30% missing). Try different thresholds!  
missing_data = snps_unind@meta[snps_unind@meta$missing_PC > 30,]
#create a data frame for graphical output
df = as.data.frame(tapply(missing_data$local_ID,missing_data$local_ID,length))
df$local = rownames(df)
colnames(df) = c("sample_size", "locality")
df= df[!is.na(df$sample_size),]
#analyze in a graph where are from the samples you will remove. 
pie(df$sample_size, df$locality, main="Missing Data per Locality")
#analyze in a table where are from the samples you will remove.
df

#B. After testing different thresholds, set the missing data threshold and remove the individuals:
value = 70
indtokeep <- snps_unind@meta[snps_unind@meta$missing_PC <= value,]
snps_ind_pc <- Subset(snps_unind, samples = indtokeep$sample_name)
snps_ind_pc@meta
VCFsummary(snps_unind) ## 277 individuals and 35214 SNPs.
removed_missing_SZ = capture.output(VCFsummary(snps_ind_pc)) 
removed_missing_SZ ## 277 individuals and 35214 SNPs.


#------------------------------------------------------------------------------
#                          2. Filtering Adaptive Loci 
#------------------------------------------------------------------------------
###2.1. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE:
#A. From dataset with low missing by individuals
snps_fil_low = Filter(snps_ind_pc, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200)) 
VCFsummary(snps_fil_low) #277 individuals and 25241 SNPs.

#B. From dataset with all individuals
snps_fil_high <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200)) 
VCFsummary(snps_fil_high)

#C. Choose one of the datasets (A ou B):
snps_fil = snps_fil_low
adaptativeFilter_SZ = capture.output(VCFsummary(snps_fil))

#D. Define R² value:
r2 = 0.4


###2.2. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:
#A. remove snps with R² value
#ld_within <- Linkage(snps_fil, type="geno-r2", linkageOptions(min.r2= r2)) 
ld_within <- read.csv(paste0("adapt_var_mapping/Filtering/ld_within_", r2 , "_test2.csv")) #load this file if it was saved before
head(ld_within)
hist(ld_within$R.2)
#write.csv(ld_within, file=paste0("adapt_var_mapping/Filtering/ld_within_", r2 , "_test2.csv"))

#B. Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within$ID1
nold_snps <- snps_fil@site_id[!(snps_fil@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil, sites=nold_snps)
adaptativeLD_SZ = capture.output(VCFsummary(snps_fil_ld))
adaptativeLD_SZ #277 individuals and 19025 SNPs.


###2.3. VERIFY THE REAL MISSING DATA VALUE IN THE DATASET AFTER REMOVING CORRELATED SNPS:
#A. Missing per individual:
Missing_ind <- apply(GenotypeMatrix(snps_fil_ld),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 50.16%
hist(Missing_ind)

#B. Missing per locus:
Missing <- apply(GenotypeMatrix(snps_fil_ld), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 29.96%
hist(Missing)


###2.4. SAVE THE .VCF FILE FILTERED WITH ADAPTATIVE SNPS
Save(snps_fil_ld, paste0("vcf/", project_name,"_filtered_adaptative.vcf"))


#------------------------------------------------------------------------------
#                        3. Filtering Neutral Loci - Part I
#------------------------------------------------------------------------------
###3.1. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE, AND HARDY-WEINBERG EQUILIBRIUM (HWE):
#A. From dataset with low missing by individuals
snps_fil_hwe_low <- Filter(snps_ind_pc, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
VCFsummary(snps_fil_hwe_low) #277 individuals and 13958 SNPs.

#B. From dataset with all individuals
snps_fil_hwe_high <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
VCFsummary(snps_fil_hwe_high)  ## 277 individuals and 13958 SNPs.

#C. Choose one of the datasets (A ou B):
snps_fil_hwe = snps_fil_hwe_low
neutralFilter_SZ = capture.output(VCFsummary(snps_fil_hwe))

#D. Define R² value:
r2 = 0.4


###3.2. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:
## If you have more than one population, you can use this to identify SNPs deviating from HW equilibrium within each population, and then removes those SNPs that are in desequilibrium in all populations. You just need to subset your samples:

#A. remove snps with R² value
#ld_within <- Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=r2))
ld_within <- read.csv(paste0("adapt_var_mapping/Filtering/ld_within_",r2, "_hwe_test2.csv")) #load this file if it was saved before
head(ld_within)
hist(ld_within$R.2)
#write.csv(ld_within, file=paste0("adapt_var_mapping/Filtering/ld_within_",r2, "_hwe_test2.csv"))

#B. Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within$ID1
nold_snps <- snps_fil_hwe@site_id[!(snps_fil_hwe@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil_hwe, sites=nold_snps) # Keep snps that are not in LD.
neutralLDWithin_SZ = capture.output(VCFsummary(snps_fil_ld))
neutralLDWithin_SZ ##277 individuals and 11104 SNPs.


###3.3. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) BETWEEN CONTIGS:
#A. remove snps with R²
#ld_between <- Linkage(snps_fil_ld, type="interchrom-geno-r2", linkageOptions(min.r2=r2)) 
ld_between <- read.csv(paste0("adapt_var_mapping/Filtering/ld_between_", r2, "_hwe_test2.csv")) #load this file if it was saved before
head(ld_between)
hist(ld_between$R.2)
#write.csv(ld_between, file= paste0("adapt_var_mapping/Filtering/ld_between_", r2, "_hwe_test2.csv"))

#B. Select one set of the correlated snps (ID1 or ID2)
ld2_snps <- ld_between$ID1
nold2_snps <- snps_fil_ld@site_id[!(snps_fil_ld@site_id %in% ld2_snps)]
snps_fil_ldF <- Subset(snps_fil_ld, sites=nold2_snps) # Keep snps that are not in LD.
neutralLDBetween_SZ = capture.output(VCFsummary(snps_fil_ldF)) 
neutralLDBetween_SZ ##277 individuals and 6602 SNPs.


###3.4. VERIFY THE QUALITY OF THE FILTERED DATASET:
#A. Quality indexes:
site.depth2 <- Query(snps_fil_ldF, "site-mean-depth")
quality2 <- Query(snps_fil_ldF, "site-quality")
HWE2 <- Query(snps_fil_ldF, type="hardy")

hist(site.depth2$MEAN_DEPTH, breaks=30)
hist(site.depth2$MEAN_DEPTH[site.depth2$MEAN_DEPTH <200])
hist(quality2$QUAL)
hist(HWE2$P_HWE[HWE2$P_HWE<200])

#B. Verify the real missing data per individual:
Missing_ind <- apply(GenotypeMatrix(snps_fil_ldF),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 47.73%
hist(Missing_ind)

#C. Verify the real missing data per locus:
Missing <- apply(GenotypeMatrix(snps_fil_ldF), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 29.96%
hist(Missing)

VCFsummary(snps_unind) #277 individuals and 35214 SNPs.
VCFsummary(snps_fil_ldF) ##277 individuals and 6602 SNPs.

###3.5. SAVE THE .VCF FILE WITH ONLY NEUTRAL SNPS:
Save(snps_fil_ldF, paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"))


#------------------------------------------------------------------------------
#             4. Testing the Missing Data Influence on Genetic Structure 
#------------------------------------------------------------------------------
###4.1. TEST IF MISSING DATA IS AFFECTING THE GENETIC STRUCTURE ESTIMATION: 
#A. Load neutral and adaptative .vcf files from filtering step:
vcf_adaptative = read.vcfR(paste0("vcf/", project_name,"_filtered_adaptative.vcf"), verbose = FALSE)
vcf_neutral = read.vcfR(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), verbose = FALSE)

#B. Create a data frame with missing data values per individual based on snps_unind or snps_ind_pc if you removed some individual
#if you removed individuals in Action #6
missing_values = snps_ind_pc@meta$missing_PC
#if you did not remove individuals in Action #6
#missing_values = snps_unind@meta$missing_PC
#verify
length(missing_values)

#C. Convert "VCF" to "GENIND" for both datasets
genind_adap = convert_vcf2genind(vcf_adaptative)
genind_neut = convert_vcf2genind(vcf_neutral)

#D. Remove missing data in both datasets:
genind_adap_0 = missingno(genind_adap)
genind_adap_0 #verify number of SNPs = 2,921 SNPs

genind_neut_0 = missingno(genind_neut)
genind_neut_0 #verify number of SNPs =  927 SNPs

#E. Defining all datasets, names, and the vector with missing values. Define a threshold in % for missing data. Samples > threshold for missing data will be plot in red.
datasets = c(genind_adap, genind_adap_0, genind_neut, genind_neut_0)
names = c("Adapative_all","Adapative_0", "Neutral_all", "Neutral_0")
threshold = 20
missing_values = missing_values

#F. Run PCA by dataset
pca_missing (datasets, names, threshold, missing_values)

#G. Verify in the PCA graphs if there are clusters in your samples, if the number of clusters changed with or without missing data and if any individual in red has changed its position among the clusters. Greater spreading in data with more missing data and SNPs is common. Be careful! Changes could be due to a drop in the number of SNPs. Data with 1,000 SNPs are more robust to keep the genetic structure pattern. If necessary, re-filter and remove unstable individuals or reduce the allowed missing data threshold in Action #6, #7, and #11.


#------------------------------------------------------------------------------
#                       5. Filtering Neutral Loci - Part II - PCadapt 
#------------------------------------------------------------------------------
# This code serves to carry out PCAdapt to outlier detection.
# This script was first developed by Brenna Forester and adapted to this data.
# And also using the tutorial of https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
# PCadapt is a differentiation-based outlier method (no environmental data required)
# It can be run on individuals or pool-seq data.
# It does not require the identification of discrete populations, and performs well
# even in cases of hierarchical population structure and individual admixture.
# Like almost all differentiation-based methods, it is a univariate test which means
# we need to correct for multiple testing.
# You can use data with missing value, do not need to impute the missing data.
#---------------------------------------------------------------------------------

###5.1. CREATE A PCadapt MATRIX:
#A1. Load the vcf file "neutral_partial":
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
#A2. Or the snps_fil_ldF if it is on Global Environment
snps = snps_fil_ldF
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#B. Convert VCF in a Matrix
genmat = GenotypeMatrix(snps)
genmat[1:10,1:10]

#C. Amount of missing data in the matrix:
genmat[genmat == "-1"] <- NA
dim(genmat)
((sum(is.na(genmat)))/(dim(genmat)[1]*dim(genmat)[2]))*100
#12.65 %  #Amount of missing data

#D. Convert Matrix to LFMM
lfmm_input = Lfmm(snps, output.file=paste0("./Results_Filters/PCadapt_FST/", project_name, "_filtered_snps_lfmm_format.lfmm"))
lfmm_input = read.lfmm(paste0("./Results_Filters/PCadapt_FST/", project_name, "_filtered_snps_lfmm_format.lfmm"))
class(lfmm_input)

#E. Convert LFMM to a PCadapt matrix
gen.pcadapt = read.pcadapt(lfmm_input, type="lfmm")


###5.2. RUN PCadapt:
#A. Choosing the number K of Principal Components. Number of K is the number of PC and possible clusters:
pcadapt.test = pcadapt(gen.pcadapt, K=10, ploidy=2)
#"Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)
plot(pcadapt.test, option="screeplot") #2 PCs
plot(pcadapt.test, option="scores") #PC1 and PC2
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4
plot(pcadapt.test, option = "scores", i = 5, j = 4) #PC5 and PC4
# With PC3, PC4, PC5 the groups have changed. So we can use 2 or 1PC, K = 2 or 1
plot(pcadapt.test, option = "scores", i = 5, j = 6) #PC5 and PC6

#B. Save figure as pdf files:
pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_screenplot.pdf"), onefile = T)
plot(pcadapt.test, option="screeplot")
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_scores_PC1-PC2.pdf"), onefile = T)
plot(pcadapt.test, option="scores") #PC1 and PC2
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_scores_PC2-PC3.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_scores_PC3-PC4.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_scores_PC4-PC5.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 4, j = 5) #PC4 and PC5
dev.off()

#C. K = 2. Run PCadapt again. Test K and see the p-values graphs
pcadapt.test = pcadapt(gen.pcadapt, K=2, ploidy=2, min.maf=0.05, method="mahalanobis")
summary(pcadapt.test)
# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space. "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.

#D. Graphical tools.
plot(pcadapt.test, option = "manhattan")
plot(pcadapt.test, option = "qqplot")
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.test, option = "stat.distribution")

#save graphs
pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_manhattan.pdf"), onefile = T)
plot(pcadapt.test, option = "manhattan")
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_qqplot.pdf"), onefile = T)
plot(pcadapt.test, option = "qqplot")
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_p-values.pdf"), onefile = T)
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()

pdf(paste0("./Results_Filters/PCadapt_FST/", project_name, "_stat_distribution.pdf"), onefile = T)
plot(pcadapt.test, option = "stat.distribution")
dev.off()


###5.3. CHOOSING A METHOD TO IDENTIFY OUTLIERS:
#A. Choosing a cutoff for outlier detection
pcadapt.test$gif # 1.38718

#-----------------------------------------------------------------------------------------
# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.
#
# GIF of 1=well calibrated, >1 =liberal (too many small p-values), <1=conservative (too few small p-values) # Note: GIFs > 2 indicate poor test calibration.
#
#For a given alpha (real valued number between 0 and 1), SNPs with q-values less than alpha will be considered as outliers with an expected false discovery rate bounded by alpha. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs, for an expected false discovery rate lower than 10%.
#-----------------------------------------------------------------------------------

#A1. q-values
#qval = qvalue(pcadapt.test$pvalues)$qvalues
#alpha = 0.1
#outliers = which(qval < alpha)
#length(outliers) #It will be eliminated XXX SNPs

#A2. Benjamini-Hochberg Procedure
padj = p.adjust(pcadapt.test$pvalues, method="BH")
alpha = 0.1
outliers = which(padj < alpha)
length(outliers) #It will be eliminated 9 SNPs

#A3. Bonferroni correction
#padj2 = p.adjust(pcadapt.test$pvalues, method="bonferroni")
#alpha = 0.1
#outliers = which(padj2 < alpha)
#length(outliers) #It will be eliminated XXX SNPs

# q-values > BH > bonferroni. Choose one for outliers

#B. Selecting candidate outlier 
outlier_snps = Subset(snps, sites=outliers)
outlier_snps@site_id ##These are all the Outlier SNPs.
length(outlier_snps@site_id)

#C. Exclude candidate outlier
snps_candidate = Subset(snps, sites= outlier_snps@site_id)
snps_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_candidate)

candidates_fst = snps_candidate@site_id
All_snp = snps@site_id
N_snp = All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_neutral = Subset(snps, sites=N_snp)
VCFsummary(snps_neutral)
#verify if the exclusion was correct
length(N_snp)
length(snps@site_id)-length(outliers) 
length(snps_neutral@site_id)

#D. Save neutral snp dataset
neutral_after_fst_PCadapt = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_PCadapt #"277 individuals and 6217 SNPs."

Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_PCadapt.vcf"))


#------------------------------------------------------------------------------
#                      6. Filtering Neutral Loci - Part II - sNMF 
#------------------------------------------------------------------------------
# sNMF - sparse Nonnegative matrix factorization (NMF) for Estimation of Individual Ancestry Coefficients similar to STRUCTURE and ADMIXTURE.
# STRUCTURE and ADMIXTURE assumptions include absence of genetic drift and Hardy–Weinberg and linkage equilibrium in ancestral populations. sNMF was more appropriate to deal with inbred lineages. sNMF enabled the estimation of homozygote and heterozygote frequencies and avoided Hardy–Weinberg equilibrium assumptions. Cross-entropy criterion indicated better predictive results for sNMF than for ADMIXTURE. You can use impute() to predict Ancestry Coefficients for missing values.
#-------------------------------------------------------------------------------

###6.1. CREATE A geno OBJECT:
#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#B. Convert VCF to geno object. You need to specify the output file. It will automatically subset the vcf file and assign it as a new object.
snps = Geno(snps, output.file = paste0("./Results_Filters/sNMF_FST/", project_name, "_filtered_neutral_partial.geno"))
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#C. Choose the 6 values of alpha. Alpha is the value of the regularization parameter (by default: 10). The results can depend on the value of this parameter, especially for small data sets. Less than 10,000 SNPs you can use values from 1000 to 10000. More than 10,000 SNP between 1 to 2,000. You can test different values.
# >= 10,000 SNPs (lower values of alpha)
#alpha_values = c(1, 10, 100, 500, 1000, 2000)
# < 10,000 SNPs (higher values of alpha)
#alpha_values = c(1000, 2000, 4000, 6000, 8000, 10000)
# costumize values of alpha
alpha_values = c(10, 100, 500, 1000, 2000, 4000)


#D. Create folders for alpha values and copy .geno object in each folder:
for (i in alpha_values){
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i)
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("./Results_Filters/sNMF_FST/", project_name, "_filtered_neutral_partial.geno"), path )
}

#E. Set parameters to run SNMF (LEA) using different alpha values.
K = c(1:10) # set the number of K to be tested
replications = 5 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel

###6.2. RUN sNMF (LEA)
#A.Run a loop for all alpha values 
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop), pro_snmf)
}

#B. To load the SNMF projects. This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop), pro_snmf)
}

#C. Summary of the project
summary(project1)
summary(project2)
summary(project3)
summary(project4)
summary(project5)
summary(project6)

#D. View cross-Entropy plot
PlotK(project1) #5
PlotK(project2) #5
PlotK(project3) #5
PlotK(project4) #4
PlotK(project5) #5
PlotK(project6) #5
#another way to view the results
#plot(project1, lwd = 5, col = "red", pch=1)
#plot(project2, lwd = 5, col = "red", pch=1)
#plot(project3, lwd = 5, col = "red", pch=1)
#plot(project4, lwd = 5, col = "red", pch=1)
#plot(project5, lwd = 5, col = "red", pch=1)
#plot(project6, lwd = 5, col = "red", pch=1)

#E. Save Cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results_Filters/sNMF_FST/Cross_Entropy_sNMF_Alpha_",  i, ".pdf"), onefile = F)
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
}

#F. Select optimal K value
optimal_K = 5

#G. ATTENTION, if your dataset is K = 1 force a K = 2 to be able to filter SNPs with FST outliers.

#H. Select best run (lowest cross-entropy)  
best_run = Best.run(nrep=replications, optimalK=optimal_K, p1=project1, p2=project2, p3=project3, p4=project4, p5=project5, p6=project6)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("./Results_Filters/sNMF_FST/Alpha", alpha_values[as.numeric(best_run_split[6])],"/", project_name, "_filtered_neutral_partial.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])


###6.3. IDENTIFY AND REMOVE FST OUTLIERS
#--------------------------------------------------------------------------
#We need to find a number of k to perform analyses and remove the adaptative loci.
#incorrect model specification can generate large numbers of false-positive associations
# A key condition for applying false discovery rates (FDR) control algorithms is that the test P-values behave as uniformly distributed random variables when the null hypothesis is correct.
#In this case, the empirical distribution of P-values will consist of a mixture of a uniform distribution and a peaky distribution showing a peak close to 0. The uniform distribution corresponds to neutral loci, whereas the interesting loci are found under the peak (Efron 2004).
#According to standard population genetic theory, the definition of FST , which is based on Wright’s work, corresponds to the amount of variance in allele frequency that can be explained by population structure (Wright 1951). FST correspond to estimates of correlation coefficients.
#For appropriately calibrated methods, the inflation factor is expected to be close to one, and the histogram of test significance values is expected to display a flat shape.
#It should be noted that assessing inflation factors is always subjective, as the proximity to the value of one will vary, usually by increasing with the sample size and the background levels of FST or R2 . To be justified in recalibrating the significance values, it is advisable to report inflation factors prior to adjustment and histograms following adjustment.
#-------------------------------------------------------------------------

#A. Compute the FST statistics using best run
FST = fst(project, run, optimal_K) # you need at least 2 populations for a population-based test, so K>1.

#B. Compute the GIF - genomic inflation factor
lambda <- GIF(project, run, optimal_K, fst.values=FST)
lambda

#C. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

#D. Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/3.5, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
hist(adj.p.values, col = "green")
C_fst <- candidates(alpha=0.1, adj.p.values)
ManPlot(adj.p.values, C_fst,"Fst")

#E. after you choose one, save the result
pdf("./Results_Filters/sNMF_FST/sNMF_FST_Outliers_adj_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()

#F. Candidate loci for FDR control: Benjamini-Hochberg at level q
C_fst <- candidates(alpha=0.1, adj.p.values)

#G. save the Manhatan plot
pdf("./Results_Filters/sNMF_FST/sNMF_FST_Outliers_Manhatan_plot.pdf", onefile = T) #add the number of snps removed in the name of pdf
ManPlot(adj.p.values, C_fst, paste0("Fst - Removing ",  length(C_fst), " SNPs"))
dev.off()

#H. Exclude candidate FST outlier
snps_fil_ldF_candidate <- Subset(snps, sites=C_fst)
snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_neutral <- Subset(snps, sites=N_snp)
VCFsummary(snps_neutral) #277 individuals and 4941 SNPs.
length(N_snp)
length(snps@site_id)-length(C_fst) 
length(snps_neutral@site_id)

neutral_after_fst_sNMF = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_sNMF #277 individuals and 4941 SNPs.

#I. Save the vcf
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_sNMF.vcf"))


#------------------------------------------------------------------------------
#                    7. Filtering Neutral Loci - Part II - TESS3 
#------------------------------------------------------------------------------
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html

###7.1. INPUT FILES FOR TESS:
#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#B. Create a Genotype matrix
genotypes = GenotypeMatrix(snps) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing;
class(genotypes)
genotypes = replace(genotypes, genotypes == -1, NA)

#C. Create a Matrix with long and lat 
coordinates = snps@meta[,4:5]
class(coordinates)
coordinates = data.matrix(coordinates, rownames.force = NA)
class(coordinates)
#verify the coords
plot(coordinates, pch = 19, cex = .5, xlab = "Longitude", ylab = "Latitude")


###7.2. RUNNING THE TESS3R FUNCTION:
#A. Costumize values of lambda
lambda_values = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)

#B.Run a loop for all alpha values
set.seed(13)
for (i in lambda_values){
    tess3.ls = tess3(genotypes, coord = coordinates, K = 1:10, mask = 0.05, lambda = i,
                   method = "projected.ls", max.iteration = 5000, rep = 5,
                   ploidy = 2, openMP.core.num = 4)
    save(tess3.ls, file = paste0("./Results_Filters/TESS_FST/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
    assign(paste0("TESS3_Lambda_", i), tess3.ls)
    pdf(paste0("./Results_Filters/TESS_FST/TESS3_PlotK_CrossValidation_Lambda_", i, ".pdf"), onefile =F)
    plot.new()
    plot(tess3.ls, pch = 19, col = "blue", type="b",lwd=2,
         crossvalid = T, crossentropy = T,
         xlab = "Number of ancestral populations", cex.lab=1.5,
         ylab = "Cross-validation score")
    dev.off()
}

#C. Plot results for all lambda values:
#Lambda 0: k = 6
plot(TESS3_Lambda_0, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 0.25: k = 5
plot(TESS3_Lambda_0.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 0.5: k = 5
plot(TESS3_Lambda_0.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 0.75: k = 5
plot(TESS3_Lambda_0.75, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 1: k = 4
plot(TESS3_Lambda_1, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 1.25: k = 7
plot(TESS3_Lambda_1.25, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#Lambda 1.5: k = 5
plot(TESS3_Lambda_1.5, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")

#D. Choose best K:
optimal_K = 5
#ATTENTION, if your dataset is K = 1 force a K = 2 to be able to filter SNPs with FST outliers.

min(TESS3_Lambda_0[[optimal_K]]$crossvalid.crossentropy) #0.3630269
min(TESS3_Lambda_0.25[[optimal_K]]$crossvalid.crossentropy) #0.3599267
min(TESS3_Lambda_0.75[[optimal_K]]$crossvalid.crossentropy) #0.3626061
min(TESS3_Lambda_0.5[[optimal_K]]$crossvalid.crossentropy) #0.3619763
min(TESS3_Lambda_1[[optimal_K]]$crossvalid.crossentropy) #0.3598105
min(TESS3_Lambda_1.25[[optimal_K]]$crossvalid.crossentropy) #0.3625805
min(TESS3_Lambda_1.5[[optimal_K]]$crossvalid.crossentropy) #0.3615824

#Best = Lambda_1

###7.3. FST OUTLIER DETECTION
#A. Use the best lambda project to retrieve results
res = Gettess3res(TESS3_Lambda_1, K=5)

#B. Compute the FST statistics using best run
FST = res$Fst 

#C. Compute the GIF - genomic inflation factor
lambda = res$gif
lambda # 3.610677

#D. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

#E. Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/2, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
hist(adj.p.values, col = "green")
C_fst <- candidates(alpha=0.1, adj.p.values)
ManPlot(adj.p.values, C_fst,"Fst")

#F. after you choose one, save the result
pdf("./Results_Filters/TESS_FST/TESS_FST_Outliers_adj_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()

#G. Candidate loci for FDR control: Benjamini-Hochberg at level q 10%
C_fst <- candidates(alpha=0.1, adj.p.values)

#H. save the Manhatan plot
pdf("./Results_Filters/TESS_FST/TESS_PlotK_Manhattan_pvalues.pdf", onefile = T) #add the number of snps removed in the name of pdf
ManPlot(adj.p.values, C_fst, paste0("Fst - Removing ",  length(C_fst), " SNPs"))
dev.off()

#I. Exclude candidate FST outlier
snps_fil_ldF_candidate <- Subset(snps, sites=C_fst)
snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_neutral <- Subset(snps, sites=N_snp)
VCFsummary(snps_neutral) #277 individuals and 4936 SNPs.
length(N_snp)
length(snps@site_id)-length(C_fst) 
length(snps_neutral@site_id)

neutral_after_fst_TESS = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_TESS #277 individuals and 4936 SNPs.

#J. Save the vcf
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_TESS.vcf"))


#------------------------------------------------------------------------------
#                       8. Saving Filtering Results 
#------------------------------------------------------------------------------
###8.1. SAVE THE FILTERING RESULTS BY DATASETS:
#A. Combine all results
datasets = c(raw_SZ, filtered_SZ, removed_missing_SZ, adaptativeFilter_SZ, adaptativeLD_SZ, neutralFilter_SZ, neutralLDWithin_SZ, neutralLDBetween_SZ, neutral_after_fst_PCadapt, neutral_after_fst_sNMF, neutral_after_fst_TESS)

#B. Create a matrix to save the results.
results_snps = matrix("NA", 11, 2)
colnames(results_snps) = c("Individuals", "SNPs")
rownames(results_snps) = c("Raw", "Biallelic SNPs", "Low Missing Data per Individual", "Quality Filter for Adaptive SNPs", "Linkage Desequilibrium Filter for Adaptive SNPs (Within Contigs)", "Quality Filter for Neutral SNPs", "Linkage Desequilibrium Filter for Neutral SNPs (Within Contigs)", "Linkage Desequilibrium Filter for Neutral SNPs (Between Contigs)", "Neutral SNPs (PCadapt FST outliers)", "Neutral SNPs (sNMF FST outliers)", "Neutral SNPs (TESS FST outliers)" )

#C. Save the results in matrix
for (i in 1:length(datasets)){
  text = scan(text = datasets[i], what = "")
  results_snps[i,] = c(text[1], text[4])
}

#D. Verify the results
as.data.frame(results_snps)

#E. Save result as .csv
write.csv(as.data.frame(results_snps), paste0("./Results_Metafiles/Results_SNPs_datasets_", project_name, ".csv"))

###8.2. CHOOSE BEST WAY TO FILTER LOCI OUTLIERS:
#A. Tess3 in my example:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_TESS.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#B. Save as final VCF:
Save(snps, paste0("vcf/", project_name, "_filtered_neutral.vcf"))


#------------------------------------------------------------------------------
#                       9. Creating RDA and LFMM2 inputs 
#------------------------------------------------------------------------------
#Some analysis as RDA and LFMM2 does not allow missing data, so we used LFMM to
#impute genetic missing data, based on the population that each individual
#belongs, using the package LEA
#------------------------------------------------------------------------------

###9.1. LOAD THE .VCF FILE WITH ADAPTATIVE SNPS:
snps <-  vcfLink(paste0("vcf/", project_name,"_filtered_adaptative.vcf"), overwriteID=T)
snps@meta
VCFsummary(snps) #277 individuals and 19025 SNPs.


###9.2. REPLACE THE MISSING VALUE BY "NA"
gen <- GenotypeMatrix(snps)
gen[gen==-1] <- NA
sum(is.na(gen))
dim(gen)


###9.3. CREATE THE LFMM FILE:
#A. Run sNMF to find clusters
snps_fil_lfmm <- Lfmm(snps, paste0("vcf/", project_name, "_filtered_adaptative.lfmm"))
project.snmf = snmf( paste0("vcf/", project_name, "_filtered_adaptative.lfmm"), K = 1, 
                    entropy = TRUE, repetitions = 10,
                    project = "new")
  
#B. Select the run with the lowest cross-entropy value
best = which.min(LEA::cross.entropy(project.snmf, K = 1))

#C. Impute the missing genotypes by clusters
impute(project.snmf,  paste0("vcf/", project_name, "_filtered_adaptative.lfmm"), method = 'mode', K = 1, run = best)

#D. Convert and save lfmm to geno and save
lfmm2geno(paste0("vcf/", project_name, "_filtered_adaptative.lfmm_imputed.lfmm"), output.file = paste0("vcf/", project_name,"_filtered_adaptative_imputed.geno"), force = TRUE)

##END
