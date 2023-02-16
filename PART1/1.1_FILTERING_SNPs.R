
####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
########################## STEP 01: FILTERING SNPS ############################




### Script prepared by Jeronymo Dalapicolla, Carolina S. Carvalho, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###




#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. ALIGNED AND CLEANED ".vcf" OR ".gvcf" FILE AFTER POPULATIONS STEP IN STACKS OR FROM ANOTHER SOFTWARE.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IT IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline

#C. THE FILE ".CSV" WITH THE GEOGRAPHICAL INFORMATION OF GENETIC SAMPLES IN DECIMALS DEGREES PREFERENCIALLY AND ONE COLUMN WITH LOCALITIES IDs. DIFFERENT INDIVIDUALS MUST HAVE DIFFERENT COORDINATES. IF TWO OR MORE SAMPLES WERE COLLECTED IN THE SAME LOCALITY, CHANGE THE VALUE OF 4TH OR 5TH DECIMAL NUMBER TO SET A DIFFERENCE OF 1 OR 5 METERS. SOME ANALYSES COULD RETURN ERRORS IF THE DIFFERENCES IN METERS BETWEEN INDIVIDUALS IS 0. 



  
##2. GOALS FOR THIS STEP ----
#A. FILTER RAW SNPS TWO DATASETS: ONE WITH NEUTRAL SNPS AND ANOTHER DATASET FOR LOCAL ADAPTATION ANALYSES.
  



##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

#B. CREATE A FOLDER NAMED "vcf" INSIDE YOUR WORKING DIRECTORY AND COPY THE .vcf OR .gvcf FILE THERE.

#C. CREATE A FOLDER NAMED "Inputs" INSIDE YOUR WORKING DIRECTORY AND COPY THERE, THE FILE .csv WITH GEOGRAPHICAL\TRAITS INFORMATION OF THE SAMPLES
  
  



##3. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
rm(list=ls())





##4. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP ----
#You can download it in https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
source("functions_LanGen.R")





##5. INSTALL AND LOAD THE PACKAGES ----
#For r2vcftools do you need install VCFTools in you computer:https://vcftools.github.io/index.html
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
if (!require('r2vcftools'))   remotes::install_github("nspope/r2vcftools");          library('r2vcftools')
if (!require('LEA'))          BiocManager::install("LEA");                           library('LEA')
if (!require('qvalue'))       BiocManager::install("qvalue");                        library('qvalue')
if (!require('tess3r'))       devtools::install_github("bcm-uga/TESS3_encho_sen");   library('tess3r')
if (!require('SNPRelate'))    BiocManager::install("SNPRelate");                     library('SNPRelate')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('poppr'))        install.packages("poppr");             library('poppr')
if (!require('pcadapt'))      install.packages("pcadapt");           library('pcadapt')
if (!require('psych'))        install.packages("psych");             library('psych')
if (!require('VennDiagram'))  install.packages("VennDiagram");       library('VennDiagram')
if (!require('sf'))           install.packages("sf");                library('sf')

#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, tidyverse, vcfR, SNPRelate, dartR, poppr, adegenet, LEA, pcadapt, qvalue, psych, tess3r, VennDiagram, sf)




##6. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS ----
create_dir(c("./Results/Step01/Linkage", "./Results/Step01/Metafiles", "./Results/Step01/PCA", "./Results/Step01/PCadapt", "./Results/Step01/sNMF", "./Results/Step01/TESS"))
  




#### ANALYSIS ---- 
  

#### 1. LOAD FILES -----
###1.1. CHOOSE A NAME FOR THE PROJECT AND THE PATH FOR VCF/GVCF AND COORDS FILES:
#A. vcf/gvcf file path:
vcf_file = "vcf/pilocarpusredo.vcf"
  
#B. Project name:
project_name = "pilocarpus"
  
#C. Coordenates file path and load it:
coord_file = "Inputs/pilocarpus_coords.csv"
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
aggregate = 5
  

#E. Load the VCF file
snps_raw =  vcfLink(vcf_file, overwriteID=T)
snps_raw

#F.Check the vcf file
snps_raw@sample_id
snps_raw@site_id
snps_raw@meta
Chrom(snps_raw) ##Chromosome, positions, and IDs

#G. Basic stats sample size (SZ)
raw_SZ = capture.output(VCFsummary(snps_raw)) 
raw_SZ #291 individuals and 36672 SNPs.




#### 2. BASIC FILTERS BY BIALLELIC SNPS AND BY INDIVIDUALS -----
#A. Remove loci that are not SNPs
snps = Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2), indels="remove")

# If you get any warning like this:
#"Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Ref+, Ref-..."
#It's just a warning that vcftools doesn't know how to handle the comma within the Description tag. If you remove that comma in the description, the warning will go away. Otherwise, it can generally be ignored. https://github.com/vcftools/vcftools/issues/134
#r2vcftools was designed for R3 so for R4 we have a issue with comma in Description tag and fucntion scan(), so if you have this problem, delete all ##INFO lines in this step. Keep a copy of the original vcf with INFO description tag for future check if you need!

#B. Basic stats sample size (SZ)
filtered_SZ = capture.output(VCFsummary(snps)) 
filtered_SZ #291 individuals and 35214 SNPs.


#C. Remove one individual
snps@meta$sample_name
UNIND = snps@sample_id[snps@sample_id != "9202_rep_sorted"]
snps_unind = Subset(snps, samples=UNIND)
snps_unind@meta
VCFsummary(snps_unind) ## 290 individuals and 35214 SNPs.


#D. Or remove more than one individual. In the case Granito locality, the rep individual (9202_rep), and individuals (9541, 9561):
snps@meta$sample_name
remove_ind = c("9202_rep_sorted", "9541_sorted", "9561_sorted", "9521_sorted","9522_sorted", "9523_sorted", "9526_sorted","9530_sorted","9531_sorted", "9532_sorted","9533_sorted","9537_sorted","9538_sorted","9540_sorted")
UNIND = snps@sample_id[!snps@sample_id %in% remove_ind]
snps_unind = Subset(snps, samples=UNIND)


#E. Verify sample size
length(snps@meta$sample_name)-length(remove_ind)
length(snps_unind@meta$sample_name)
VCFsummary(snps_unind) ## 277 individuals and 35214 SNPs.


#F. If you do not need to remove any individuals:
snps_unind = snps
VCFsummary(snps_unind) ## 277 individuals and 35214 SNPs.





#### 3. VERIFY DATA QUALITY -----
#A. Convert to a genotype matrix
genotypes = GenotypeMatrix(snps_unind)
genotypes[1:10, 1:10] ## -1 is missing; otherwise, gives the number of derived alleles in the genotype -- ie. a 0 is homozygotes for reference allele; 1 heterozygotes; and 2 is homozygotes for alternative allele. 
  

#B. Look at depth, quality, HWE, HE, allele frequencies, and Pi
site.depth = Query(snps_unind, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 63.40 / Median = 47.30
hist(site.depth$MEAN_DEPTH, breaks=30) 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <100]) ### >20 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <250]) ### <200

quality = Query(snps_unind, type="site-quality") 
summary(quality$QUAL) #Mean = 46.47 / Median = 48.63
hist(quality$QUAL) ## >30
  
PI = Query(snps_unind, type="site-pi")
mean(PI$PI) ## Mean nucleotide divergence per-site ###0.241
hist(PI$PI)
  
HWE = Query(snps_unind, type="hardy")
summary(HWE$P_HWE) #Mean = 0.129 / Median = 0.0012
hist(HWE$P_HWE)
hist(HWE$P_HWE[HWE$P_HWE<0.0001])
  
HE = Query(snps_unind, type="het") #Calculates a measure of heterozygosity on a per-individual basis.
hist(HE$O.HOM) # Observed homozygosity
hist(HE$E.HOM) #Expected homozygosity
hist(HE$N_SITES) # N_SITES
hist(HE$F) ## Inbreeding coefficient

#Heterozygosity estimation is implemented in GenDiv function in "function_LanGen.R" file
#(N_Sites - O(HOM))/ N_Sites = O(HET)
#(N_Sites - E(HOM))/ N_Sites = E(HET)
  
freq = Query(snps_unind, type="freq2")
hist(freq$X.FREQ.) #reference allele
hist(freq$X.FREQ.1) #alternative allele
head(freq)


#C. Missing per locus
Missing = apply(GenotypeMatrix(snps_unind), 2, function(x) sum(x < 0)/length(x)*100)
summary(Missing) ## Max missing = 54.9%
hist(Missing)


#D. Missing per individual
Missing_ind = apply(GenotypeMatrix(snps_unind),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 50.3%
hist(Missing_ind)


  
  
### 4. ADD IN THE VCF METAFILE THE GEOGRAPHIC AND MISSING DATA INFORMATION ----
#A. Save missing data information as data frame
Missingind.df = Missing_ind %>%
  as.data.frame() %>%
  rownames_to_column(., var = "ID_original") %>%
  mutate (ID = str_remove(ID_original, pattern="_sorted")) %>%
  setNames(., c("ID_original", "Missing_data", "ID"))


#B. Filter coords according to remaining individuals
set_coords = coords %>%
  as.data.frame() %>%
  setNames(., c("Local", "ID", "Longitude", "Latitude", "Aggregate"))
set_coords = semi_join(set_coords, Missingind.df)
set_coords = inner_join(Missingind.df, set_coords)
head(set_coords)


#C. Convert coordinates in UTM if you need. My original data is WGS and the UTM data will be in SAD69 UTM zone 22S. Search in internet the EPSG code for that area. For example https://epsg.io/ #In the case is EPSG:32722.
set_coords = st_as_sf(set_coords, #data frame
         coords = c("Longitude", "Latitude"), # for point data
         remove = F, # don't remove these lat/lon cols from df
         crs = 4326) %>% # add original projection (this is WGS84)
  st_transform(crs =st_crs(32722)) %>% #choose a new EPSG code
  mutate(X_UTM = sf::st_coordinates(.)[,1], Y_UTM = sf::st_coordinates(.)[,2]) %>%
  as.data.frame() %>%
  setNames(., c("ID_original", "Missing_data", "ID", "Local_ID", "Longitude", "Latitude", "Aggregate", "Geometry", "X_UTM", "Y_UTM")) %>%
  select(., -Geometry)
head(set_coords)


#D. Check the order:
identical(as.character(set_coords$ID_original), as.character(snps_unind@meta$sample_name))
setdiff(as.character(set_coords$ID_original), as.character(snps_unind@meta$sample_name))
setdiff(as.character(snps_unind@meta$sample_name), as.character(set_coords$ID_original))


#E. Add the information in the metafiles. You can add more columns if you need, like ecological traits or other clustering classification
head(snps_unind@meta) #we have sample_num and sample_name
snps_unind@meta = inner_join(snps_unind@meta, set_coords, by = c("sample_name" = "ID_original"))
head(snps_unind@meta)

  
  
###5. REMOVE INDIVIDUALS WITH LARGE AMOUNT OF MISSING DATA ----
#A. Test individuals with large amounts of missing (set to max 70% missing). Try different thresholds!
value = 70

ind_remove = snps_unind@meta %>%
  filter(Missing_data > value) %>%
  group_by(Local_ID) %>%
  count() %>%
  as.data.frame()

#B. Create a data frame for graphical output
pie(ind_remove$n, ind_remove$Local_ID, main="Missing Data per Locality")

#C. After testing different thresholds, set the missing data threshold and remove the individuals:
indtokeep = snps_unind@meta %>%
  filter(Missing_data < value) %>%
  select(sample_name)

snps_ind_pc = r2vcftools::Subset(snps_unind, samples = indtokeep$sample_name)
removed_missing_SZ = capture.output(VCFsummary(snps_ind_pc)) 
removed_missing_SZ ## 277 individuals and 35214 SNPs.







  

#### 6. FILTERING SNPS FOR ADAPTATIVE ANALYSES ----
###6.1. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE:
#A. From dataset with low missing by individuals. If the VCF does not have Quality information remove minQ argument.
snps_fil_low = Filter(snps_ind_pc, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200)) 
VCFsummary(snps_fil_low) #277 individuals and 25241 SNPs.
  
#B. From dataset with all individuals. If the VCF does not have Quality information remove minQ argument.
snps_fil_high = Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200)) 
VCFsummary(snps_fil_high) ##277 individuals and 25241 SNPs.

#C. Choose one of the datasets (A ou B):
snps_fil = snps_fil_low

#D. After filters, metafile is not save, so put it again in the new vcf
snps_fil@meta = snps_unind@meta
row.names(snps_fil@meta) = snps_unind@meta$sample_name
head(snps_fil@meta)
adaptativeFilter_SZ = capture.output(VCFsummary(snps_fil))


###6.2. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:
#A. Define threshold R² value for linkage:
r2 = 0.4

#B. Select SNPs by R² value with contigs:
ld_within = Linkage(snps_fil, type="geno-r2", linkageOptions(min.r2= r2)) 
write.csv(ld_within, file=paste0("Results/Step01/Linkage/Filtering/ld_within_", r2 , "_test2.csv"))

# or load a saved result
ld_within = read.csv(paste0("Results/Step01/Linkage/ld_within_", r2 , "_test2.csv"))
head(ld_within)
hist(ld_within$R.2)

  
#C. Subset the VCF keeping one of the sets of the correlated snps (ID1 or ID2)
nold_snps = anti_join(as.data.frame(snps_fil@site_id), ld_within, by = c("snps_fil@site_id" = "ID1"))
head(nold_snps)

snps_fil_ld = Subset(snps_fil, sites = nold_snps$`snps_fil@site_id`)
adaptativeLD_SZ = capture.output(VCFsummary(snps_fil_ld))
adaptativeLD_SZ #277 individuals and 19025 SNPs.



###6.3. VERIFY THE REAL MISSING DATA AND COVERAGE DEPTH VALUES IN THE DATASET AFTER REMOVING CORRELATED SNPS:
#A. Mean coverage in all dataset
site.depth = Query(snps_fil_ld, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 58.61 / Median = 47.17
hist(site.depth$MEAN_DEPTH, breaks=30)

#B. Mean coverage by individual
coverage_ind = c()
for (p in 1:length(snps_fil_ld@sample_id)){
  beta = Query(Subset(snps_fil_ld, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
head(coverage_ind)

#C. Save results in the VCF Metafile
snps_fil_ld@meta$Coverage_depth = coverage_ind
head(snps_fil_ld@meta)

#D. Missing per individual:
Missing_ind = apply(GenotypeMatrix(snps_fil_ld),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 50.16%
hist(Missing_ind)
  
#E. Missing per locus:
Missing = apply(GenotypeMatrix(snps_fil_ld), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 29.96%
hist(Missing)

#F. Replace values of missing data after filters in the VCF Metafile
snps_fil_ld@meta$Missing_data = as.data.frame(Missing_ind)[,1]
head(snps_fil_ld@meta)

#G. Save metafile information:
write_csv(snps_fil_ld@meta, paste0("./Results/Step01/Metafiles/Metafile_adaptive_vcf_filtered_", project_name, ".csv"))
  
#H. Save the VCF
Save(snps_fil_ld, paste0("vcf/", project_name,"_filtered_adaptative.vcf"))




  
#### 7. FILTERING FOR NEUTRAL SNPS - PART 1 (QUALITY & LINKAGE) -----

#7.1. FILTER DATASET BY QUALITY, MISSING, ALLELIC FREQUENCY, MIN AND MAX COVERAGE, AND HARDY-WEINBERG EQUILIBRIUM (HWE):
#A. From dataset with low missing by individuals
snps_fil_hwe_low = Filter(snps_ind_pc, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
VCFsummary(snps_fil_hwe_low) #277 individuals and 13958 SNPs.
  
#B. From dataset with all individuals
snps_fil_hwe_high = Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.7, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
VCFsummary(snps_fil_hwe_high)  ## 277 individuals and 13958 SNPs.
  
#C. Choose one of the datasets (A ou B):
snps_fil_hwe = snps_fil_hwe_low
neutralFilter_SZ = capture.output(VCFsummary(snps_fil_hwe))
neutralFilter_SZ

#D. After filters, metafile is not save, so put it again in the new vcf
snps_fil_hwe@meta = snps_unind@meta
row.names(snps_fil_hwe@meta) = snps_unind@meta$sample_name
head(snps_fil_hwe@meta)




###7.2. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN AND BETWEEN CONTIGS:
## If you have more than one population, you can use this to identify SNPs deviating from HW equilibrium within each population, and then removes those SNPs that are in desequilibrium in all populations. You just need to subset your samples:

#A. Define R² value:
r2 = 0.4  


#B. Select SNPs by R² value within contigs:
ld_within = Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=r2))
write.csv(ld_within, file=paste0("Results/Step01/Linkage/ld_within_",r2, "_hwe_test2.csv"))

# or load a saved result
ld_within = read.csv(paste0("Results/Step01/Linkage/ld_within_",r2, "_hwe_test2.csv"))
head(ld_within)
hist(ld_within$R.2)


#C. Subset the VCF keeping one of the sets of the correlated snps (ID1 or ID2)
nold_snps = anti_join(as.data.frame(snps_fil_hwe@site_id), ld_within, by = c("snps_fil_hwe@site_id" = "ID1"))
head(nold_snps)

snps_fil_ldn = Subset(snps_fil_hwe, sites = nold_snps$`snps_fil_hwe@site_id`)
neutralLDWithin_SZ = capture.output(VCFsummary(snps_fil_ldn))
neutralLDWithin_SZ ##277 individuals and 11104 SNPs.


#D. Select SNPs by R² value between contigs:
ld_between = Linkage(snps_fil_ldn, type="interchrom-geno-r2", linkageOptions(min.r2=r2))
write.csv(ld_between, file= paste0("Results/Step01/Linkage/ld_between_", r2, "_hwe_test2.csv"))

#or load a saved results
ld_between = read.csv(paste0("Results/Step01/Linkage/ld_between_", r2, "_hwe_test2.csv"))
head(ld_between)
hist(ld_between$R.2)


#E. Subset the VCF keeping one of the sets of the correlated snps (ID1 or ID2)
nold2_snps = anti_join(as.data.frame(snps_fil_ldn@site_id), ld_between, by = c("snps_fil_ldn@site_id" = "ID1"))
head(nold2_snps)

snps_fil_ldF = Subset(snps_fil_ldn, sites = nold2_snps$`snps_fil_ldn@site_id`)
neutralLDBetween_SZ = capture.output(VCFsummary(snps_fil_ldF)) 
neutralLDBetween_SZ ##277 individuals and 6602 SNPs.
  




###7.3. VERIFY THE REAL MISSING DATA AND COVERAGE DEPTH VALUES IN THE DATASET AFTER REMOVING CORRELATED SNPS:
#A. Mean coverage in all dataset
site.depth = Query(snps_fil_ldF, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 58.78 / Median = 48.75
hist(site.depth$MEAN_DEPTH, breaks=30)

#B. Mean coverage by individual
coverage_ind = c()
for (p in 1:length(snps_fil_ldF@sample_id)){
  beta = Query(Subset(snps_fil_ldF, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
head(coverage_ind)

#C. Save results in the VCF Metafile
snps_fil_ldF@meta$Coverage_depth = coverage_ind
head(snps_fil_ldF@meta)

#D. Missing per individual:
Missing_ind = apply(GenotypeMatrix(snps_fil_ldF),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 47.73%
hist(Missing_ind)

#E. Missing per locus:
Missing = apply(GenotypeMatrix(snps_fil_ldF), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 29.96%
hist(Missing)

#F. Replace values of missing data after filters in the VCF Metafile
snps_fil_ldF@meta$Missing_data = as.data.frame(Missing_ind)[,1]
head(snps_fil_ldF@meta)
VCFsummary(snps_fil_ldF) ##277 individuals and 6602 SNPs.

#G. Save metafile information:
write_csv(snps_fil_ldF@meta, paste0("./Results/Step01/Metafiles/Metafile_neutral_partial_vcf_filtered_", project_name, ".csv"))

#H. Save the VCF
Save(snps_fil_ldF, paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"))
  







#### 8. TESTING THE MISSING DATA INFLUENCE ON GENETIC STRUCTURE ----
#A. Load neutral and adaptive .vcf files from filtering step:
vcf_adaptative = read.vcfR(paste0("vcf/", project_name,"_filtered_adaptative.vcf"), verbose = FALSE)
vcf_neutral = read.vcfR(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), verbose = FALSE)


#B. Convert "VCFs" to "GENIND"
genind_adap = vcfR2genind(vcf_adaptative)
genind_neut = vcfR2genind(vcf_neutral)


#C. Remove missing data in both datasets:
genind_adap_0 = missingno(genind_adap)
genind_adap_0 #verify number of SNPs = 2,921 SNPs

genind_neut_0 = missingno(genind_neut)
genind_neut_0 #verify number of SNPs =  927 SNPs


#D. Defining all datasets, names, and the vector with missing values. Define a threshold in % for missing data. Samples > threshold for missing data will be plot in red on the PCA.
datasets = c(genind_adap, genind_adap_0, genind_neut, genind_neut_0)
names = c("Adapative_all","Adapative_0", "Neutral_all", "Neutral_0")
threshold = 20
missing_values = snps_fil_ldF@meta$Missing_data #between adaptive and neutral dataset the amount of missing did not change so much. However you can choose the adaptive values for missing data, the individuals must be the same!

#E. Run PCA by dataset
pca_missing (datasets, names, threshold, missing_values)
  
#F. Verify in the PCA graphs if there are clusters in your samples, if the number of clusters changed with or without missing data and if any individual in red has changed its position among the clusters. Greater spreading in data with more missing data and SNPs is common. Be careful! Changes could be due to a drop in the number of SNPs. Data with 1,000 SNPs are more robust to keep the genetic structure pattern. If necessary, re-filter and remove unstable individuals or reduce the allowed missing data threshold in Action #1.6, #2.1, and #3.1.
 








#### 9. FILTERING FOR NEUTRAL SNPS - PART 2 (OUTLIER SNPS) -----

### 9.1. PCadapt ----
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

#A. Load the vcf file "neutral_partial":
snps = vcfLink(paste0("./vcf/", project_name,"_filtered_neutral_partial.vcf"), meta = NULL, overwriteID=F)
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
lfmm_input = Lfmm(snps, output.file=paste0("./Results/Step01/PCadapt/", project_name, "_filtered_snps_lfmm_format.lfmm"))
lfmm_input = read.lfmm(paste0("./Results/Step01/PCadapt/", project_name, "_filtered_snps_lfmm_format.lfmm"))
class(lfmm_input)
  
#E. Convert LFMM to a PCadapt matrix and run the analyses
pcadapt.test = lfmm_input %>%
  read.pcadapt(., type="lfmm") %>%
  pcadapt(., K=10, ploidy=2)
  
#F. "Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)
plot(pcadapt.test, option="screeplot") #2 PCs
plot(pcadapt.test, option="scores") #PC1 and PC2 #3 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3 # 4 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4 # 1-2 clsuters
plot(pcadapt.test, option = "scores", i = 5, j = 4) #PC5 and PC4 # 1-2 clsuters
plot(pcadapt.test, option = "scores", i = 5, j = 6) #PC5 and PC6 # 2 clusters
#With PC2 is reached the maximum number of clusters #4
#With PC3, PC4, PC5 the groups have changed to 1-2.
#So we can use 2 PCs, K = 2.
  
#Save figure as pdf files:
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_screenplot.pdf"), onefile = T)
plot(pcadapt.test, option="screeplot")
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_scores_PC1-PC2.pdf"), onefile = T)
plot(pcadapt.test, option="scores") #PC1 and PC2
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_scores_PC2-PC3.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_scores_PC3-PC4.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_scores_PC4-PC5.pdf"), onefile = T)
plot(pcadapt.test, option = "scores", i = 4, j = 5) #PC4 and PC5
dev.off()


#G. K = 2. Run PCadapt again. Test K and see the p-values graphs
pcadapt.test = pcadapt(gen.pcadapt, K=2, ploidy=2, min.maf=0.05, method="mahalanobis")
summary(pcadapt.test)
# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space. "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.
  
#H. Graphical tools.
plot(pcadapt.test, option = "manhattan")
plot(pcadapt.test, option = "qqplot")
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.test, option = "stat.distribution")
  
#save graphs
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_manhattan.pdf"), onefile = T)
plot(pcadapt.test, option = "manhattan")
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_qqplot.pdf"), onefile = T)
plot(pcadapt.test, option = "qqplot")
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_p-values.pdf"), onefile = T)
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()
  
pdf(paste0("./Results/Step01/PCadapt/", project_name, "_stat_distribution.pdf"), onefile = T)
plot(pcadapt.test, option = "stat.distribution")
dev.off()
  


#I. Choosing a cutoff for outlier detection
pcadapt.test$gif #  1.214431
  

# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.
#
# GIF of 1=well calibrated, >1 =liberal (too many small p-values), <1=conservative (too few small p-values) # Note: GIFs > 2 indicate poor test calibration.
#
#For a given alpha (real valued number between 0 and 1), SNPs with q-values less than alpha will be considered as outliers with an expected false discovery rate bounded by alpha. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs, for an expected false discovery rate lower than 10%.

  
#I1. q-values
qval = qvalue(pcadapt.test$pvalues)$qvalues
alpha = 0.1
outliers1 = which(qval < alpha)
length(outliers1) #It will be eliminated 152 SNPs
  
#I2. Benjamini-Hochberg Procedure
padj = p.adjust(pcadapt.test$pvalues, method="BH")
alpha = 0.1
outliers2 = which(padj < alpha)
length(outliers2) #It will be eliminated 152 SNPs
  
#I3. Bonferroni correction
padj2 = p.adjust(pcadapt.test$pvalues, method="bonferroni")
alpha = 0.1
outliers3 = which(padj2 < alpha)
length(outliers3) #It will be eliminated 20 SNPs
  

#J. Choose one approach to eliminate outlier SNPs. In the case Benjamini-Hochberg Procedure the "outliers2"
snps_neutral = Subset(snps, sites=snps@site_id[-outliers2])
VCFsummary(snps_neutral)
#verify if the exclusion was correct
length(snps@site_id)-length(outliers2) 
length(snps_neutral@site_id)
  
#K. Save neutral SNP dataset
neutral_after_fst_PCadapt = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_PCadapt #"277 individuals and 6450 SNPs."

write_csv(as.data.frame(snps@site_id[outliers2]), paste0("./Results/Step01/PCadapt/Metafile_names_outlier_SNPs_PCadapt_", project_name, ".csv"))
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_PCadapt.vcf"))




### 9.2. sNMF ----
# sNMF - sparse Nonnegative matrix factorization (NMF) for Estimation of Individual Ancestry Coefficients similar to STRUCTURE and ADMIXTURE.
# STRUCTURE and ADMIXTURE assumptions include absence of genetic drift and Hardy–Weinberg and linkage equilibrium in ancestral populations. sNMF was more appropriate to deal with inbred lineages. sNMF enabled the estimation of homozygote and heterozygote frequencies and avoided Hardy–Weinberg equilibrium assumptions. Cross-entropy criterion indicated better predictive results for sNMF than for ADMIXTURE. You can use impute() to predict Ancestry Coefficients for missing values.

#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.
  
#B. Convert VCF to geno object. You need to specify the output file. It will automatically subset the vcf file and assign it as a new object.
snps = Geno(snps, output.file = paste0("./Results/Step01/sNMF/", project_name, "_filtered_neutral_partial.geno"))
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
  path = paste0("./Results/Step01/sNMF/Alpha", i)
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("./Results/Step01/sNMF/", project_name, "_filtered_neutral_partial.geno"), path )
  }
  
#E. Set parameters to run SNMF (LEA) using different alpha values.
K = c(1:10) # set the number of K to be tested
replications = 5 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel
  

#F.Run a loop for all alpha values in RUN sNMF (LEA)
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results/Step01/sNMF/Alpha", i,"/", project_name, "_filtered_neutral_partial.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop), pro_snmf)
  }
  
#G. To load the SNMF projects. This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results/Step01/sNMF/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop), pro_snmf)
  }
  
#H. Summary of the project
summary(project1)
summary(project2)
summary(project3) #0.5497001
summary(project4)
summary(project5)
summary(project6)
  
#I. View cross-Entropy plot
PlotK(project1) #5
PlotK(project2) #5
PlotK(project3) #5
PlotK(project4) #4
PlotK(project5) #4
PlotK(project6) #5
#another way to view the results
#plot(project1, lwd = 5, col = "red", pch=1)
#plot(project2, lwd = 5, col = "red", pch=1)
#plot(project3, lwd = 5, col = "red", pch=1)
#plot(project4, lwd = 5, col = "red", pch=1)
#plot(project5, lwd = 5, col = "red", pch=1)
#plot(project6, lwd = 5, col = "red", pch=1)
  
#J. Save Cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results/Step01/sNMF/Cross_Entropy_sNMF_Alpha_",  i, ".pdf"), onefile = F)
  path = paste0("./Results/Step01/sNMF/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
  }
  
#K. Select optimal K value
optimal_K = 5
  
#L. ATTENTION, if your dataset is K = 1 force a K = 2 to be able to filter SNPs with FST outliers.
  
#M. Select best run (lowest cross-entropy)  
best_run = Best.run(nrep=replications, optimalK=optimal_K, p1=project1, p2=project2, p3=project3, p4=project4, p5=project5, p6=project6)
best_run
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("./Results/Step01/sNMF/Alpha", alpha_values[as.numeric(best_run_split[6])],"/", project_name, "_filtered_neutral_partial.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])
  
  
#N. Compute the FST and GIF - genomic inflation factor statistics using best run 
FST = fst(project, run, optimal_K) # you need at least 2 populations for a population-based test, so K>1.
lambda = GIF(project, run, optimal_K, fst.values=FST)
lambda #6.224445
  
#O. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")
  
pdf("./Results/Step01/sNMF/sNMF_FST_Outliers_raw_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()
  
#P. Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/3, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
hist(adj.p.values, col = "green")
C_fst = candidates(alpha=0.05, adj.p.values) #Candidate loci for FDR control: Benjamini-Hochberg at level q
ManPlot(adj.p.values, C_fst,"Fst")
  
#after you choose one, save the result
pdf("./Results/Step01/sNMF/sNMF_FST_Outliers_adj_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()

pdf("./Results/Step01/sNMF/sNMF_FST_Outliers_Manhatan_plot.pdf", onefile = T)
ManPlot(adj.p.values, C_fst, paste0("Fst - Removing ",  length(C_fst), " SNPs"))
dev.off()


#Q. Eliminate outlier SNPs
snps_neutral = Subset(snps, sites=snps@site_id[-C_fst])
VCFsummary(snps_neutral)
#verify if the exclusion was correct
length(snps@site_id)-length(C_fst) 
length(snps_neutral@site_id)

#K. Save neutral SNP dataset
neutral_after_fst_sNMF = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_sNMF #277 individuals and 4922 SNPs.

#I. Save the results
write_csv(as.data.frame(snps@site_id[C_fst]), paste0("./Results/Step01/sNMF/Metafile_names_outlier_SNPs_snmf_", project_name, ".csv"))
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_sNMF.vcf"))











### 9.3. TESS3 ----
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html


#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.


#B. Create a Genotype matrix
genotypes = GenotypeMatrix(snps) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing;
class(genotypes)
genotypes = replace(genotypes, genotypes == -1, NA)


#C. Create a Matrix with long and lat 
coordinates = snps@meta %>%
  select(Longitude, Latitude) %>%
  data.matrix(., rownames.force = NA)
#verify the coords
plot(coordinates, pch = 19, cex = .5, xlab = "Longitude", ylab = "Latitude")


#D. Costumize values and run a loop for all alpha values
lambda_values = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5) #test lambda values around 1.
K = c(1:10) # set the number of K to be tested
replications = 5 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel
mask = 0.05 #proportion of masked values

set.seed(13)
for (i in lambda_values){
    tess3.ls = tess3(genotypes, coord = coordinates, K = K, mask = mask, lambda = i,
                   method = "projected.ls", max.iteration = 5000, rep = replications,
                   ploidy = ploidy, openMP.core.num = CPU)
    save(tess3.ls, file = paste0("./Results/Step01/TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
  
}

 
#E. Choose the lambda with the minimum cross-validation value:
#create a matrix to save results
cross_value = matrix(NA, max(K)*length(lambda_values), 3)
colnames(cross_value) = c("K", "Lambda", "Crossvalidation")
  
loop = 0 #Always set loop = 0.
for (i in lambda_values){
  load(file = paste0("./Results/Step01/TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
    for (j in 1:max(K)){
      loop=loop+1
      res = Gettess3res(tess3.ls, K=j)
      cross_value[loop,1] = j
      cross_value[loop,2] = i
      cross_value[loop,3] = min(res$crossvalid.crossentropy)}
  }

#save as csv
write.csv(cross_value,paste0("./Results/Step01/TESS/Tess3ls_Crossvalidation_values_", project_name, ".csv"))
#choose best lambda
lambda_tess = as.vector(cross_value[cross_value[,3] == min(cross_value[,3]), ][2])
lambda_tess #0.5


#F. Choose best K:
#load the best lambda project:
load(file = paste0("./Results/Step01/TESS/Tess3ls_Lambda_", lambda_tess,"_", project_name, ".RData"))
#plot results
pdf(paste0("./Results/Step01/TESS/TESS3_PlotK_CrossValidation_Lambda_", lambda_tess, ".pdf"), onefile =F)
plot.new()
plot(tess3.ls, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")
dev.off()
#best K
optimal_K = 5
#ATTENTION, if your dataset is K = 1 force a K = 2 to be able to filter SNPs with FST outliers.
  
  
#G. Compute results for the best K
res = Gettess3res(tess3.ls, K=optimal_K)
FST = res$Fst #Select FST statistics for the best run
lambda = res$gif #Compute the GIF - genomic inflation factor
lambda # 3.750238

  
#H. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

  
#I. Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/2, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
hist(adj.p.values, col = "green")
C_fst = candidates(alpha=0.05, adj.p.values)
ManPlot(adj.p.values, C_fst,"Fst")

#after you choose one, save the result
pdf("./Results/Step01/TESS/TESS_FST_Outliers_adj_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()
  
#save the Manhatan plot
pdf("./Results/Step01/TESS/TESS_PlotK_Manhattan_pvalues.pdf", onefile = T) #add the number of snps removed in the name of pdf
ManPlot(adj.p.values, C_fst, paste0("Fst - Removing ",  length(C_fst), " SNPs"))
dev.off()


#J. Eliminate outlier SNPs
snps_neutral = Subset(snps, sites=snps@site_id[-C_fst])
VCFsummary(snps_neutral)
#verify if the exclusion was correct
length(snps@site_id)-length(C_fst) 
length(snps_neutral@site_id)


#K. Save the vcf files
neutral_after_fst_TESS = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_TESS #277 individuals and 5268 SNPs.

write_csv(as.data.frame(snps@site_id[C_fst]), paste0("./Results/Step01/TESS/Metafile_names_outlier_SNPs_TESS_", project_name, ".csv"))
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_TESS.vcf"))
  





  
  
#### 10. SAVING FILTERING RESULTS ---- 

#A. Number of SNPs and Individuals in the filtering steps
datasets = c(raw_SZ, filtered_SZ, removed_missing_SZ, adaptativeFilter_SZ, adaptativeLD_SZ, neutralFilter_SZ, neutralLDWithin_SZ, neutralLDBetween_SZ, neutral_after_fst_PCadapt, neutral_after_fst_sNMF, neutral_after_fst_TESS)

#create a matrix to save the results.
results_snps = matrix("NA", 11, 2)
colnames(results_snps) = c("Individuals", "SNPs")
rownames(results_snps) = c("Raw", "Biallelic SNPs", "Low Missing Data per Individual", "Quality Filter for Adaptive SNPs", "Linkage Desequilibrium Filter for Adaptive SNPs (Within Contigs)", "Quality Filter for Neutral SNPs", "Linkage Desequilibrium Filter for Neutral SNPs (Within Contigs)", "Linkage Desequilibrium Filter for Neutral SNPs (Between Contigs)", "Neutral SNPs (PCadapt FST outliers)", "Neutral SNPs (sNMF FST outliers)", "Neutral SNPs (TESS FST outliers)" )
  
#Save the results in matrix
for (i in 1:length(datasets)){
  text = scan(text = datasets[i], what = "")
  results_snps[i,] = c(text[1], text[4])
}
as.data.frame(results_snps)
  
#Save result as .csv
write.csv(as.data.frame(results_snps), paste0("./Results/Step01/Metafiles/Results_SNPs_datasets_", project_name, ".csv"))




#B. Creating a Venn Diagram to verify shared SNPs in different filters approach and save it as pdf
A_pcadapt = read.csv("./Results/Step01/PCadapt/Metafile_names_outlier_SNPs_PCadapt_pilocarpus.csv")
A_pcadapt = as.character(A_pcadapt$snps.site_id.outliers2.)
A_pcadapt

B_snmf = read.csv("./Results/Step01/sNMF/Metafile_names_outlier_SNPs_snmf_pilocarpus.csv")
B_snmf = as.character(B_snmf$snps.site_id.C_fst.)
B_snmf

C_TESS = read.csv("./Results/Step01/TESS/Metafile_names_outlier_SNPs_TESS_pilocarpus.csv")
C_TESS = as.character(C_TESS$snps.site_id.C_fst.)
C_TESS

pdf(paste0("./Results/Step01/Metafiles/VennDiagram_Outliers_SNPs.pdf"), onefile = T)
plot.new()
grid.draw(w1 <- venn.diagram(list(PCadapt=A_pcadapt, sNMF=B_snmf, TESS=C_TESS),
                             lty = c("blank", "blank", "blank"),
                             fill = c("red", "blue", "yellow"),
                             alpha = c(0.5, 0.5, 0.5), cat.cex = 1.2, cex= 1.5,  cat.pos = 0, 
                             filename=NULL ))

dev.off()



#C. Choose best approach to filter outlier SNPs. I chose TESS3 in my example:
snps_final = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_TESS.vcf"), overwriteID=T)
VCFsummary(snps_final) ##277 individuals and 5268 SNPs.


#D. Verify real missing data and coverage depth in the final dataset:
#A. Mean coverage in all dataset
site.depth = Query(snps_final, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 57.34 / Median = 47.84
hist(site.depth$MEAN_DEPTH, breaks=30)

#B. Mean coverage by individual
coverage_ind = c()
for (p in 1:length(snps_final@sample_id)){
  beta = Query(Subset(snps_final, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
head(coverage_ind)

#C. Save results in the VCF Metafile
snps_final@meta$Coverage_depth = coverage_ind
head(snps_final@meta)

#D. Missing per individual:
Missing_ind = apply(GenotypeMatrix(snps_final),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 47.63%
hist(Missing_ind)

#E. Missing per locus:
Missing = apply(GenotypeMatrix(snps_final), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 29.96%
hist(Missing)

#F. Replace values of missing data after filters in the VCF Metafile
snps_final@meta$Missing_data = as.data.frame(Missing_ind)[,1]
head(snps_final@meta)

#G. Save metafile information:
write_csv(snps_final@meta, paste0("./Results/Step01/Metafiles/Metafile_neutral_final_vcf_filtered_", project_name, ".csv"))

#H. Save the VCF
Save(snps_final, paste0("vcf/", project_name,"_filtered_neutral.vcf"))

##END
