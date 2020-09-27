###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
###########   STEP 04: GENETIC DIVERSITY, DISTANCE, AND TAJIMA D   ############

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
#A. CALCULATE GENETIC DIVERSITY INTRA AND INTER POPULATIONS/CLUSTERS
#A. CALCULATE GENETIC DISTANCE AMONG POPULATIONS/CLUSTERS AND INDIVIDUALS
#B. CALCULATE TAJIMA D FOR POPULATION EXPANSION


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
if("vcfR" %in% rownames(installed.packages()) == FALSE){install.packages("vcfR")
} else {print (paste0("'vcfR' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("reshape2" %in% rownames(installed.packages()) == FALSE){install.packages("reshape2")
} else {print (paste0("'reshape2' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("mmod" %in% rownames(installed.packages()) == FALSE){install.packages("mmod")
} else {print (paste0("'mmod' has already been installed in library"))}
if("poppr" %in% rownames(installed.packages()) == FALSE){install.packages("poppr")
} else {print (paste0("'poppr' has already been installed in library"))}
if("ecodist" %in% rownames(installed.packages()) == FALSE){install.packages("ecodist")
} else {print (paste0("'ecodist' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vegan, ecodist, vcfR, adegenet, poppr, mmod, reshape2, ggplot2, dartR)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Diversity", "./Results_Distance", "./Results_TajimaD"))


#------------------------------------------------------------------------------
#                        1. Loading Files
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "pilocarpus"

#B. Choose the number of clusters according to step 2. In my case TESS3. 

###1.2. LOAD VCF FILES AND GEOGRAPHICAL INFORMATIONS: 
#A. Load neutral .vcf file with geographical information and genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 6217 SNPs.
names(snps_neutral@meta) #verify col names in metafile

#B. Number of cluster and method:
optimal_K = 4
method = "TESS"

#B. Position of samples by populion by DAPC approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_tess))){
  pop = which(snps_neutral@meta$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}


###1.3. VERIFY THE MEAN COVERAGE DEPTH: 
#A. Mean coverage in all dataset
site.depth = Query(snps_neutral, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 83.81 / Median = 71.49
hist(site.depth$MEAN_DEPTH, breaks=30)

#B. Mean coverage by individual
coverage_ind = c()
for (p in 1:length(snps_neutral@sample_id)){
beta = Query(Subset(snps_neutral, samples = p), type="site-mean-depth")
coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps_neutral@meta$coverage = coverage_ind

write.csv(coverage_ind, file=paste0("Results_Metafiles/coverage_", project_name, "_by_individual.csv"))




#------------------------------------------------------------------------------
#                            2. Genetic Diversity
#------------------------------------------------------------------------------
###2.1. ESTIMATE GENETIC DIVERSITY:
#A. By species/taxon (all samples together)
Overall = GenDiv(snps_neutral)
write.csv(Overall, file=paste0("Results_Diversity/Diversity_Overall_neutral_", project_name, ".csv"))

#B. By individual
ind = GenDiv_IND(snps_neutral)
write.csv(ind, file=paste0("Results_Diversity/Diversity_Individual_neutral_", project_name , ".csv"))

#C. By cluster. ONLY IF YOUR K >= 2
#define a list of positions for each cluster
clusters = list(pop_TESS_1, pop_TESS_2, pop_TESS_3, pop_TESS_4)
#subset
for (i in 1:optimal_K){
  UNIND = snps_neutral@sample_id[clusters[[i]]]
  pop = Subset(snps_neutral, samples=UNIND)
  assign(paste0("pop_", i), pop)
}
#genetic diversity by cluster
pops = list(pop_1, pop_2, pop_3, pop_4)

for (j in 1:optimal_K){
cluster_div = GenDiv(pops[[j]])
write.csv(cluster_div, file=paste0("Results_Diversity/Diversity_Cluster_neutral_", project_name, "_POP", j, "_", method, ".csv"))
}

####2.2. COMPARE GENETIC DIVERSITY BETWEEN POPULATION OR SPECIES:
#A. Estimated diversity by individual in each population
ind_ste1 = GenDiv_IND(pop_1)
ind_ste1$SPE = c("STE_1")
ind_ste2 = GenDiv_IND(pop_2)
ind_ste2$SPE = c("STE_2")
ind_ste3 = GenDiv_IND(pop_3)
ind_ste3$SPE = c("STE_3")
ind_ste4 = GenDiv_IND(pop_4)
ind_ste3$SPE = c("STE_4")

#B.Create a data frame to analysis
table_diversity = rbind(ind_ste1,ind_ste2,ind_ste3,,ind_ste4)
variaveis = table_diversity[, 1:4]
species = as.factor(table_diversity[ ,5])
plan_anova = data.frame(species,variaveis)
classes = as.factor(plan_anova$species)

#C. Calculate ANOVA
sink(paste0("Results_Diversity/Ind_Diversity_By_Population_ANOVA_TUKEYS_", project_name, ".doc"))
for(i in 2:ncol(plan_anova)){
  column<- names (plan_anova[i])
  result_anova<- aov(plan_anova[,i]~classes, data= plan_anova)
  result_anova2<- summary(aov(plan_anova[,i]~classes, data= plan_anova))
  tk<- TukeyHSD (result_anova)
  print(column)
  print(result_anova2)
  print(tk)
}
sink()

#D. Calculate general p-value
sink(paste0("Results_Diversity/Ind_Diversity_By_Population_ANOVA_p-values_", project_name, ".doc"))
for(i in 2:ncol(plan_anova)){
  column<- names (plan_anova[i])
  result_anova<- summary(aov(plan_anova[,i]~classes, data= plan_anova))[[1]][["Pr(>F)"]]
  print(column)
  print(result_anova) 
}
sink()


#### 2.3. BOXPLOT TO COMPARE VALUES AMONG SPECIES:
#A. Reshape the dataframe:
dat.m = melt(plan_anova,id.vars="species")

#B. Create a graph
boxplot_diversity = ggplot(dat.m)+
  geom_boxplot(aes(x=species, y=value))+
  facet_wrap(~variable, scales="free")+
  geom_point(aes(species, value), size = 1, pch=19)+
  theme_bw()+
  labs(y= NULL, x = NULL)

#C.Verify
boxplot_diversity

#D. Save as pdf:
pdf(paste0("Results_Diversity/Diversity_Genetic_Boxplot_Populations_", project_name, ".csv"))
boxplot_diversity
dev.off()



#------------------------------------------------------------------------------
#                         3. Genotype Imputation
#------------------------------------------------------------------------------
###3.1 GENOTYPE IMPUTATION USING LEA. The genotypic matrix completion is based on estimated ancestry coefficients and ancestral genotype frequencies
#A. Convert VCF to LFMM
lfmm = Lfmm(snps_neutral, output.file=paste0("Results_Distance/Impute_LEA_", project_name, ".lfmm"))

#B. Run sNMF to calculate ancestry coefficients and ancestral genotype frequencies
project.snmf = snmf(paste0("Results_Distance/Impute_LEA_", project_name, ".lfmm"), K = optimal_K, entropy = TRUE, repetitions = 10, project = "new", seed=123)

#C. You may load project file if it was saved in previous run
#project.snmf = load.snmfProject(paste0("Results_Distance/Impute_LEA_", project_name, ".snmfProject"))

#D. select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.snmf, K = optimal_K))

#E. Impute the missing genotypes
impute(project.snmf, paste0("Results_Distance/Impute_LEA_", project_name, ".lfmm"), method = 'mode', K = optimal_K, run = best)

#F. Load imputed file
lfmmformat_imp = read.lfmm(paste0("Results_Distance/Impute_LEA_", project_name, ".lfmm_imputed.lfmm"))
lfmmformat_imp[1:10,1:10]
dim(lfmmformat_imp)

#G. Rename row and columns
length(snps_neutral@sample_id)
length(snps_neutral@site_id)
rownames(lfmmformat_imp) = snps_neutral@sample_id
colnames(lfmmformat_imp) = snps_neutral@site_id
lfmmformat_imp[1:10,1:10]


#------------------------------------------------------------------------------
#                     4. Genetic Distance by Individual
#------------------------------------------------------------------------------
###4.1 PCA-DISTANCE BASED ON BROKEN STICK RULE
#PCA-based distance and Dps: distance values among all pairs are calculated twice (below AND above diagonal) without distance within individuals
#A. Calcule PCA
PCAprcomp_imp = prcomp(lfmmformat_imp, center=TRUE, scale=TRUE)
str(PCAprcomp_imp)
summary(PCAprcomp_imp)

#B. Broken Stick PC values
screeplot(PCAprcomp_imp, bstick=TRUE, type="lines")
screeplot(PCAprcomp_imp, bstick=TRUE, type="barplot")
#4 PCs
summary(PCAprcomp_imp)
#4PCs are 10.662% of variance
n_pcs = 4

#C. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp = distance(PCAprcomp_imp$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp)
class(PC_distLEA_prcomp)
#convert to matrix
t_prcompLEA = as.matrix(PC_distLEA_prcomp)
t_prcompLEA[1:10,1:10]

#D. Save results
write.csv(t_prcompLEA, file=paste0("Results_Distance/PCA_Distance_BSR_IND_neutral_", project_name, ".csv"))


###4.2. PCA-DISTANCE BASED ON 95% OF VARIANCE
#A. Verify the number of PC that explain 95% of variance
summary(PCAprcomp_imp)
#248PCs are 94.893% of variance or see Step 2 in DAPC action #E
n_pcs = 248

##B. Calculate PCA-based distance based on 95% of variance
PC_distLEA_95 = distance(PCAprcomp_imp$x[,1:n_pcs], method = "euclidean")
head(PC_distLEA_95)
class(PC_distLEA_95)
t_prcomp95LEA = as.matrix(PC_distLEA_95)
t_prcomp95LEA[1:10,1:10]

#C. Save results
write.csv(t_prcomp95LEA, file=paste0("Results_Distance/PCA_Distance_95Var_IND_neutral_", project_name, ".csv"))


###4.3. RELATEDNESS IN ALL SAMPLES
#Yang's Relatedness: distance values among all pairs are calculated once (below diagonal) with distance within individuals
#A. Calculate Relatedness
REL_YANG = Relatedness(snps_neutral, type = "yang",verbose = TRUE)
head(REL_YANG)
nrow(REL_YANG) 
colnames(REL_YANG)<-c("INDV1","INDV2","RELATEDNESS_AJK_Yang")

#B. Save results:
write.csv(REL_YANG, file=paste0("Results_Distance/Yang_Reletedness_IND_neutral_", project_name, ".csv"))


###4.4. Dps DISTANCE
#A. Convert VCF to genind
snps =  read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), verbose = T)

#B. Convert data into genind
geninddata = vcfR2genind(snps)
geninddata@tab[1:10,1:10]

#C. Calculate proportions of shared alleles between pairs of individuals
dps_shared = propShared(geninddata)
head(dps_shared)

#D. Save results:
write.csv(dps_shared, file=paste0("Results_Distance/Dps_Shared_IND_neutral_", project_name, ".csv"))


#------------------------------------------------------------------------------
#                     5. Genetic Distance by Clusters
#------------------------------------------------------------------------------
###5.1. RELATEDNESS BY CLUSTERS
#A. Clusters were defined in #2.C
pops

#B. Calclate relatedness by cluster
for (j in 1:optimal_K){
  REL_POP = Relatedness(pops[[j]], type = "yang",verbose = TRUE)
  colnames(REL_POP) = c("INDV1","INDV2","RELATEDNESS_AJK_Yang")
  write.csv(REL_POP, file=paste0("Results_Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP", j, ".csv"))
}


###5.2. FST. ONLY IF YOUR K >= 2:
#A. Convert vcfR to genlight:
genlightdata = vcfR2genlight (snps)

#B. Define pop as POP_ID from DAPC
genlightdata@pop = as.factor(snps_neutral@meta$PopID_tess)

#C. Run FST
fst = gl.fst.pop(genlightdata, nboots = 100, percent = 95, nclusters = 1)
#save results
write.csv(fst$Fsts, file=paste0("Results_Distance/FST_clusters_Fstas_", project_name , "_", method, ".csv"))
write.csv(fst$Pvalues, file=paste0("Results_Distance/FST_clusters_pvalues_", project_name , "_", method, ".csv"))
write.csv(fst$Bootstraps, file=paste0("Results_Distance/FST_clusters_bootstrap_", project_name , "_", method, ".csv"))


#5.3. GST, D, AND OTHERS
#A.Define pop as POP_ID from TESS
geninddata
pop(geninddata) = snps_neutral@meta$PopID_tess
geninddata

#B. Define Inputs
datasets = c(geninddata)
names = c("TESS")

#C. Calculate summary statistic for dataset
summary_stat (datasets, names)


#------------------------------------------------------------------------------
#                              6. Tajima D
#------------------------------------------------------------------------------
###6.1. FILTER TO A SINGLE SNP PER CONTIG
#A. Load the VCF file:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.

#B. A single SNP per contig. This thins SNPs to a given distance in bp from one another. Setting the distance higher than the length of the contig ensures that you'll have a single SNP per contig.
snps_thin = Filter(snps_neutral, filterOptions(thin=300)) 
VCFsummary(snps_thin) #277 individuals and 3412 SNPs.


###6.2. SUBSET BY GENETIC CLUSTERS/POPULATIONS. SAME AS #3.C
#define a list of positions for each cluster
clusters = list(pop_TESS_1, pop_TESS_2, pop_TESS_3, pop_TESS_4)
#subset
for (i in 1:optimal_K){
  UNIND = snps_thin@sample_id[clusters[[i]]]
  pop = Subset(snps_thin, samples=UNIND)
  assign(paste0("popD_", i), pop)
}
#verify
VCFsummary(popD_1) #47 individuals and 3412 SNPs.
VCFsummary(popD_2) #47 individuals and 3412 SNPs.
VCFsummary(popD_3) #65 individuals and 3412 SNPs.
VCFsummary(popD_4) #118 individuals and 3412 SNPs.


###6.3. ESTIMATE TAJIMA'S D, BIAS-CORRECTED FOR MAF
tajd_p1 = TajimaD(popD_1, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_p1$results
str(tajd_p1$simulations)

tajd_p2 = TajimaD(popD_2, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_p2$results
str(tajd_p2$simulations)

tajd_p3 = TajimaD(popD_3, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_p3$results
str(tajd_p3$simulations)

tajd_p4 = TajimaD(popD_4, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_p4$results
str(tajd_p4$simulations)


###6.4. SAVE AND LOAD TAJIMA'S D RESULTS
#A. Save
save(tajd_p1, file = "./Results_TajimaD/tajd_p1.Rdata")
save(tajd_p2, file = "./Results_TajimaD/tajd_p2.Rdata")
save(tajd_p3, file = "./Results_TajimaD/tajd_p3.Rdata")
save(tajd_p4, file = "./Results_TajimaD/tajd_p4.Rdata")

#B. Load 
load("./Results_TajimaD/tajd_p1.Rdata")
load("./Results_TajimaD/tajd_p2.Rdata")
load("./Results_TajimaD/tajd_p3.Rdata")
load("./Results_TajimaD/tajd_p4.Rdata")


###6.5. PLOT OBSERVED TAJIMA'D AGAINST THE NULL DISTRIBUTION
#A.The null distribution (histogram) is shown next to the observed Tajima's D value (red line)

#B.POP1
pdf(paste0("./Results_TajimaD/TajimaD_POP1.pdf"), onefile = F)
ggplot(data.frame(x=tajd_p1$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p1$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw() + ggtitle("Simulations for POP1") +labs(y= "Frequency", x = "Tajima's D") 
dev.off()

#C.POP2
pdf(paste0("./Results_TajimaD/TajimaD_POP2.pdf"), onefile = F)
ggplot(data.frame(x=tajd_p2$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p2$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw() + ggtitle("Simulations for POP2") +labs(y= "Frequency", x = "Tajima's D") 
dev.off()

#D.POP3
pdf(paste0("./Results_TajimaD/TajimaD_POP3.pdf"), onefile = F)
ggplot(data.frame(x=tajd_p3$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p3$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw() + ggtitle("Simulations for POP3") +labs(y= "Frequency", x = "Tajima's D") 
dev.off()

#E.POP4
pdf(paste0("./Results_TajimaD/TajimaD_POP4.pdf"), onefile = F)
ggplot(data.frame(x=tajd_p4$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p4$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw() + ggtitle("Simulations for POP4") +labs(y= "Frequency", x = "Tajima's D") 
dev.off()


#------------------------------------------------------------------------------
#               7. Converting VCF to Genepop to run NeEstimator 2.14
#------------------------------------------------------------------------------
###5.1. LOAD VCF
snpsR = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), verbose = T)


###5.2. DEFINING POPULATIONS USING THE  VCF METAFILE:
snps = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps)
population = as.factor(snps@meta$PopID_tess)


##5.3. CONVERTING FILES:
#A. VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#B. Adding strata (pops) into Genind
snps_genind@pop = population

#C. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

#D. Converting Gtypes to GENEPOP to run in NeEstimator and save it:
genepopWrite(snps_gtypes, "pilocarpus_genepop.txt")

#E. Load this file in NeEstimator to run Ne analyses.

##END

##END
