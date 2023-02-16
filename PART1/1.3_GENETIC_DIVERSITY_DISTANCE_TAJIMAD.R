
####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
######### STEP 03: GENETIC DIVERSITY, DISTANCE, TAJIMA D, AND NE  #############




### Script prepared by Jeronymo Dalapicolla, Carolina S. Carvalho, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###




#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. THE FILE ".VCF" CLEANED AFTER FILTERING AND WITH DELIMITED GENETIC CLUSTERS, STEP 2.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline



##2. GOALS FOR THIS STEP:
#A. CALCULATE GENETIC DIVERSITY INTRA AND INTER POPULATIONS/CLUSTERS
#B. CALCULATE GENETIC DISTANCE AMONG POPULATIONS/CLUSTERS AND INDIVIDUALS
#C. CALCULATE TAJIMA D FOR POPULATION EXPANSION
#D. CONVERT INPUT FOR NEESTIMATOR 2.1



##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H


##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())




##5. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("functions_LanGen.R")


##6. INSTALL AND LOAD THE PACKAGES ----
#For r2vcftools do you need install VCFTools in you computer:https://vcftools.github.io/index.html
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
if (!require('r2vcftools'))   remotes::install_github("nspope/r2vcftools");          library('r2vcftools')
if (!require('LEA'))          BiocManager::install("LEA");                           library('LEA')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('ggplot2'))      install.packages("ggplot2");           library('ggplot2')
if (!require('mmod'))         install.packages("mmod");              library('mmod')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2')
if (!require('poppr'))        install.packages("poppr");             library('poppr')
if (!require('ecodist'))      install.packages("ecodist");           library('ecodist')
if (!require('vegan'))        install.packages("vegan");             library('vegan')
if (!require('gghighlight'))  install.packages("gghighlight");       library('gghighlight')
if (!require('ggpubr'))       install.packages("ggpubr");            library('ggpubr')
if (!require('usedist'))      install.packages("usedist");           library('usedist')
if (!require('strataG'))      install.packages("strataG");           library('strataG')

#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, LEA, vegan, ecodist, vcfR, adegenet, poppr, mmod, reshape2, ggplot2, dartR, gghighlight, ggpubr, usedist, strataG)



##7. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS ----
create_dir(c("./Results/Step03/Diversity", "./Results/Step03/Distance", "./Results/Step03/TajimaD", "./Results/Step03/Ne" ))



##8. CREATE A PATTERN FOR GRAPHIC FOLDERS S TO SAVE THE RESULTS ----
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
#A. Project name:
project_name = "pilocarpus"

#B. Load neutral .vcf file with geographical information and genetic clusters ID, choosen in step 2. 
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.
names(snps_neutral@meta) #verify col names in metafile

#C. Number and name of cluster and method:
optimal_K = 4
name_clusters = c("A", "B", "C", "D")

#D.Position of samples by population/genetic clusters:
for (i in name_clusters){
  pop = which(snps_neutral@meta$POP_ID == i)
  assign(paste0("pop_POSI_", i), pop)
}








#### 2. GENETIC DIVERSITY BY CLUSTERS ----
#A. By species/taxon (all samples together)
Overall = GenDiv(snps_neutral)
write.csv(Overall, file=paste0("Results/Step03/Diversity/Diversity_Overall_", project_name, ".csv"))


#B. By cluster. ONLY IF YOUR K >= 2
#define a list of positions for each cluster
clusters = list(pop_POSI_A, pop_POSI_B, pop_POSI_C, pop_POSI_D)
#subset
for (i in 1:optimal_K){
  UNIND = snps_neutral@sample_id[clusters[[i]]]
  pop = Subset(snps_neutral, samples=UNIND)
  assign(paste0("pop_", name_clusters[i]), pop)
}
#genetic diversity by cluster
pops = list(pop_A, pop_B, pop_C, pop_D)

for (j in 1:optimal_K){
  cluster_div = GenDiv(pops[[j]])
  write.csv(cluster_div, file=paste0("Results/Step03/Diversity/Diversity_Clusters_", project_name, "_POP_", name_clusters[j], ".csv"))
}






#### 3. GENETIC DIVERSITY BY INDIVIDUALS ----
#A. All individuals
ind = GenDiv_IND(snps_neutral)
write.csv(ind, file=paste0("Results/Step03/Diversity/Diversity_Individual_", project_name , ".csv"))



#B. By Individual in each cluster 
pops = list(pop_A, pop_B, pop_C, pop_D)

for (j in 1:optimal_K){
  cluster_div = GenDiv_IND(pops[[j]])
  assign(paste0("ind_diver_pop_", name_clusters[j]), cluster_div)
  write.csv(cluster_div, file=paste0("Results/Step03/Diversity/Diversity_Individual_Clusters_", project_name, "_POP_", name_clusters[j], ".csv"))
}



#C. Comparing genetic diversity between clusters by ANOVA/Tukey's test:
# creating a df for analyses:
ind_diver_pop_A = ind_diver_pop_A %>% mutate(POP = "A", .before = Ho_O)
ind_diver_pop_B = ind_diver_pop_B %>% mutate(POP = "B", .before = Ho_O)
ind_diver_pop_C = ind_diver_pop_C %>% mutate(POP = "C", .before = Ho_O)
ind_diver_pop_D = ind_diver_pop_D %>% mutate(POP = "D", .before = Ho_O)
  
df_diversity =
  full_join(ind_diver_pop_A, ind_diver_pop_B) %>%
  full_join(ind_diver_pop_C) %>%
  full_join(ind_diver_pop_D) %>%
  mutate(POP = as.factor(POP))

head(df_diversity)
str(df_diversity)

#Perform ANOVA and Tukey
for(i in 2:ncol(df_diversity)){
  column =  names(df_diversity[i])
  result_anova = aov(df_diversity[,i] ~ POP, data= df_diversity)
  result_anova2 = summary(aov(df_diversity[,i] ~ POP, data= df_diversity))
  tk = TukeyHSD (result_anova)
  df = tibble::rownames_to_column(as.data.frame(tk$POP), "POP")
  write.csv(df, file=paste0("Results/Step03/Diversity/Comparing_Individual_Diversity_Clusters_", project_name, "_", column, ".csv"))
  write.csv(as.matrix(result_anova2[[1]]),file=paste0("Results/Step03/Diversity/Comparing_Individual_Diversity_Clusters_p-value_", project_name, "_", column, ".csv"))
  
  graph = ggplot(data=df, aes(y=POP, x=diff, xmin=lwr, xmax=upr))+
    geom_point() +
    geom_errorbarh(height=.2) +
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=1)+
    theme_bw() +
    theme_genetics +
    xlab("Differences") + ylab("Comparing Clusters")
  
  ggsave(filename = paste0("Results/Step03/Diversity/Comparing_Individual_Diversity_Clusters_", project_name, "_", column, ".pdf"), plot = graph, device = "pdf", dpi = 300)
  
}





#D. Plot bars graphs with highlights to representing clusters genetic metrics
#creating inputs for graphics:
df_highligth = df_diversity %>%
  melt(., id.vars=c("POP"))
facets_ind = c("Ho_O","He_O","He_E", "FIS")
colors_pop = c('#ffff00','#ffc0cb', "#ff0000", '#0000ff') #color for the genetic cluster
labels_pop = c(A = "POP A", B = "POP B", C = "POP C", D = "POP D") #name for the genetic cluster


#Graphics
ho_O = ggplot(df_highligth[df_highligth$variable %in% facets_ind[1],], aes(value, fill = POP))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ POP, nrow = 1, labeller = labeller(POP = labels_pop)) +
  scale_fill_manual (values = colors_pop) +
  xlab(expression(bold(Observed~Homozygosity~(Ho[OBS])))) + ylab("Frequency") +
  theme_bw()+
  theme_genetics
ho_O

he_O = ggplot(df_highligth[df_highligth$variable %in% facets_ind[2],], aes(value, fill = POP))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ POP, nrow = 1) +
  scale_fill_manual (values = colors_pop) +
  xlab(expression(bold(Observed~Heterozygosity~(He[OBS])))) + ylab("Frequency") +
  theme_bw()+
  theme_genetics +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
he_O


he_E = ggplot(df_highligth[df_highligth$variable %in% facets_ind[3],], aes(value, fill = POP))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ POP, nrow = 1) +
  scale_fill_manual (values = colors_pop) +
  xlab(expression(bold(Expected~Heterozygosity~(He[EXP])))) + ylab("Frequency") +
  theme_bw()+
  theme_genetics+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
he_E


fis = ggplot(df_highligth[df_highligth$variable %in% facets_ind[4],], aes(value, fill = POP))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ POP, nrow = 1) +
  scale_fill_manual (values = colors_pop) +
  xlab(expression(bold(Inbreeding~Coefficient~(F[IS])))) + ylab("Frequency") +
  theme_bw()+
  theme_genetics +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
fis


pdf("Results/Step03/Diversity/Histogram_IND_panel.pdf", onefile =F, width = 7 , height = 10)
ggarrange(ho_O, he_O, he_E, fis,
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 16, face = "bold", color ="black"),
          #vjust = 1, hjust = 1,
          ncol = 1, nrow = 4)
dev.off()







#### 4. GENETIC DISTANCE BY CLUSTERS----
#A. FST. ONLY IF YOUR K >= 2:
#Load as vcfR 
snps = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), verbose = T)

#Convert vcfR to genlight:
genlightdata = vcfR2genlight (snps)

# Define clusters as POP_ID
genlightdata@pop = as.factor(snps_neutral@meta$POP_ID)

# Run FST
fst = gl.fst.pop(genlightdata, nboots = 100, percent = 95, nclusters = 1)

#save results
write.csv(fst$Fsts, file=paste0("Results/Step03/Distance/FST_clusters_Fstas_", project_name, ".csv"))
write.csv(fst$Pvalues, file=paste0("Results/Step03/Distance/FST_clusters_pvalues_", project_name , ".csv"))
write.csv(fst$Bootstraps, file=paste0("Results/Step03/Distance/FST_clusters_bootstrap_", project_name , ".csv"))




#B. Jost's D, Gst's Nei, and Gst's Hedrick
#Convert vcfR to genind and define POP ID:
geninddata = vcfR2genind (snps)
pop(geninddata) = snps_neutral@meta$POP_ID
geninddata

#Jost's D 
D = pairwise_D(geninddata, linearized = FALSE)
write.csv(as.matrix(D), file = paste0("Results/Step03/Distance/D_clusters_", project_name,".csv"))

#Gst's Nei
GST_N = pairwise_Gst_Nei(geninddata, linearized = FALSE)
write.csv(as.matrix(GST_N), file = paste0("Results/Step03/Distance/GST_Nei_clusters_",project_name,".csv"))

#Gst's Hedrick
GST_H = pairwise_Gst_Hedrick(geninddata, linearized = FALSE)
write.csv(as.matrix(GST_H), file = paste0("Results/Step03/Distance/GST_Hend_clusters_", project_name, ".csv"))








#### 5. GENETIC DISTANCE BY INDIVIDUALS ----
#A.Load vcf if you need it
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.


#B. Convert VCF to LFMM
lfmm = Lfmm(snps_neutral, output.file=paste0("Results/Step03/Distance/Impute_LEA_", project_name, ".lfmm"))

#C. Run sNMF to calculate ancestry coefficients and ancestral genotype frequencies
optimal_K = 4
project.snmf = snmf(paste0("Results/Step03/Distance/Impute_LEA_", project_name, ".lfmm"), K = optimal_K, entropy = TRUE, repetitions = 10, project = "new", seed=123)

#D. You may load project file if it was saved in previous run
project.snmf = load.snmfProject(paste0("Results/Step03/Distance/Impute_LEA_", project_name, ".snmfProject"))

#E. select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.snmf, K = optimal_K))

#F. Impute the missing genotypes
impute(project.snmf, paste0("Results/Step03/Distance/Impute_LEA_", project_name, ".lfmm"), method = 'mode', K = optimal_K, run = best)

#G. Load imputed file
lfmmformat_imp = read.lfmm(paste0("Results/Step03/Distance/Impute_LEA_", project_name, ".lfmm_imputed.lfmm"))
lfmmformat_imp[1:10,1:10]
dim(lfmmformat_imp)

#H. Rename row and columns
length(snps_neutral@sample_id)
length(snps_neutral@site_id)
rownames(lfmmformat_imp) = snps_neutral@sample_id
colnames(lfmmformat_imp) = snps_neutral@site_id
lfmmformat_imp[1:10,1:10]




### 5.1. PCA-DISTANCE BASED ON BROKEN STICK RULE ----
#PCA-based distance and Dps: distance values among all pairs are calculated twice (below AND above diagonal) without distance within individuals
#A. Calculate PCA
PCAprcomp_imp = prcomp(lfmmformat_imp, center=TRUE, scale=TRUE)
str(PCAprcomp_imp)
summary(PCAprcomp_imp)

#B. Broken Stick PC values
screeplot(PCAprcomp_imp, bstick=TRUE, type="lines")
screeplot(PCAprcomp_imp, bstick=TRUE, type="barplot")
#3 PCs
summary(PCAprcomp_imp)
#4PCs are 9.22% of variance
n_pcs = 3

#C. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp = distance(PCAprcomp_imp$x[,1:n_pcs], method = "mahalanobis") # "mahalanobis" ou "euclidean", for multivariate space + than 2 PC, "mahalanobis" is better! 
#verify
PC_distLEA_prcomp
class(PC_distLEA_prcomp)
#convert to matrix
t_prcompLEA = as.matrix(PC_distLEA_prcomp)
t_prcompLEA[1:10,1:10]

#D. Save results
write.csv(t_prcompLEA, file=paste0("Results/Step03/Distance/PCA_Distance_BSR_IND_mahalanobis_", project_name, ".csv"))





###5.2. PCA-DISTANCE BASED ON 95% OF VARIANCE ----
#A. Verify the number of PC that explain 95% of variance
summary(PCAprcomp_imp)
#248PCs are 94.893% of variance or see Step 2 in DAPC action #E
n_pcs = 249

##B. Calculate PCA-based distance based on 95% of variance
PC_distLEA_95 = distance(PCAprcomp_imp$x[,1:n_pcs], method = "mahalanobis")
PC_distLEA_95
class(PC_distLEA_95)
t_prcomp95LEA = as.matrix(PC_distLEA_95)
t_prcomp95LEA[1:10,1:10]

#C. Save results
write.csv(t_prcomp95LEA, file=paste0("Results/Step03/Distance/PCA_Distance_95Var_IND_mahalanobis_", project_name, ".csv"))




###5.3. RELATEDNESS IN ALL SAMPLES ----
#Yang's Relatedness: distance values among all pairs are calculated once (below diagonal) with distance within individuals
#A. Calculate Relatedness
REL_YANG = Relatedness(snps_neutral, type = "yang",verbose = TRUE)
head(REL_YANG)
nrow(REL_YANG) 
colnames(REL_YANG) = c("INDV1","INDV2","RELATEDNESS_AJK_Yang")

#B. Exclude Relatedness between same individual
REL_YANG = REL_YANG[REL_YANG$INDV1 != REL_YANG$INDV2,]
#verify
head(REL_YANG)
nrow(REL_YANG) #1485 comparisons


#B. Save results:
write.csv(REL_YANG, file=paste0("Results/Step03/Distance/Yang_Reletedness_IND_neutral_", project_name, ".csv"))



###5.5. Dps DISTANCE ----
#A. Convert VCF to genind
snps =  read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), verbose = T)

#B. Convert data into genind
geninddata = vcfR2genind(snps)
geninddata@tab[1:10,1:10]

#C. Calculate proportions of shared alleles between pairs of individuals
dps_shared = propShared(geninddata)
head(dps_shared)

#D. Save results:
write.csv(dps_shared, file=paste0("Results/Step03/Distance/Dps_Shared_IND_neutral_", project_name, ".csv"))





###5.5. GRAPHICS FOR DISTANCE METRICS ----
#A. Metrics in a matrix:
df = as.data.frame(dps_shared[lower.tri(dps_shared, diag = F)])
names(df) = "Dps"
names(df)

p = ggplot(data=df, mapping = aes(x= Dps)) +
  geom_histogram(binwidth = 0.0005) +
  theme_classic()+
  theme_genetics+
  xlab ("Proportion of Shared Alleles") + ylab ("Frequency")

#B. Metrics in a df:
q = ggplot(data=REL_YANG, mapping = aes(x= RELATEDNESS_AJK_Yang)) +
  geom_histogram(binwidth = 0.0005) +
  theme_classic()+
  theme_genetics+
  xlab ("Relatedness Coefficient") + ylab ("Frequency")


#C. Choose one and save it as pdf:
pdf("Results/Step03/Distance/Histogram_relatedness_all.pdf", onefile = F)
q
dev.off()


#D. Genetic distance by Clusters
#Relatedness
pops #Clusters were defined

# relatedness by cluster
for (j in 1:optimal_K){
  REL_POP = Relatedness(pops[[j]], type = "yang",verbose = TRUE)
  colnames(REL_POP) = c("INDV1","INDV2","RELATEDNESS_AJK_Yang")
  REL_POP = REL_POP[REL_POP$INDV1 != REL_POP$INDV2,]
  write.csv(REL_POP, file=paste0("Results/Step03/Distance/Yang_Reletedness_CLUSTERS_neutral_", project_name, "_POP_", name_clusters[j], ".csv"))
  
  q2 = ggplot(data=REL_POP, mapping = aes(x= RELATEDNESS_AJK_Yang)) +
    geom_histogram(binwidth = 0.0005) +
    theme_classic()+
    theme_genetics+
    xlab ("Relatedness Coefficient") + ylab ("Frequency")
  
  ggsave(filename = paste0("Results/Step03/Distance/Yang_Reletedness_CLUSTERS_", project_name, "_POP_", name_clusters[j], ".pdf"), plot = q2, device = "pdf", dpi = 300)

}

#Dps or another matrix metric:
#define a list of positions for each cluster
clusters = list(pop_POSI_A, pop_POSI_B, pop_POSI_C, pop_POSI_D)

class(dps_shared)

for (j in 1:optimal_K){
  positions = clusters[[j]]
  mt = dist_subset(dps_shared, snps_neutral@sample_id[positions])
  df = as.data.frame(mt[lower.tri(mt, diag = F)])
  names(df) = "Dps"

  p2 = ggplot(data=df, mapping = aes(x= Dps)) +
    geom_histogram(binwidth = 0.0005) +
    theme_classic()+
    theme_genetics+
    xlab ("Proportion of Shared Alleles") + ylab ("Frequency")
  
  ggsave(filename = paste0("Results/Step03/Distance/Dps_Shared_CLUSTERS_", project_name, "_POP_", name_clusters[j], ".pdf"), plot = p2, device = "pdf", dpi = 300)
}









#### 6. TAJIMA D BY CLUSTERS ----
#A. Load the VCF file:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.

#B. Filter one SNP per contig. This thins SNPs to a given distance in bp from one another. Setting the distance higher than the length of the contig ensures that you'll have a single SNP per contig.
snps_thin = Filter(snps_neutral, filterOptions(thin=300)) 
VCFsummary(snps_thin) #277 individuals and 3412 SNPs.


#C. Subset vcf by cluster
clusters = list(pop_POSI_A, pop_POSI_B, pop_POSI_C, pop_POSI_D)
optimal_K = 4
name_clusters = c("A", "B", "C", "D")

#subset
for (i in 1:optimal_K){
  UNIND = snps_thin@sample_id[clusters[[i]]]
  pop = Subset(snps_thin, samples=UNIND)
  assign(paste0("popD_", name_clusters[i]), pop)
}
#verify
VCFsummary(popD_A) #47 individuals and 3412 SNPs.
VCFsummary(popD_B) #47 individuals and 3412 SNPs.
VCFsummary(popD_C) #65 individuals and 3412 SNPs.
VCFsummary(popD_D) #118 individuals and 3412 SNPs.


#D. ESTIMATE TAJIMA'S D, BIAS-CORRECTED FOR MAF
tajd_pA = TajimaD(popD_A, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_pA$results
str(tajd_pA$simulations)
save(tajd_pA, file = "./Results/Step03/TajimaD/tajd_pA.Rdata")

tajd_pB = TajimaD(popD_B, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_pB$results
str(tajd_pB$simulations)
save(tajd_pB, file = "./Results/Step03/TajimaD/tajd_pB.Rdata")

tajd_pC = TajimaD(popD_C, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_pC$results
str(tajd_pC$simulations)
save(tajd_pC, file = "./Results/Step03/TajimaD/tajd_pC.Rdata")

tajd_pD = TajimaD(popD_D, nboot=10000, maf=0.05, use_vcftools_D=FALSE)
tajd_pD$results
str(tajd_pD$simulations)
save(tajd_pD, file = "./Results/Step03/TajimaD/tajd_pD.Rdata")


#E. Load results 
load("./Results/Step03/TajimaD/tajd_pA.Rdata")
load("./Results/Step03/TajimaD/tajd_pB.Rdata")
load("./Results/Step03/TajimaD/tajd_pC.Rdata")
load("./Results/Step03/TajimaD/tajd_pD.Rdata")


#F. Save results in a table
A_res = cbind.data.frame(tajd_pA$results$`Tajima D, bias-corrected`, tajd_pA$results$`Tajima D, bias-corrected 95CI`[1], tajd_pA$results$`Tajima D, bias-corrected 95CI`[2], tajd_pA$results$`Pr(abs(D) > abs(D_observed)|Null)`) %>%
  setNames(c("Tajima's D", "Lower", "Upper", "p-value"))

B_res = cbind.data.frame(tajd_pB$results$`Tajima D, bias-corrected`, tajd_pB$results$`Tajima D, bias-corrected 95CI`[1], tajd_pB$results$`Tajima D, bias-corrected 95CI`[2], tajd_pB$results$`Pr(abs(D) > abs(D_observed)|Null)`)  %>%
  setNames(c("Tajima's D", "Lower", "Upper", "p-value"))

C_res = cbind.data.frame(tajd_pC$results$`Tajima D, bias-corrected`, tajd_pC$results$`Tajima D, bias-corrected 95CI`[1], tajd_pC$results$`Tajima D, bias-corrected 95CI`[2], tajd_pC$results$`Pr(abs(D) > abs(D_observed)|Null)`) %>%
  setNames(c("Tajima's D", "Lower", "Upper", "p-value"))

D_res = cbind.data.frame(tajd_pD$results$`Tajima D, bias-corrected`, tajd_pD$results$`Tajima D, bias-corrected 95CI`[1], tajd_pD$results$`Tajima D, bias-corrected 95CI`[2], tajd_pD$results$`Pr(abs(D) > abs(D_observed)|Null)`) %>%
  setNames(c("Tajima's D", "Lower", "Upper", "p-value"))

res = rbind(A_res, B_res, C_res, D_res)
rownames(res) = c("A", "B", "C", "D")
head(res)
write.csv(res, "./Results/Step03/TajimaD/Tajima_results.csv")



#G. Plot null distribution (histogram) is shown next to the observed Tajima's D value (red line)
options(scipen = 9999)
dfA = 
  cbind.data.frame(tajd_pA$simulations$`Null, bias-corrected`, tajd_pA$simulations$`Bootstrap, bias-corrected`) %>%
  setNames(., c("Nullc", "Boostrapc")) %>%
  mutate(., POP = "A")
head(dfA)

dfB = 
  cbind.data.frame(tajd_pB$simulations$`Null, bias-corrected`, tajd_pB$simulations$`Bootstrap, bias-corrected`) %>%
  setNames(., c("Nullc", "Boostrapc")) %>%
  mutate(., POP = "B")
head(dfB)

dfC = 
  cbind.data.frame(tajd_pC$simulations$`Null, bias-corrected`, tajd_pC$simulations$`Bootstrap, bias-corrected`) %>%
  setNames(., c("Nullc", "Boostrapc")) %>%
  mutate(., POP = "C")
head(dfC)

dfD = 
  cbind.data.frame(tajd_pD$simulations$`Null, bias-corrected`, tajd_pD$simulations$`Bootstrap, bias-corrected`) %>%
  setNames(., c("Nullc", "Boostrapc")) %>%
  mutate(., POP = "D")
head(dfD)

df_taj = 
  full_join (dfA, dfB) %>%
  full_join(dfC) %>%
  full_join(dfD)

plot_taj = ggplot(df_taj, aes(x = Nullc, fill = POP)) +
  geom_histogram(binwidth=0.005) +
  geom_vline(xintercept=0) +
  geom_vline(data=filter(df_taj, POP=="A"), aes(xintercept=mean(Boostrapc)), lty=2, col="red") + 
  geom_vline(data=filter(df_taj, POP=="B"), aes(xintercept=mean(Boostrapc)), lty=2, col="red") + 
  geom_vline(data=filter(df_taj, POP=="C"), aes(xintercept=mean(Boostrapc)), lty=2, col="red") + 
  geom_vline(data=filter(df_taj, POP=="D"), aes(xintercept=mean(Boostrapc)), lty=2, col="red") + 
  facet_wrap(~ POP, nrow = 1) +
  scale_fill_manual (values = c("gray", "gray", "gray", "gray"))+
  theme_bw() +
  theme_genetics +
  theme(legend.position = "none") +
  ylab("Frequency") +xlab ("Tajima's D")
plot_taj 

pdf("./Results/Step03/TajimaD/TajimaD_simulations_graphs.pdf", onefile = F)
plot_taj
dev.off()






#### 7. CONVERTING VCF TO GENEPOP TO RUN NeEstimator -----
#A. LOAD VCF
snpsR = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), verbose = T)

#B. DEFINING POPULATIONS USING THE VCF METAFILE:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.
population = as.factor(snps_neutral@meta$POP_ID)

#C. Converting VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#D. Adding strata (pops) into Genind
snps_genind@pop = population

#E. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

#F. Converting Gtypes to GENEPOP to run in NeEstimator and save it:
genepopWrite(snps_gtypes, "Results/Step03/Ne/Genepop_NeEstimator_Neutral")


#G. Running Ne in StrataG. StrataG did not consider missing data.
Ne = ldNe(snps_gtypes, maf.threshold = 0, by.strata = TRUE, ci = 0.95, drop.missing = TRUE, num.cores = 4)
Ne

#H. Save the results:
write.csv(Ne, paste0("Results/Step03/Ne/NeLD_Nomissing_", project_name,  ".csv"))


##END;
