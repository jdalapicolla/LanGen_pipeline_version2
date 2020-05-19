###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
#######################   STEP 02: GENETIC STRUCTURE   ########################

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
#A. ESTIMATE THE NUMBER OF GENETIC CLUSTERS IN THE DATASET, USING:
  #A1: CLUSTERING METHOD - sNMF (STRUCTURE-LIKE METHOD)
  #A2: NAIVE CLUSTERING MODEL-FREE METHOD - DAPC
  #A3: MODEL-FREE METHOD - PCA
  #A4: SPATIALLY-EXPLICITY CLUSTERING METHODS - TESS

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
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")
} else {print (paste0("'tidyverse' has already been installed in library"))}
if("tess3r" %in% rownames(installed.packages()) == FALSE){devtools::install_github("bcm-uga/TESS3_encho_sen")
} else {print (paste0("'tess3r' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("maps" %in% rownames(installed.packages()) == FALSE){install.packages("maps")
} else {print (paste0("'maps' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vcfR, dartR, adegenet, tidyverse, tess3r, ggplot2, raster, maps)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_snmf", "./Results_DAPC", "./Results_PCA", "./Results_TESS"))


#------------------------------------------------------------------------------
#                       1. Genetic Structure using sNMF 
#------------------------------------------------------------------------------
##1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "pilocarpus"

###1.2. POPULATION ASSIGNMENT ANALYSIS USING sNMF WITHOUT FST OUTLIERS:
#Now we will carry out population assignment again, but using only neutral loci dataset.
#A1. Load neutral .vcf file after remove the outlier SNPs
snps_fil_ldF_neutral <- vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID = T) 
VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.

##B. Convert to geno object.You need to specify the output file. It will automatically subset the vcf file and assign it as a new object
snps_fil_ldF_neutral <- Geno(snps_fil_ldF_neutral, output.file = paste0("vcf/", project_name, "_filtered_neutral.geno"))
VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.

#C. Create folders for alpha values and copy .geno object in each folder:
alpha_values = c(10, 100, 500, 1000, 2000, 4000)
for (i in alpha_values){
  path = paste0("./Results_snmf/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("vcf/", project_name, "_filtered_neutral.geno"), path )
}

#D. Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # numeber of replication by K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores

#E. Run sNMF (LEA) using different alpha.
set.seed(123)
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}

#F. To load the SNMF projects in a new R session (after quitting R), use:  project = load.snmfProject("Alpha1000//pilocarpus_filtered_ld_hw.snmfProject") ##This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop, "n"), pro_snmf)
}

#G. summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
summary(project5n)
summary(project6n)

#H. View Cross-Entropy plots
PlotK(project1n) #4
PlotK(project2n) #4
PlotK(project3n) #4
PlotK(project4n) #4
PlotK(project5n) #3
PlotK(project6n) #4
#another way to view the results
#plot(project1n, lwd = 5, col = "red", pch=1)
#plot(project2n, lwd = 5, col = "red", pch=1)
#plot(project3n, lwd = 5, col = "red", pch=1)
#plot(project4n, lwd = 5, col = "red", pch=1)
#plot(project5n, lwd = 5, col = "red", pch=1)
#plot(project6n, lwd = 5, col = "red", pch=1)

#I. Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results_snmf/Cross_Entropy_sNMF_Alpha_",  i, "n.pdf"), onefile = F)
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
  }

#J. Select optimal K value
optimal_K = 4

#K. Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n, p4=project4n, p5=project5n, p6=project6n)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("Results_snmf/Alpha", alpha_values[as.numeric(best_run_split[6])],"n/", project_name, "_filtered_neutral.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])

#L. Barplot replace best run information
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
order = barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
write.table(order[[1]], "./Results_Metafiles/snmf_bestK_pipegraph_order_neutral.txt")

pdf("./Results_snmf/sNMF_Pipegraph_neutral_version2.pdf", onefile = F)
plot.new()
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
dev.off()

pdf("./Results_snmf/sNMF_Pipegraph_neutral.pdf", onefile = F)
my.colors <- rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1, cex.axis = .3)
dev.off()

#M. Add admixture coefficient and replace the population ID to vcf file
Qmat <- as.data.frame(Q(project, run=run, K=optimal_K))
head(Qmat)

columns = c() 

for (i in 1:ncol(Qmat)){
  columns[i] = paste0("Adx_Coeff_", i)
}

colnames(Qmat) = columns
head(Qmat)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:optimal_K){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat[i]
  }
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#N. Verify individuals and maximum Adx_Coeff and add to metafile:
popIds = apply(Qmat, 1, which.max)
snps_fil_ldF_neutral@meta$PopID_snmf <- popIds
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#O. Save new vcf file with sNMF results and pop ID
Save(snps_fil_ldF_neutral, paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf")) #save vcf
write.csv(snps_fil_ldF_neutral@meta, paste0("Results_Metafiles/ancestry_coef_LEA_", project_name, ".csv"), quote = F) #save result as table
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_snmf), file= paste0("Results_Metafiles/", project_name, "_neutral_popsID_LEA.csv")) #save only pop ID from sNMF

VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.


#------------------------------------------------------------------------------
#                       2. Genetic Structure using DAPC 
#------------------------------------------------------------------------------    
#DA Discriminant Analysis focus on between-group variability, while neglecting 
#within-group variation. this method also allows for a probabilistic assignment
#of individuals to each group, as in Bayesian clustering methods. the method
#requires the number of variables (alleles) to be less than the number of
#observations (individuals). This condition is generally not fulfilled in Single
#Nucleotide Polymorphism (SNP) or re-sequencing datasets. Second, it is hampered
#by correlations between variables, which necessarily occur in allele frequencies
#due to the constant-row sum constraint [i.e., compositional data].
#Uncorrelated variables will be even more blatant in the presence of linkage disequilibrium
#
#DAPC relies on data transformation using PCA as a prior step to DA, which ensures
#that variables submitted to DA are perfectly uncorrelated, and that their number is
#less than that of analysed individuals (1:10 proportion). Without implying a necessary 
#loss of genetic information, this transformation allows DA to be applied to any genetic data.
#
#K-means relies on the same model as DA to partition genetic variation into a between-group
#and a within-group component, and attempts to find groups that minimize the latter. We use
#Bayesian Information Criterion (BIC) to assess the best supported model, and therefore the
#number and nature of clusters.
#-------------------------------------------------------------------------------------

###2.1. POPULATION ASSIGNMENT ANALYSIS USING DAPC:
#A. Load neutral .vcf file after removing outlier SNPs:
vcf = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral.vcf"), verbose = FALSE)

#B. Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input

#C. Perform a PCA to choose the number of PC in the DAPC:
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE)

#D. % of PC variation
pc_pca = as.data.frame(pca_input$eig)
pc_pca[,2] = (pc_pca/sum(pc_pca))*100

#D1. Rule of 100% of variance:
index_100 = length(rownames(input@tab))-1
index_100 # number of PC to reach to 100% of explained variance

#D2. Rule of at least 95% of variance:
index_95 = length(which(cumsum(pc_pca[,2]) <= 95))
index_95

#D3. Rule of at least 70% of variance:
index_70 = length(which(cumsum(pc_pca[,2]) <= 70))
index_70 

#D4. Rule of minimum variance per PCs:
variance_pc = 100/(nrow(input@tab)-1)
variance_pc #PCs that increase the explained variance bellow this threshould will be removed
#calculate number of PCs
index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
index_min

#E. Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust=40)
index_100 #276 PCs - k = 4 #For me, 100% of variation is more coherent.
index_95  #248 PCs - k = 4
index_70  #147 PCs - k = 5
index_min #88 PCs - k = 7
#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size

#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
pdf("./Results_DAPC/Bestk_DAPC.pdf", onefile = T)
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Best K")
dev.off()

#G. ATTENTION, if your dataset is k = 1 go to PCA analysis (Action #4)! Graphs for DAPC need datasets with K >= 2

#H. Choose the best PCs number to recover correctly the clusters. The input should be a Genind object without missing data. Maximum number of PCs is number of individuals -1. Replication by deafaul is 30. Save automatically the graph
pdf("./Results_DAPC/Best_PCs_Number_DAPC.pdf", onefile = T)
number_PCs = xvalDapc(tab(input, NA.method = "mean"), grp$grp, scale = T, n.pca.max = (nrow(input@tab)-1))
dev.off()

#I. Verify the number of PCs and DA used and summary of DAPC
number_PCs$DAPC$n.pca
number_PCs$DAPC$n.da
summary(number_PCs$DAPC)

#J. Verify scatter plot and define the group colors and names:
# color and names
col_dapc= c('red', 'blue', 'yellow', "darkslategrey")
legend_dapc = c("POP1", "POP2", "POP3", "POP4")
#plot graph
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = col_dapc, txt.leg = legend_dapc,
        clabel = FALSE, posi.leg = "bottomright", scree.pca = TRUE, pch=19:20,
        posi.pca = "bottomleft", posi.da = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

#K. Save scatterplots. Please edit the position of caption
pdf("./Results_DAPC/scatter_DAPC.pdf", onefile = T)
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = col_dapc, txt.leg = legend_dapc,
        clabel = FALSE, posi.leg = "bottomright", scree.pca = TRUE, pch=19:20,
        posi.pca = "bottomleft", posi.da = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dev.off()

#L. Save a STRUCTURE-like graph 
pdf("./Results_DAPC/pipegraph_DAPC.pdf", onefile = T)
compoplot(number_PCs$DAPC, col=col_dapc, legend=TRUE, txt.leg=legend_dapc, cleg=.8)
dev.off()

#M. Add DAPC clusters ID in the VCF file
snps_fil_ldF_neutral@meta$PopID_DAPC = as.character(grp$grp)
head(snps_fil_ldF_neutral@meta)

#N. Add DAPC posterior probabilities in the VCF file:
Qmat_2 <- as.data.frame(number_PCs$DAPC$posterior)
head(Qmat_2)
columns = c() 
for (i in 1:ncol(Qmat_2)){
  columns[i] = paste0("DAPC_Posterior_", i)
}
colnames(Qmat_2) = columns
head(Qmat_2)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:length(grp$size)){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat_2[i]
}
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#O. Save new vcf file and pop ID file
Save(snps_fil_ldF_neutral, paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC.vcf"))
write.csv(snps_fil_ldF_neutral@meta, paste0("Results_Metafiles/ancestry_coef_LEA_DAPC_", project_name, ".csv"))
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_DAPC), file= paste0("Results_Metafiles/", project_name , "_neutral_popsID_DAPC.csv")) #edit column position by the number of coefficients.


#------------------------------------------------------------------------------
#                       3. Genetic Structure using TESS3 
#------------------------------------------------------------------------------
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html

###3.1. INPUT FILES FOR TESS:
#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 5268 SNPs.

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


###3.2. RUNNING THE TESS3R FUNCTION:
#A. Costumize values for the run
lambda_values = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5) #test lambda values around 1.
K = c(1:10) # set the number of K to be tested
replications = 5 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel
mask = 0.05 #porportion of masked values

#B.Run a loop for all alpha values
set.seed(13)
for (i in lambda_values){
  tess3.ls = tess3(genotypes, coord = coordinates, K = K, mask = mask, lambda = i,
                   method = "projected.ls", max.iteration = 5000, rep = replications,
                   ploidy = ploidy, openMP.core.num = CPU)
  save(tess3.ls, file = paste0("./Results_TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
  
}

#C. Choose the lambda with the minimum cross-validation value:
#create a matrix to save results
cross_value = matrix(NA, max(K)*length(lambda_values), 3)
colnames(cross_value) = c("K", "Lambda", "Crossvalidation")

loop = 0 #Always set loop = 0.
for (i in lambda_values){
  load(file = paste0("./Results_TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
  for (j in 1:max(K)){
    loop=loop+1
    res = Gettess3res(tess3.ls, K=j)
    cross_value[loop,1] = j
    cross_value[loop,2] = i
    cross_value[loop,3] = min(res$crossvalid.crossentropy)}
}
#save as csv
write.csv(cross_value,paste0("./Results_TESS/Tess3ls_Crossvalidation_values_", project_name, ".csv"))
#choose best lambda
lambda_tess = as.vector(cross_value[cross_value[,3] == min(cross_value[,3]), ][2])
lambda_tess

#D. Choose best K:
#load the best lambda project:
load(file = paste0("./Results_TESS/Tess3ls_Lambda_", lambda_tess,"_", project_name, ".RData"))
#plot results
pdf(paste0("./Results_TESS/TESS3_PlotK_CrossValidation_Lambda_", lambda_tess, ".pdf"), onefile =F)
plot.new()
plot(tess3.ls, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")
dev.off()

#E.Choose the Best K and create the Q-matrix
#Best K
optimal_K = 5
res_tess3 = Gettess3res(tess3.ls, K=optimal_K)

#F. Structure-like barplot for the Q-matrix
my.colors = c('red', 'pink', "darkslategrey", 'blue', 'yellow') #choose the colors according to the number of K
my.palette = CreatePalette(my.colors, 9)
pdf("./Results_TESS/TESS_Ancestry.pdf", onefile =F)
plot.new()
barplot(res_tess3$Q, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col.palette = my.palette) -> bp
axis(1, at = 1:nrow(res_tess3$Q), labels = bp$order, las = 3, cex.axis = .4)
dev.off()

#G. Spatial interpolation of ancestry coefficient
pdf("./Results_TESS/TESS_MAP.pdf", onefile =F)
plot.new()
plot(res_tess3$Q, coordinates, method = "map.max", interpol = FieldsKrigModel(10),
     main = "Ancestry coefficients",cex.main=1.5,
     xlab = "Longitude", ylab = "Latitude",cex.lab=1.3,
     resolution = c(700,700), cex = 1.5,
     col.palette = my.palette)
dev.off()

#H. Add admixture coefficient and replace the population ID to vcf file
Qmat_tess = as.data.frame(res_tess3$Q)
head(Qmat_tess)

columns = c() 

for (i in 1:ncol(Qmat_tess)){
  columns[i] = paste0("Adx_Coeff_TESS_", i)
}

colnames(Qmat_tess) = columns
head(Qmat_tess)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:optimal_K){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat_tess[i]
}
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#verify individuals and maximum Adx_Coeff and add to metafile:
popIds = apply(Qmat_tess, 1, which.max)
snps_fil_ldF_neutral@meta$PopID_tess <- popIds
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#I. Save new vcf file with TESS results and pop ID
Save(snps_fil_ldF_neutral, paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf")) #save vcf
write.csv(snps_fil_ldF_neutral@meta, paste0("Results_Metafiles/ancestry_coef_LEA_DAPC_TESS_", project_name, ".csv"), quote = F) #save result as table
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_snmf), file= paste0("Results_Metafiles/", project_name, "_neutral_popsID_LEA_DAPC_TESS.csv")) #save only pop ID from sNMF

VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.


#------------------------------------------------------------------------------
#                       4. Genetic Structure using PCA 
#------------------------------------------------------------------------------
#PCA is its ability to identify genetic structures in very large datasets within
#negligible computational time, and the absence of any assumption about the
#underlying population genetic model. PCA aims to summarize the overall variability
#among individuals, which includes both the divergence between groups
#(i.e., structured genetic variability), and the variation occurring within groups
#(‘random’ genetic variability).
#--------------------------------------------------------------------------------

###4.1. RUN THE PCA:
#A. You ran the PCA before for DPAC:
pca_input

#B. % of contribution for three first PCs
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")

#C. Min and max values for axes
pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))

#D. save results as .CSV
#coordenates in the PCA by axis
write.csv(pca_input$li, file = "Results_PCA/Axis_coord.csv")
#% of PC variation
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "Results_PCA/contribution_pc_eig.csv")


###4.2. PCA GRAPHS
#open the PDF
pdf(paste0("./Results_PCA/PCA_", project_name,".pdf"))
plot.new()
#margins of a single graph (unit = lines width)
par(mar=c(5,5,1,1))
#size of the graph
plot.window(xlim=pc1_range, ylim=pc2_range)

#number of samples 
samples = nrow(input@tab)
#plot all points in blue. You can change the color in bg = 'color'
points(pca_input$li[1:samples,1], pca_input$li[1:samples,2], col = 'white', bg = "black", cex = 1.5, pch=21) 

#draw the axis; put the limits similar to PC Min and Max values
axis(1, at=seq((round(pc1_range[1],-1)-10), (round(pc1_range[2],-1)+10), by=30), cex.axis=1.15);
axis(2, at=seq((round(pc2_range[1],-1)-10), (round(pc2_range[2],-1)+10), by=30), cex.axis=1.15, las=1);

#complete the name of each axis
mtext(side=1, text=PC1, line=2.5, cex=1)
mtext(side=2, text=PC2, line=2.8, cex=1)

#labeling all individuals on the plot
pop.names = row.names(input@tab)
#if you need to remove some pattern in the individuals name
#pop.names = str_remove(row.names(input@tab), pattern="_sorted")
#plot the text
#text(pca_input$li[1:samples,1], pca_input$li[1:samples,2], labels = pop.names[1:samples], pos = 4, cex = 0.7)

#finishing the graph
dev.off()


###4.3. TRACY-WIDOM TEST FOR EIGENVALUES
#Tracy CA and Widom H. (1994). Level spacing distributions and the bessel kernel. Commun Math Phys. 161 :289--309. Patterson N, Price AL and Reich D. (2006). Population structure and eigenanalysis. PLoS Genet. 2 :20.
#A. Load the VCF file after DAPC step
snps_fil_ldF_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC.vcf"), overwriteID = T)
VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.
#verify the metadata 
snps_fil_ldF_neutral@meta

#B. Convert VCF to LFMM
lfmm_pca = Lfmm(snps_fil_ldF_neutral, output.file=paste0("./Results_PCA/PCA_LFMM_", project_name, "_filtered_neutral.lfmm"))

#C. Create a pcaProject object:
pc = LEA::pca(paste0("./Results_PCA/PCA_LFMM_", project_name, "_filtered_neutral.lfmm"), scale = TRUE, center=TRUE)

#C. Perform Tracy-Widom tests
tw = tracy.widom(pc)

#D. Display the p-values for the Tracy-Widom tests. 
tw$pvalues
#24 PCs are <0.05. K = PC -1. So my K = 23 

#E. Save results
write.csv(tw,paste0("./Results_PCA/PCA_Tracy-Widom_", project_name, ".csv"))

###4.4. PLOT PC1 VALUES IN A QQPLOT TO TEST NORMALITY
pdf("./Results_PCA/PCA_Normality_PC1_qqplot.pdf", onefile =F)
plot.new()
qqnorm(pca_input$li$Axis1, pch = 1, frame = FALSE)
qqline(pca_input$li$Axis1, col = "black", lwd = 2)
dev.off()


#------------------------------------------------------------------------------
#                       5. Comparing Genetic Clusters 
#------------------------------------------------------------------------------
###5.1. DO THE CLUSTERS (K) FROM DIFFERENT ANALYSES HAVE THE SAME COMPOSITION? ONLY IF BOTH ANALYSES WERE K >= 2:
#A. Position of samples by population by sNMF approach:
for (i in 1:length(unique(snps_fil_ldF_neutral@meta$PopID_snmf))){
  pop = which(snps_fil_ldF_neutral@meta$PopID_snmf == i)
  assign(paste0("pop_snmf_", i), pop)
}

#B. Position of samples by populion by DAPC approach:
for (i in 1:length(unique(snps_fil_ldF_neutral@meta$PopID_DAPC))){
  pop = which(snps_fil_ldF_neutral@meta$PopID_DAPC == i)
  assign(paste0("pop_DAPC_", i), pop)
}

#C. Position of samples by populion by TESS approach:
for (i in 1:length(unique(snps_fil_ldF_neutral@meta$PopID_tess))){
  pop = which(snps_fil_ldF_neutral@meta$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}

#D. PCA approach do not have a assignment tools but you can analyze the clusters by the pdf saved in action #5 if the number of sample is small.

#E. Compare the sample size. Add the number of pop that you have:
pop_size_snmf = c(length(pop_snmf_1), length(pop_snmf_2), length(pop_snmf_3), length(pop_snmf_4))
pop_size_tess = c(length(pop_TESS_1), length(pop_TESS_2), length(pop_TESS_3), length(pop_TESS_4), length(pop_TESS_5))

pop_size_snmf # 118  65  47  47
grp$size #65 118  47  47
pop_size_tess #65  15  47 103  47

#F. Check if identical samples in each group. If they are corrected, function will return "TRUE"
#POPs with the same sample size #BLUE
identical(pop_snmf_1, pop_DAPC_2, sort(c(pop_TESS_2, pop_TESS_4))) #TRUE

#POPs with the same sample size in three or more analyses #RED
identical(pop_snmf_2, pop_DAPC_1, pop_TESS_1) #TRUE

#POPs with the same sample size #GRAY
identical(pop_snmf_3, pop_DAPC_4, pop_TESS_3) #TRUE

#POPs with the same sample size #YELLOW
identical(pop_snmf_4, pop_DAPC_3, pop_TESS_5) #TRUE

#G. Define same colors for the same POP:
col_dapc= c('red', 'blue', 'yellow', "darkslategrey")
legend_dapc = c("POP1", "POP2", "POP3", "POP4")

col_snmf= c('blue', 'red', 'darkslategrey', "yellow")
legend_snmf = c("POP1", "POP2", "POP3", "POP4")

col_tess= c('red', 'pink', "yellow", 'blue', 'darkslategrey')
legend_tess = c("POP1", "POP2", "POP3", "POP4", "POP5")


#H. If the results are false, check the difference between files, which samples are present in a file and not in other. Invert the order of population. If the 1st pop is inside the 2nd pop function will return "character (0)"
setdiff(pop_snmf_2, pop_DAPC_3)
setdiff(pop_DAPC_3, pop_snmf_2)

# Choose one number of clusters (K) for next steps!

##END
