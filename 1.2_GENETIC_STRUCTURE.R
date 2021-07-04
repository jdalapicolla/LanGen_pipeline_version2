
####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
######################### STEP 02: GENETIC STRUCTURE ############################




### Script prepared by Jeronymo Dalapicolla, Carolina S. Carvalho, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###




#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline



##2. GOALS FOR THIS STEP:
#A. ESTIMATE THE NUMBER OF GENETIC CLUSTERS IN THE DATASET, USING:
#A1: CLUSTERING METHOD - sNMF (STRUCTURE-LIKE METHOD)
#A2: NAIVE CLUSTERING MODEL-FREE METHOD - DAPC
#A3: MODEL-FREE METHOD - PCA
#A4: SPATIALLY-EXPLICITY CLUSTERING METHODS - TESS


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
if (!require('tess3r'))       devtools::install_github("bcm-uga/TESS3_encho_sen");   library('tess3r')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('maps'))         install.packages("maps");              library('maps')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2')
if (!require('ggrepel'))      install.packages("ggrepel");           library('ggrepel')
if (!require('rworldmap'))    install.packages("rworldmap");         library('rworldmap')



#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, LEA, vcfR, dartR, adegenet, tidyverse, tess3r, ggplot2, raster, maps, reshape2, ggrepel, rworldmap)

##7. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS ----
create_dir(c("./Results/Step02/Metafiles", "./Results/Step02/PCA", "./Results/Step01/DAPC", "./Results/Step02/sNMF", "./Results/Step02/TESS"))



##8. CREATE A PATTERN FOR GRAPHIC FOLDERS S TO SAVE THE RESULTS ----
theme_genetics = theme(axis.text=element_text(size=10, color="black"), #text in ticks axes
                       axis.title=element_text(size=12, face="bold"), #label axes
                       axis.line = element_line(colour = "black", size = 1, linetype = "solid"), #line on axes
                       axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"), #line on ticks
                       axis.ticks.length = unit(.25, "cm"), #ticks length
                       axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), #space between axis and label
                       axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), #space between axis and label
                       panel.grid.major = element_blank(), # remove grids
                       panel.grid.minor = element_blank(), # remove grids
                       panel.background = element_blank(), # remove background
                       panel.border = element_blank()) # remove borders)  





#### ANALYSIS ---- 


#### 1. LOAD FILES -----
##1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP01:
#A. Project name:
project_name = "pilocarpus"


#B. Load neutral .vcf file after remove the outlier SNPs
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID = T) 
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.
head(snps_neutral@meta)








#### 2.GENETIC STRUCTURE USING sNMF ---- 
##A. Convert to genotype object. You need to specify the output file. It will automatically subset the vcf file and assign it as a new object
snps_fil_ldF_neutral = Geno(snps_neutral, output.file = paste0("vcf/", project_name, "_filtered_neutral.geno"))
VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.


#B. Create folders for alpha values and copy .geno object in each folder:
alpha_values = c(10, 100, 500, 1000, 2000, 4000)
for (i in alpha_values){
  path = paste0("./Results/Step02/sNMF/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("vcf/", project_name, "_filtered_neutral.geno"), path )
}


#C. Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # number of replication by K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores


#D. Run sNMF (LEA) using different alpha.
set.seed(123)
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results/Step02/sNMF/Alpha", i,"n/", project_name, "_filtered_neutral.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}


#E. To load the SNMF projects in a new R session (after quitting R).
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results/Step02/sNMF/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop, "n"), pro_snmf)
}

#F. summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
summary(project5n)
summary(project6n)

#G. View Cross-Entropy plots
PlotK(project1n) #6
PlotK(project2n) #6
PlotK(project3n) #4-5
PlotK(project4n) #4
PlotK(project5n) #4
PlotK(project6n) #4
#another way to view the results
#plot(project1n, lwd = 5, col = "red", pch=1)
#plot(project2n, lwd = 5, col = "red", pch=1)
#plot(project3n, lwd = 5, col = "red", pch=1)
#plot(project4n, lwd = 5, col = "red", pch=1)
#plot(project5n, lwd = 5, col = "red", pch=1)
#plot(project6n, lwd = 5, col = "red", pch=1)


#H. Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results/Step02/sNMF/Cross_Entropy_sNMF_Alpha_",  i, "n.pdf"), onefile = F)
  path = paste0("Results/Step02/sNMF/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
  }


#I. Select optimal K value
optimal_K = 4


#J. Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n, p4=project4n, p5=project5n, p6=project6n)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("Results/Step02/sNMF/Alpha", alpha_values[as.numeric(best_run_split[6])],"n/", project_name, "_filtered_neutral.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])


#K. Barplot replace best run information
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
order = barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
write.table(order[[1]], "./Results/Step02/Metafiles/snmf_bestK_pipegraph_order_neutral.txt")

pdf("./Results/Step02/sNMF/sNMF_Pipegraph_neutral_version2.pdf", onefile = F)
plot.new()
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
dev.off()

pdf("./Results/Step02/sNMF/sNMF_Pipegraph_neutral.pdf", onefile = F)
my.colors = rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1, cex.axis = .3)
dev.off()


#L. Add admixture coefficient and replace the population ID to vcf file
Qmat =
  project %>%
  Q(., run=run, K=optimal_K) %>%
  as.data.frame() %>%
  setNames(c("Adx_Coeff_1", "Adx_Coeff_2", "Adx_Coeff_3","Adx_Coeff_4")) %>%
  mutate(PopID_snmf = apply(., 1, which.max))
head(Qmat)

snps_neutral@meta = cbind(snps_neutral@meta, Qmat)



#M. Barplot in ggplot2:
#define colors and name for populations, in this case according to Monteiro et al. 2021:
colors_pilocarpus = c('#0000ff','#ffc0cb', "#ff0000", '#ffff00') #color for the four genetic cluster
labels_pilocarpus = c("POP D", "POP B", "POP C", "POP A") #name for the four genetic cluster
parcels = unique(snps_neutral@meta$Local_ID)

#create a dataframe (df) for ggplot2, all localities
df = snps_neutral@meta %>%
  dplyr::select(ID, Aggregate, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4) %>%
  dplyr::arrange(Aggregate) %>%
  melt(., id.vars=c("ID", "Aggregate"))
head(df)

p = ggplot(data=df, mapping = aes(x=reorder(factor(ID), Aggregate),
                                y= value*100,
                                fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_pilocarpus,
                    labels= labels_pilocarpus) +
  theme_minimal() +
  xlab("") + ylab("")

pdf(paste0("./Results/Step02/sNMF/Pipe_ggplot2_snmf_AllParcels_", project_name , ".pdf"), onefile =F)
plot(p)
dev.off()

##create a dataframe (df) for ggplot2, by localities
df = snps_neutral@meta %>%
  dplyr::select(ID, Local_ID, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4)
head(df)


for (i in 1:length(parcels)) {
  test = df %>%
    dplyr::filter(Local_ID == parcels[i]) %>%
    dplyr::select(- Local_ID) %>%
    melt(., id.vars=c("ID"))
  
  p = ggplot(data=test, mapping = aes(x=factor(ID),
                                      y= value*100,
                                      fill = factor(variable))) +
    geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
    scale_fill_manual(values= colors_pilocarpus,
                      labels= labels_pilocarpus) +
    theme_minimal()+
    xlab("") + ylab("")
  
  
  pdf(paste0("./Results/Step02/sNMF/Pipe_ggplot2_snmf_", parcels[i] , ".pdf"), onefile =F)
  plot(p)
  dev.off()
}

#N. Save new vcf file with sNMF results and pop ID
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf")) #save vcf
write.csv(snps_neutral@meta, paste0("Results/Step02/Metafiles/ancestry_coef_LEA_snmf_", project_name, ".csv"), quote = F) #save result as table







#### 3.GENETIC STRUCTURE USING DAPC ---- 
#DA Discriminant Analysis focus on between-group variability, while neglecting within-group variation. this method also allows for a probabilistic assignment of individuals to each group, as in Bayesian clustering methods. the method requires the number of variables (alleles) to be less than the number of observations (individuals). This condition is generally not fulfilled in Single Nucleotide Polymorphism (SNP) or re-sequencing datasets. Second, it is hampered by correlations between variables, which necessarily occur in allele frequencies due to the constant-row sum constraint [i.e., compositional data]. Uncorrelated variables will be even more blatant in the presence of linkage disequilibrium
#DAPC relies on data transformation using PCA as a prior step to DA, which ensures that variables submitted to DA are perfectly uncorrelated, and that their number is less than that of analysed individuals (1:10 proportion). Without implying a necessary loss of genetic information, this transformation allows DA to be applied to any genetic data.
#K-means relies on the same model as DA to partition genetic variation into a between-group and a within-group component, and attempts to find groups that minimize the latter. We use Bayesian Information Criterion (BIC) to assess the best supported model, and therefore the number and nature of clusters.

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
variance_pc #PCs that increase the explained variance bellow this threshold will be removed
#calculate number of PCs
index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
index_min

#E. Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust = 40)
#If you see a plateau on the graph, you can choose a number of PCs in the beginning of this plateau. In Pilocarpus there's no plateau. We will test different number of PC's
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain. Complete by the results
index_100 #276 PCs - k = 4 #For me, 100% of variation is more coherent.
index_95  #248 PCs - k = 4
index_70  #147 PCs - k = 6
index_min #88 PCs - k = 7

#Choose the number of PCs and cluster and verify Group (DPCA) of the first 10 individuals and the Size of Groups
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
head(grp$grp, 10)
grp$size
#[1]  65 118  47  47


#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
dapc_df = as.data.frame(grp$Kstat) %>%
  mutate(K = c(1:10)) %>%
  setNames(c("BIC", "K"))
  
p = ggplot(data=dapc_df, mapping = aes(x=K,y= BIC)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  scale_x_continuous(breaks = c(1:10)) +
  theme_bw() +
  theme_genetics +
  xlab("Best K") + ylab("Bayesian Information Criterion (BIC)") #set labels

#check the graphics
p
pdf("./Results/Step02/DAPC/Bestk_DAPC.pdf", onefile = T)
p
dev.off()



#G. ATTENTION, if your dataset is k = 1 go to PCA analysis (Action #4)! Graphs for DAPC need datasets with K >= 2

#H. Choose the best PCs number to recover correctly the clusters. The input should be a Genind object without missing data. Maximum number of PCs is number of individuals -1. Replication by default is 30. Save automatically the graph
pdf("./Results/Step02/DAPC/Best_PCs_Number_DAPC.pdf", onefile = T)
number_PCs = xvalDapc(tab(input, NA.method = "mean"), grp$grp, scale = T, n.pca.max = (nrow(input@tab)-1))
dev.off()

#I. Verify the number of PCs and DA used and summary of DAPC, and percentage explained by the three first dapc
number_PCs$DAPC$n.pca
number_PCs$DAPC$n.da
summary(number_PCs$DAPC)
(number_PCs$DAPC$eig[1]/sum(number_PCs$DAPC$eig))*100
(number_PCs$DAPC$eig[2]/sum(number_PCs$DAPC$eig))*100
(number_PCs$DAPC$eig[3]/sum(number_PCs$DAPC$eig))*100


#J. Verify plots and define the group colors and names:
#plot graph
table(grp$grp)
#[1]  65 118  47  47
snps_neutral@meta[which(grp$grp == 1),] #POP B in snmf
snps_neutral@meta[which(grp$grp == 2),] #POP D
snps_neutral@meta[which(grp$grp == 3),] #POP A
snps_neutral@meta[which(grp$grp == 4),] #POP C
#define colors and name for populations, in this case according to Monteiro et al. 2021:
colors_dapc = c('#ffc0cb', '#0000ff', '#ffff00',  "#ff0000" ) #color for the four genetic cluster
labels_dapc = c("POP B", "POP D", "POP A", "POP C") #name for the four genetic cluster

#ggplot graph:
df_dapc = number_PCs$DAPC$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("IND") %>%
  mutate(POP = as.data.frame(grp$grp)[,1])
  
d_gg = ggplot(df_dapc, aes(x=LD1, y=LD2, fill=POP))+
  geom_point(size=4, pch=21)+
  scale_fill_manual(values= colors_dapc,
                    labels= labels_dapc) +
  theme_bw()+
  theme_genetics +
  xlab("DAPC1 (40.14%)") + ylab("DAPC2 (36.51%)") #set labels

d_gg

pdf("./Results/Step02/DAPC/scatter_DAPC_ggplot.pdf", onefile = T)
d_gg
dev.off()


#classic plot graph
pdf("./Results/Step02/DAPC/scatter_DAPC_classical.pdf", onefile = T)
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = colors_dapc, txt.leg = labels_dapc,
        clabel = FALSE, posi.leg = "center", scree.pca = TRUE, pch=19:20,
        posi.pca = "bottomleft", posi.da = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dev.off()



#STRUCTURE-like graph using posterior probability
pdf("./Results/Step02/DAPC/pipegraph_DAPC.pdf", onefile = T)
compoplot(number_PCs$DAPC, col=col_dapc, legend=TRUE, txt.leg=legend_dapc, cleg=.8)
dev.off()



#K. Add DAPC clusters ID in the VCF file
identical(row.names(as.data.frame(grp$grp)), as.character(snps_neutral@meta$sample_name))
setdiff(row.names(as.data.frame(grp$grp)), as.character(snps_neutral@meta$sample_name))
setdiff(as.character(snps_neutral@meta$sample_name), row.names(as.data.frame(grp$grp)))

snps_neutral@meta$PopID_DAPC = as.character(grp$grp)
head(snps_neutral@meta)


#L. Save new vcf file and pop ID file
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_DAPC.vcf"))
write.csv(snps_neutral@meta, paste0("Results/Step02/Metafiles/DAPC_clusters", project_name, ".csv"))








#### 4.GENETIC STRUCTURE USING TESS3 ---- 
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html

#A. Load the .VCF file with only neutral SNPs:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps_neutral) ##277 individuals and 5268 SNPs.

#B. Create a Genotype matrix
genotypes = GenotypeMatrix(snps) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing;
class(genotypes)
genotypes = replace(genotypes, genotypes == -1, NA)

#C. Create a Matrix with long and lat 
coordinates = snps@meta %>%
  dplyr::select(Longitude, Latitude) %>%
  data.matrix(., rownames.force = NA)
#verify the coords
plot(coordinates, pch = 19, cex = .5, xlab = "Longitude", ylab = "Latitude")


#D. Customize values for run TESS3
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
  save(tess3.ls, file = paste0("./Results/Step02/TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
  
}

#C. Choose the lambda with the minimum cross-validation value:
#create a matrix to save results
cross_value = matrix(NA, max(K)*length(lambda_values), 3)
colnames(cross_value) = c("K", "Lambda", "Crossvalidation")

loop = 0 #Always set loop = 0.
for (i in lambda_values){
  load(file = paste0("./Results/Step02/TESS/Tess3ls_Lambda_", i,"_", project_name, ".RData"))
  for (j in 1:max(K)){
    loop=loop+1
    res = Gettess3res(tess3.ls, K=j)
    cross_value[loop,1] = j
    cross_value[loop,2] = i
    cross_value[loop,3] = min(res$crossvalid.crossentropy)}
}
#save as csv
write.csv(cross_value,paste0("./Results/Step02/TESS/Tess3ls_Crossvalidation_values_", project_name, ".csv"))
#choose best lambda
lambda_tess = as.vector(cross_value[cross_value[,3] == min(cross_value[,3]), ][2])
lambda_tess


#D. Choose best K by graphics:
#load the best lambda project:
load(file = paste0("./Results/Step02/TESS/Tess3ls_Lambda_", lambda_tess,"_", project_name, ".RData"))
#plot results, best K is more distance from dashed line. 4 or 5 clusters.
pdf(paste0("./Results/Step02/TESS/TESS3_RE_PlotK_Lambda_", lambda_tess, ".pdf"), onefile =F)
plot.new()
plot(tess3.ls, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-Entropy score")
segments(x0 = 1, y0 = mean(tess3.ls[[1]]$crossvalid.crossentropy),
         x1 = 10, y1 = mean(tess3.ls[[10]]$crossvalid.crossentropy), lty =2, lwd = 2)
dev.off()


#E. Choose best K by assessement greater than 80%, following Bastos et al. 2016 doi: 10.1371/journal.pone.0145165:
#create the dataframe:
results_acess = as.data.frame(matrix(NA, 10, 3))
colnames(results_acess) = c("K", "H_Assessment", "Percentage")
for (l in K){
  res_tess = Gettess3res(tess3.ls, K=l)
  z = res_tess$Q > 0.80
  z = table(z)["TRUE"]
  pc = (z/length(snps_neutral@sample_id))*100
  results_acess[l,] = c(l, z, pc)
}
#plot results. After K=2 which K has the assessement > 0.85 to   
pdf(paste0("./Results/Step02/TESS/TESS3_Assessment_Lambda_", lambda_tess, ".pdf"), onefile =F)
ggplot(data = results_acess, mapping = aes(x= K, y = Percentage)) + 
  geom_line() +
  geom_point(size= 2) +
  xlab("Number of Clusters (K)") + ylab("Assessement (%)") +
  scale_x_continuous(breaks = c(1:10)) +
  theme_bw() +
  theme_genetics +
  geom_label_repel(mapping = aes(x= K, y = Percentage, label=round(Percentage,1)))
dev.off()



#F.Choose the Best K and create the Q-matrix
#Best K
optimal_K = 5
res_tess3 = Gettess3res(tess3.ls, K=optimal_K)



#G. Add admixture coefficient and replace the population ID to vcf file
Qmat =
  res_tess3$Q %>%
  as.data.frame() %>%
  setNames(c("Adx_Coeff_1", "Adx_Coeff_2", "Adx_Coeff_3","Adx_Coeff_4", "Adx_Coeff_5")) %>%
  mutate(PopID_TESS3 = apply(., 1, which.max))
head(Qmat)

snps_neutral@meta = cbind(snps_neutral@meta, Qmat)
head(snps_neutral@meta)



#H. Structure-like barplot for the Q-matrix
#define colors and name for populations, in this case according to Monteiro et al. 2021:
snps_neutral@meta[snps_neutral@meta$PopID_TESS3 == 1,] #POP B in snmf
snps_neutral@meta[snps_neutral@meta$PopID_TESS3 == 2,] #Alemão 
snps_neutral@meta[snps_neutral@meta$PopID_TESS3 == 3,] #POP A
snps_neutral@meta[snps_neutral@meta$PopID_TESS3 == 4,] #POP D
snps_neutral@meta[snps_neutral@meta$PopID_TESS3 == 5,] #POP C

colors_tess = c('#ffc0cb', 'darkslategrey', '#ffff00', '#0000ff', "#ff0000") #color for the genetic cluster
labels_tess = c("POP B", "POP E", "POP A", "POP D", "POP C") #name for the genetic cluster
parcels = unique(snps_neutral@meta$Local_ID)


#create a dataframe (df) for ggplot2, all localities
df = snps_neutral@meta %>%
  dplyr::select(ID, Aggregate, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4, Adx_Coeff_5) %>%
  dplyr::arrange(Aggregate) %>%
  melt(., id.vars=c("ID", "Aggregate"))
head(df)

p = ggplot(data=df, mapping = aes(x=reorder(factor(ID), Aggregate),
                                  y= value*100,
                                  fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_tess,
                    labels= labels_tess) +
  theme_minimal() +
  xlab("") + ylab("")
#check
p

pdf(paste0("./Results/Step02/TESS/Pipe_ggplot2_snmf_AllParcels_", project_name, ".pdf"), onefile =F)
plot(p)
dev.off()

##create a dataframe (df) for ggplot2, by localities
df = snps_neutral@meta %>%
  dplyr::select(ID, Local_ID, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3, Adx_Coeff_4,Adx_Coeff_5)
head(df)


for (i in 1:length(parcels)) {
  test = df %>%
    dplyr::filter(Local_ID == parcels[i]) %>%
    dplyr::select(- Local_ID) %>%
    melt(., id.vars=c("ID"))
  
  p = ggplot(data=test, mapping = aes(x=factor(ID),
                                      y= value*100,
                                      fill = factor(variable))) +
    geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
    scale_fill_manual(values= colors_tess,
                      labels= labels_tess) +
    theme_minimal()+
    xlab("") + ylab("")
  
  
  pdf(paste0("./Results/Step02/TESS/Pipe_ggplot2_snmf_", parcels[i] , ".pdf"), onefile =F)
  plot(p)
  dev.off()
}



#I. Spatial interpolation of ancestry coefficient
pdf("./Results/Step02/TESS/TESS_MAP.pdf", onefile =F)
plot.new()
plot(res_tess3$Q, coordinates, method = "map.max", interpolation.model = FieldsKrigModel(10),
     main = "Ancestry coefficients",cex.main=1.5,
     xlab = "Longitude", ylab = "Latitude",cex.lab=1.3,
     resolution = c(700,700), cex = 1.5,
     col = colors_tess)
dev.off()



#J. Save new vcf file with TESS results and pop ID
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_TESS.vcf")) #save vcf
write.csv(snps_neutral@meta, paste0("Results/Step02/Metafiles/ancestry_coef_TESS_", project_name, ".csv"), quote = F) #save result as table







#### 5.GENETIC STRUCTURE USING PCA ---- 
#PCA is its ability to identify genetic structures in very large datasets within negligible computational time, and the absence of any assumption about the underlying population genetic model. PCA aims to summarize the overall variability among individuals, which includes both the divergence between groups (i.e., structured genetic variability), and the variation occurring within groups (‘random’ genetic variability).

#A. Load neutral .vcf file after removing outlier SNPs:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps_neutral) ##277 individuals and 5268 SNPs.

vcf = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral.vcf"), verbose = FALSE)

#B. Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input

#C. Perform a PCA to choose the number of PC in the DAPC:
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf=3)


#D. % of contribution for three first PCs
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")

#E. Min and max values for axes
pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))

#F. save results as .CSV
#coordinates in the PCA by axis
write.csv(pca_input$li, file = "Results/Step02/PCA/Axis_coord.csv")
#% of PC variation
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "Results/Step02/PCA/contribution_pc_eig.csv")


#G. PCA graphs
df_pca = pca_input$li %>%
  mutate(Local = snps_neutral@meta$Local_ID)
head(df_pca)  

pca_graph = ggplot(df_pca, aes(x=Axis1, y=Axis2))+
  geom_point(size=4, pch=21)+
  theme_bw()+
  theme_genetics +
  xlab("PC1 (3.35%)") + ylab("PC2 (3.03%)") + #set labels
  geom_text_repel(aes(label = Local), size = 1.5, max.overlaps = 50) #set labels
#check
pca_graph

#save
pdf(paste0("./Results/Step02/PCA/PCA_", project_name,".pdf"))
pca_graph
dev.off()


#H. TRACY-WIDOM TEST FOR EIGENVALUES
#Tracy CA and Widom H. (1994). Level spacing distributions and the bessel kernel. Commun Math Phys. 161 :289--309. Patterson N, Price AL and Reich D. (2006). Population structure and eigenanalysis. PLoS Genet. 2 :20.

#Load the VCF file after DAPC step
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID = T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.
#verify the metadata 
snps_neutral@meta

#convert VCF to LFMM
lfmm_pca = Lfmm(snps_neutral, output.file=paste0("./Results/Step02/PCA/PCA_LFMM_", project_name, "_filtered_neutral.lfmm"))

#create a pcaProject object:
pc = LEA::pca(paste0("./Results/Step02/PCA/PCA_LFMM_", project_name, "_filtered_neutral.lfmm"), scale = TRUE, center=TRUE)

#perform Tracy-Widom tests
options(scipen = 999)
tw = tracy.widom(pc)
#display the p-values for the Tracy-Widom tests. 
tw
#24 PCs are <0.05. K = PC -1. So my K = 23 

#Save results
write.csv(tw,paste0("./Results/Step02/PCA/PCA_Tracy-Widom_", project_name, ".csv"))

#I. PLOT PC1 and PC2 VALUES IN A QQPLOT TO TEST NORMALITY
pdf("./Results/Step02/PCA/PCA_Normality_PC1_qqplot.pdf", onefile =F)
plot.new()
qqnorm(pca_input$li$Axis1, pch = 1, frame = FALSE)
qqline(pca_input$li$Axis1, col = "black", lwd = 2)
dev.off()

pdf("./Results/Step02/PCA/PCA_Normality_PC2_qqplot.pdf", onefile =F)
plot.new()
qqnorm(pca_input$li$Axis2, pch = 1, frame = FALSE)
qqline(pca_input$li$Axis2, col = "black", lwd = 2)
dev.off()





#### 6.COMPARING METHODS ---- 
#DO THE CLUSTERS (K) FROM DIFFERENT ANALYSES HAVE THE SAME COMPOSITION? ONLY IF BOTH ANALYSES WERE K >= 2:

#A. Position of samples by population by sNMF approach:
snps_snmf = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_LEA.vcf"), overwriteID=T)
VCFsummary(snps_snmf) ##277 individuals and 5268 SNPs.
for (i in 1:length(unique(snps_snmf@meta$PopID_snmf))){
  pop = which(snps_snmf@meta$PopID_snmf == i)
  assign(paste0("pop_snmf_", i), pop)
}


#B. Position of samples by population by DAPC approach:
snps_dapc = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_DAPC.vcf"), overwriteID=T)
VCFsummary(snps_dapc) ##277 individuals and 5268 SNPs.
for (i in 1:length(unique(snps_dapc@meta$PopID_DAPC))){
  pop = which(snps_dapc@meta$PopID_DAPC == i)
  assign(paste0("pop_DAPC_", i), pop)
}

#C. Position of samples by population by TESS approach:
snps_tess = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_TESS.vcf"), overwriteID=T)
VCFsummary(snps_tess) ##277 individuals and 5268 SNPs.
for (i in 1:length(unique(snps_tess@meta$PopID_TESS3))){
  pop = which(snps_tess@meta$PopID_TESS3 == i)
  assign(paste0("pop_TESS_", i), pop)
}

#D. PCA approach do not have a assignment tools but you can analyze the clusters by the pdf saved in action #5. Our case there is 4 populations with same division than snmf and DAPC.


#E. Compare the sample size. Add the number of pop that you have:
pop_size_snmf = c(length(pop_snmf_1), length(pop_snmf_2), length(pop_snmf_3), length(pop_snmf_4))
pop_size_tess = c(length(pop_TESS_1), length(pop_TESS_2), length(pop_TESS_3), length(pop_TESS_4), length(pop_TESS_5))

pop_size_snmf # 118  65  47  47
grp$size #65 118  47  47
pop_size_tess #65  15  47 103  47 


#G. Check if identical samples in each group. If they are corrected, function will return "TRUE"
#POPs with the same sample size #BLUE
identical(pop_snmf_1, pop_DAPC_2, sort(c(pop_TESS_2, pop_TESS_4))) #Alemão, Aggregate 21, is a isolated cluster

#POPs with the same sample size in three or more analyses
identical(pop_snmf_2, pop_DAPC_1, pop_TESS_1) #TRUE

#POPs with the same sample size
identical(pop_snmf_3, pop_DAPC_4, pop_TESS_3) #TRUE

#POPs with the same sample size
identical(pop_snmf_4, pop_DAPC_3, pop_TESS_5) #TRUE

# If the results are false, check the difference between files, which samples are present in a file and not in other. Invert the order of population. If the 1st pop is inside the 2nd pop function will return "character (0)"
setdiff(pop_snmf_2, pop_DAPC_3)
setdiff(pop_DAPC_3, pop_snmf_2)


#H. Choose the number of clusters (K) to the following step! To our data is 4!!!
#Load neutral .vcf file after removing outlier SNPs:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps_neutral) ##277 individuals and 5268 SNPs.

#add information on genetic cluster in the metafile. Replace identical cols:
meta = cbind(snps_neutral@meta, snps_snmf@meta[,12:16], snps_dapc@meta[,12], snps_tess@meta[,12:17])
names(meta)[17:22] = c("PopID_DAPC", "Adx_Coeff_1_Tess", "Adx_Coeff_2_Tess", "Adx_Coeff_3_Tess", "Adx_Coeff_4_Tess", "Adx_Coeff_5_Tess")
head(meta)

snps_neutral@meta = meta
head(snps_neutral@meta)

#add information on genetic cluster of Monteiro et al. as well:
snps_neutral@meta[which(snps_neutral@meta$PopID_snmf == 1),] #POP D - Blue
snps_neutral@meta[which(snps_neutral@meta$PopID_snmf == 2),] #POP B - Pink
snps_neutral@meta[which(snps_neutral@meta$PopID_snmf == 3),] #POP C - Red
snps_neutral@meta[which(snps_neutral@meta$PopID_snmf == 4),] #POP A - Yellow

POP_ID_Monteiro = snps_neutral@meta$PopID_snmf %>%
  as.data.frame() %>%
 mutate(Monteiro = replace(., .=="1", "D") %>%
  replace(., .=="2", "B") %>%
  replace(., .=="3", "C") %>%
  replace(., .=="4", "A"))


snps_neutral@meta = cbind(snps_neutral@meta, POP_ID_Monteiro$Monteiro)
names(snps_neutral@meta)[24] = c("POP_ID")
head(snps_neutral@meta)


#J. Save new vcf file with TESS results and pop ID
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf")) #save vcf
write.csv(snps_neutral@meta, paste0("Results/Step02/Metafiles/cluster_information_", project_name, ".csv"), quote = F) #save result as table


##END
