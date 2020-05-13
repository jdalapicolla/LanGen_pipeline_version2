############################
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

##########
fgraph <- function(obj){
  # use multispati summary
  sum.obj <- summary(obj)
  # compute Imin and Imax
  L <- listw2mat(eval(as.list(obj$call)$listw))
  Imin <- min(eigen(0.5*(L+t(L)))$values)
  Imax <- max(eigen(0.5*(L+t(L)))$values)
  I0 <- -1/(nrow(obj$li)-1)
  # create labels
  labels <- lapply(1:length(obj$eig),function(i) bquote(lambda[.(i)]))
  # draw the plot
  xmax <- eval(as.list(obj$call)$dudi)$eig[1]*1.1
  par(las=1)
  var <- sum.obj[,2]
  moran <- sum.obj[,3]
  plot(x=var,y=moran,type='n',xlab='Inertia',ylab="Spatial autocorrelation (I)",
       xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n')
  text(x=var,y=moran,do.call(expression,labels))
  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),
             ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)
  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)
  title("Spatial and inertia components of the eigenvalues")
}

############################
Diference_Final <- function(All_SNPs, P1, P2, P3){ 
  Diference_1 <- intersect (All_SNPs@site_id, P1@site_id)
  Diference_2 <- intersect (All_SNPs@site_id, P2@site_id)
  Diference_3 <- intersect (All_SNPs@site_id, P3@site_id)
  
  Diference_All_1 <- intersect (Diference_1, Diference_2) 
  Diference_Final <- intersect (Diference_All_1, Diference_3)
  
  return(Diference_Final)
}

############################
create_dir = function(dir_names){
  for (i in 1:length(dir_names)){
    if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
    else (print (paste0(dir_names[i], " has already been created. Be careful you can overwrite files")))}
  
}

##################################
convert_vcf2genind = function(dataset){
    snps = dataset
    if (class(snps) == "vcfR"){
       genlight = vcfR2genlight(snps)
       genind = gl2gi(genlight)
       return(genind)
} else (print (" Dataset is not a 'vcfR' object. Please use the package 'vcfR' to laod the file")) }

#############################################
pca_missing = function(datasets, names, threshold =20, missing_values, samples_ID = row.names(input@tab) , color_regular ='blue', color_missing = 'red', plot_all_names = FALSE) {
  for (i in 1:length(datasets)){
    input =  datasets[[i]]
    if (class(input) != "genind") {print (" Dataset is not a 'genind' object. Please use the function 'vcf2genind' first")}
    else{
      input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
      pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = 3)
  
      #G. CALCULATE THE % CONTRIBUTION OF EACH PC
      PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
      PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
      PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")
  
      #I. MIN AND MAX VALUES FOR EACH AXIS
      pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
      pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
      pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))
  
      #J. SAVE THE RESULTS AS .CSV
      #coordenates in the PCA by axis
      write.csv(pca_input$li, file = paste0("Results_Filters/PCA/Axis_coord_test_", project_name, "_", names[i], "_", threshold,"_PC.csv"))
      #% of PC variation
      perc_pca = as.data.frame(pca_input$eig)
      soma = sum (perc_pca)
      perc_pca[,2] = (perc_pca/soma)*100
      colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
      write.csv(perc_pca, file = paste0("Results_Filters/PCA/contribution_pc_eig_test_", project_name,"_", names[i], "_", threshold, "_PC.csv"))
  
      #K. SAVE/VIEW THE RESULTS AS .PDF TO LEA CLUSTERS
      #for save the results
      pdf(paste0("./Results_Filters/PCA/PCA_test_", project_name,"_", names[i], "_", threshold, "_PC.pdf")) #change the name of pdf
      plot.new()
  
      #after choose, continue with these commands
      par(mar=c(5,5,1,1)); #margins of a single graph (unit = lines width)
      plot.window(xlim=pc1_range, ylim=pc2_range); #size of the graph
  
      #number of samples 
      samples = length(samples_ID)
      pop.names = samples_ID
  
      #plot all points
      points(pca_input$li[1:samples,1], pca_input$li[1:samples,2], col = 'black', bg = color_regular, cex = 1.5, pch=21) 
  
      #plot points in red with higher missing data amount
      missing_pc = which(missing_values > threshold)
      points(pca_input$li[missing_pc,1], pca_input$li[missing_pc,2], col = 'black', bg = color_missing, cex = 1.5, pch=21)
  
      #draw the axis; put the limits similar to PC Min and Max values
      axis(1, at=seq((round(pc1_range[1],-1)-10), (round(pc1_range[2],-1)+10), by=10), cex.axis=1.15)
      axis(2, at=seq((round(pc2_range[1],-1)-10), (round(pc2_range[2],-1)+10), by=10), cex.axis=1.15, las=1)
  
      #complete the name of each axis
      mtext(side=1, text=PC1, line=2.5, cex=1)
      mtext(side=2, text=PC2, line=2.8, cex=1)
      
      if(plot_all_names == TRUE){
        #labeling all individuals on the plot
        pop.names = samples_ID
        text(pca_input$li[1:samples,1], pca_input$li[1:samples,2], labels = pop.names[1:samples], pos = 3, cex = 0.7)}
  
      if(plot_all_names == FALSE){
        #labeling individuals with more than a threshold % of missing data in red letters
        text(pca_input$li[missing_pc,1], pca_input$li[missing_pc,2], labels = pop.names[missing_pc], pos = 3, col = color_missing, cex = 0.7)}
  
  #finishing the graph
  dev.off()
}
}
}

############################################
GenDiv_IND <- function(snps){
  HE.ESP <- Query(snps, type="het")
  HE.ESP$HE.ESP <- (HE.ESP$N_SITES-HE.ESP$E.HOM.)/(HE.ESP$N_SITES) ## Expected heterozygosity (HE.ESP)
  #head(HE.ESP)
  he <- (HE.ESP$HE.ESP) ## Expected heterozygosity Mean (HE.ESP_Mean)
  ho <- (1-(HE.ESP$O.HOM./HE.ESP$N_SITES)) ## Observed heterozygosity (Ho.ESP_Mean)
  Hom <-  (HE.ESP$O.HOM./HE.ESP$N_SITES)
  f <- (HE.ESP$F) ## Fis (inbreeding)
  PI <- Query(snps, type="site-pi")
  #head(PI)
  #pi <- (PI$PI) ## Nucleotide divergency
  print(paste0("OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS")) ## DIVERGENCE STATISTICS
  print(paste0("Observed homozygosity (Ho_O)"))
  print(paste0("Observed heterozygosity (He_O)"))
  print(paste0("Expected heterozygosity (He_E)"))
  print(paste0("Coefficient of inbreeding (F)"))
  #print(paste0("Measures nucleotide divergency on a per-site basis (PI)"))
  GD <- data.frame(Ho_O=Hom, He_O=ho, He_E=he, FIS=f)
  return(GD)
}

#############################
GenDiv <- function(vcf){
  HE <- Query(vcf, type="het")
  
  HE$HO <- (HE$N_SITES-HE$O.HOM.)/(HE$N_SITES) ## Observed heterozygosity (HO)
  error_ho <- qt(0.975,df=length(HE$HO)-1)*sd(HE$HO)/sqrt(length(HE$HO))
  ho <- mean(HE$HO)
  left_ho <- ho-error_ho
  right_ho <- ho+error_ho
  HO.df <- data.frame(He_O=ho, low_He_O=left_ho, up_He_O=right_ho)
  
  HE$HE <- (HE$N_SITES-HE$E.HOM.)/(HE$N_SITES) ## Expected heterozygosity (HE)
  error_he <- qt(0.975,df=length(HE$HE)-1)*sd(HE$HE)/sqrt(length(HE$HE))
  he <- mean(HE$HE)
  left_he <- he-error_he
  right_he <- he+error_he
  HE.df <- data.frame(He_E=he, low_He_E=left_he, up_He_E=right_he)
  
  error_f <- qt(0.975,df=length(HE$F)-1)*sd(HE$F)/sqrt(length(HE$F))
  f <- mean(HE$F)
  left_f <- f-error_f
  right_f <- f+error_f
  F.df <- data.frame(F=f, low_F=left_f, up_F=right_f)
  
  PI <- Query(vcf, type="site-pi")
  error_pi <- qt(0.975,df=length(PI$PI)-1)*sd(PI$PI)/sqrt(length(PI$PI))
  pi <- mean(PI$PI)
  left_pi <- pi-error_pi
  right_pi <- pi+error_pi
  PI.df <- data.frame(PI=pi, low_PI=left_pi, up_PI=right_pi)
  
  print(paste0("OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS")) ## DIVERGENCE STATISTICS
  print(paste0("Observed heterozygosity (He_O)"))
  print(paste0("Expected heterozygosity (He_E)"))
  print(paste0("Coefficient of inbreeding (F)"))
  print(paste0("Measures nucleotide divergency on a per-site basis (PI)"))
  print(paste0("The 95% confidence interval - Lower limit (low_)"))
  print(paste0("The 95% confidence interval - Upper limit (up_)"))
  RES <- cbind(HO.df,  HE.df,  F.df,  PI.df)
  return(RES)
}

#############################
PlotK <- function(snmfproject){
  
  Mk1 <- mean(cross.entropy(snmfproject, K=1))
  SDk1 <- sd(cross.entropy(snmfproject, K=1))
  Mk2 <- mean(cross.entropy(snmfproject, K=2))
  SDk2 <- sd(cross.entropy(snmfproject, K=2))
  Mk3 <- mean(cross.entropy(snmfproject, K=3))
  SDk3 <- sd(cross.entropy(snmfproject, K=3))
  Mk4 <- mean(cross.entropy(snmfproject, K=4))
  SDk4 <- sd(cross.entropy(snmfproject, K=4))
  Mk5 <- mean(cross.entropy(snmfproject, K=5))
  SDk5 <- sd(cross.entropy(snmfproject, K=5))
  Mk6 <- mean(cross.entropy(snmfproject, K=6))
  SDk6 <- sd(cross.entropy(snmfproject, K=6))
  Mk7 <- mean(cross.entropy(snmfproject, K=7))
  SDk7 <- sd(cross.entropy(snmfproject, K=7))
  Mk8 <- mean(cross.entropy(snmfproject, K=8))
  SDk8 <- sd(cross.entropy(snmfproject, K=8))
  Mk9 <- mean(cross.entropy(snmfproject, K=9))
  SDk9 <- sd(cross.entropy(snmfproject, K=9))
  Mk10 <- mean(cross.entropy(snmfproject, K=10))
  SDk10 <- sd(cross.entropy(snmfproject, K=10))
  
  CE <- data.frame(K=c(1:10), Mean = c(Mk1,Mk2,Mk3,Mk4,Mk5,Mk6,Mk7,Mk8,Mk9,Mk10),
                   SD = c(SDk1,SDk2,SDk3,SDk4,SDk5,SDk6,SDk7,SDk8,SDk9,SDk10))
  
  library(ggplot2)
  
  ggplot(CE, aes(x=K, y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2)+
    geom_line() + 
    geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of ancestral populations") + ylab("Cross-entropy")+
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=15, face ="bold" , color = "black"), axis.title.x = element_text(size=15, face="bold", color="black"),axis.title.y = element_text(size=15, face="bold", color="black")) +
    scale_x_continuous(breaks = seq(0,10, by=2))
}

#############################
Best.run <- function(nrep, optimalK, p1, p2, p3, p4, p5, p6){
  ce1 = LEA::cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = LEA::cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = LEA::cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = LEA::cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  ce5 = LEA::cross.entropy(p5, K = optimalK) # get the cross-entropy of each run for optimal K
  ce6 = LEA::cross.entropy(p6, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4, ce5, ce6)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep), rep(5, nrep), rep(6, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, " run = ", Best_run)) 
}

##################################
barplotK <- function(Qfile, Pop, Run_B){
  Q.matrix_1 <-LEA::Q(Qfile, K = Pop, run = Run_B)
  Q.matrix <- as.qmatrix(Q.matrix_1)
  barplot(Q.matrix, xlab = "Sampled individuals",
          ylab = "Ancestry coefficients",
          main = "", cex.axis = 1, cex.lab = 1.5)
}

#################################
fst = function(project, run = 1, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t}
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}

#################################
GIF <- function(project, run, K, fst.values){
  n = dim(Q(project, K, run))[1]
  fst.values[fst.values<0] = 0.000001
  z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
  lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
  print(lambda)
  return(lambda)
}

#################################
candidates <- function(alpha=0.05, adj.p.values){ ## alpha is significance level
  L = length(adj.p.values) ## Number of loci
  w = which(sort(adj.p.values) < alpha * (1:L) / L)
  candidates = order(adj.p.values)[w]
  print(length(candidates))
  return(candidates)
}

#################################
ManPlot <- function(adj.p.values, candidates, title){
  plot(-log10(adj.p.values), main=title, xlab = "Locus", cex = .7, col = "grey")
  points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")
}

#################################
summary_stat = function(datasets, names){
  
  for (i in 1:length(datasets)){
  input2 =  datasets[[i]]
  
  #for global analysis. For Phi_st No missing data allowed!!!!
  FST_stat = diff_stats(input2, phi_st = F)
  write.csv(FST_stat$global, file = paste0("./Results_Diversity/FST_global_", project_name,"_", names[i], ".csv"))
  
  #plotting the results
  per.locus = melt(FST_stat$per.locus, varnames = c("Locus", "Statistic"))
  stats = c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
  glob = data.frame(Statistic = stats, value = FST_stat$global)
  head(per.locus)
  head(glob)
  
  #create the theme for ggplot graphs
  theme_pca = theme(legend.text = element_text(face = "italic",
                                               colour="black",
                                               family = "Helvetica",
                                               size = rel(1)), 
                    axis.title = element_text(colour="black",
                                              family = "Helvetica",
                                              size = rel(1.2)), 
                    axis.text = element_text(family = "Helvetica",
                                             colour = "black",
                                             size = rel(1)), 
                    axis.line = element_line(size = 1,colour = "black"), 
                    axis.ticks = element_line(colour="black",size = rel(1)),
                    
                    panel.grid.minor = element_blank(), 
                    panel.background = element_rect(fill = "whitesmoke"), 
                    panel.grid.major = element_line(colour="black",size = rel(0.2), linetype = "dotted"),
                    legend.key = element_blank(), 
                    legend.title = element_text(colour = "black",
                                                size = rel(1.5),
                                                family = "Helvetica"), 
                    plot.title = element_text(colour = "black",
                                              face = "bold",
                                              hjust = 0.5, #alingment
                                              size = rel(1.7),
                                              family = "Helvetica"))
  
  #plot
  pp = ggplot(per.locus, aes(x = Statistic, y = value)) +
    geom_boxplot() +
    geom_point() +
    geom_point(size = rel(3), color = "red", data = glob) +
    ggtitle("Estimates of Population Differentiation") +
    theme_pca
  
  
  pdf(paste0("Results_Diversity/Estimates_Population_Differentiation_", names[i], ".pdf"), onefile = F)
  print(plot(pp))
  dev.off()
  
  
  #Hs and Ht: heterozygosity expected for this population with and without the sub-populations defined in the input data respectively. Global is a mean.
  #Gst_Nei = Gst in per locus or Gst_est in global.
  #Hedrick's G"st: is Gprime_st both areas.
  #Jost's D = D per loci and two global estimates. 'D_het' uses the averages of Hs and Ht across all loci while 'D_mean' takes the harmonic mean of all loci. It will be NA with negative values.
  #Because all of these statistics are estimated from estimators of HS and HT, it’s possible to get negative values for each of these differentiation measures. Populations can’t be negatively differentiated, so you should think of these as estimates of a number close to zero (it’s up to you and your reviewers to decide if you report the negative numbers of just zeros).
  
  #by populations. Number meanings are in 5C.
  D = pairwise_D(input2, linearized = FALSE)
  D = as.matrix(D)
  write.csv(D, file = paste0("Results_Diversity/D_groups_", project_name ,"_", names[i],".csv"))
  
  GST_N = pairwise_Gst_Nei(input2, linearized = FALSE)
  GST_N = as.matrix(GST_N)
  write.csv(GST_N, file = paste0("Results_Diversity/GST_N_groups_",project_name ,"_", names[i], ".csv"))
  
  GST_H = pairwise_Gst_Hedrick(input2, linearized = FALSE)
  GST_H = as.matrix(GST_H)
  write.csv(GST_H, file = paste0("Results_Diversity/GST_H_groups_", project_name ,"_", names[i], ".csv"))
  
  #by populations and global to other indexes. See the help for poppr() for acronyms' meanings:
  sum_stat = poppr(input2)
  write.csv(sum_stat, file = paste0("Results_Diversity/Sum_Stat_groups_", project_name ,"_", names[i], ".csv"))
}
}

################################
shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

############################
Plotntile <- function (path.mcmc, burnin, printit = FALSE, file) 
{
  fileparam <- paste(path.mcmc, "parameters.txt", sep = "/")
  param <- as.matrix(read.table(fileparam))
  thinning <- as.numeric(param[param[, 1] == "thinning", 3])
  filenpp <- paste(path.mcmc, "nuclei.numbers.txt", sep = "")
  npp <- scan(filenpp)
  if (burnin > 0) {
    sub <- -(1:burnin)
  }
  else {
    sub <- 1:length(npp)
  }
  if (printit == TRUE) {
    postscript(file)
    par(mfrow = c(1, 2))
    plot((1:length(npp)) * thinning, npp, type = "l", ylab = "Number of Tiles", 
         xlab = paste("Index of MCMC iteration", " \n Whole chain", 
                      sep = ""))
    hist(npp[sub], plot = TRUE, prob = TRUE, xlab = paste("Number of tiles along the chain \n(after a burnin of ", 
                                                          burnin, "x", thinning, " it.)", sep = ""), main = "Number of tiles along the chain  \n  after burnin")
    dev.off()
  }
  else {
    par(mfrow = c(1, 2))
    plot((1:length(npp)) * thinning, npp, type = "l", ylab = "Number of Tiles", 
         xlab = paste("Index of MCMC iteration", "\n Whole chain", 
                      sep = ""))
    hist(npp[sub], plot = TRUE, prob = TRUE, xlab = paste("Number of tiles along the chain \n(after a burnin of ", 
                                                          burnin, "x", thinning, " it.)", sep = ""), main = "Number of tiles along the chain  \n  after burnin")
  }
}
##################################
################################
mDistoLoess <- function(mat1,mat2,nperm=999){
  ltd <- mat1[lower.tri(mat1)]
  ltk <- mat2[lower.tri(mat2)]
  disto <- matrix(nrow = nperm, ncol = length(ltk))
  for(n in 1:nperm){
    perm <- sample(1:nrow(mat2))
    mat2_p <- mat2[perm,perm]
    lt <- mat2_p[lower.tri(mat2_p)]
    disto[n,] <- predict(loess(lt~ltd))
    if(n %% 100 == 0) print(n)
  }
  obs_fit <- loess(ltk~ltd)
  obs <- predict(obs_fit)
  CI <- apply(disto, 2, quantile, probs = c(0.025,0.975))
  attr(disto, "obs") <- obs
  attr(disto, "CI") <- CI
  return(list(disto=disto))
}