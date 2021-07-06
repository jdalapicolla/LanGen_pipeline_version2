### LANDSCAPE GENOMICS PIPELINE IN R ###
 
REFERENCES AND MORE INFORMATION:

#A. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
 
#B. THIS PIPELINE WAS DESIGNED IN UBUNTU 20.04 LTS, USING RSTUDIO 1.4.1717 AND R 4.1.0

#C. Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé



### PART 1 - GENETIC STRUCTURE AND GENETIC DIVERSITY ###

This suit of scripts were developed based on LanGen_pipeline from @rojaff/LanGen_pipeline

Scripts in R for different analyses of Landscape Genomics using SNPS, VCF, or Genind objects. These scripts perform the following steps:

<b>1.</b> SNPs dataset filtering using PCadapt, sNMF, and TESS3 to outlier detection.

<b>2.</b> Assessment of population genetic structure using sNMF, DAPC, TESS3, and PCA.

<b>3.</b> Estimate Tajima's D, genetic diversity and genetic distance inter and intra populations, and preparing the input for NeEstimator (Effective Population Sizes).

## PART 1 - GENETIC STRUCTURE AND GENETIC DIVERSITY https://github.com/jdalapicolla/LanGen_pipeline_version2
## PART 2 - ISOLATION BY DISTANCE AND FINE-SCALE SPATIAL GENETIC STRUCTURE https://github.com/jdalapicolla/IBD_models.R
## PART 3 - ISOLATION BY RESISTANCE USING MLPE MODELS - https://github.com/jdalapicolla/MLPE.R
## PART 4 - LOCAL ADAPTATION ANALYSES - IDENTIFICATION OF CANDIDATES LOCI UNDER SELECTION - https://github.com/jdalapicolla/LOCAL_ADAPTATION.R
    
