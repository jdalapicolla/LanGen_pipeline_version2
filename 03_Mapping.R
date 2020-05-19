###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
###########################   STEP 03: MAPPING   ##############################

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
#A. THE FILE ".VCF" AFTER STRUCTURE ESTIMATION, STEP 2, WITH GEOGRAPHICAL INFORMATIONS AND ADMIXTURE COEFFICIENTS BY INDIVIDUALS IN THE METAFILE.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline
#C. RASTERS AND SHAPEFILES FOR THE STUDY AREA TO CREATE THE MAPS AND PUT ALL IN A SINGLE FOLDER CALLED "maps", WITH SUBFOLDERS TO "shapefiles" AND "rasters".

##3. GOALS FOR THIS STEP:
#A. CREATE MAPS OF POPULATION ADMIXTURE COEFFICIENTS.

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
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("rgdal" %in% rownames(installed.packages()) == FALSE){install.packages("rgdal")
} else {print (paste0("'rgdal' has already been installed in library"))}
if("rnaturalearth" %in% rownames(installed.packages()) == FALSE){install.packages("rnaturalearth")
} else {print (paste0("'rnaturalearth' has already been installed in library"))}
if("GISTools" %in% rownames(installed.packages()) == FALSE){install.packages("GISTools")
} else {print (paste0("'GISTools' has already been installed in library"))}
if("mapplots" %in% rownames(installed.packages()) == FALSE){install.packages("mapplots")
} else {print (paste0("'mapplots' has already been installed in library"))}
if("tess3r" %in% rownames(installed.packages()) == FALSE){devtools::install_github("bcm-uga/TESS3_encho_sen")
} else {print (paste0("'tess3r' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, raster, rgdal, rnaturalearth, GISTools, tess3r, mapplots)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Maps"))


#------------------------------------------------------------------------------
#                    1. Loading Files and Configuring Maps 
#------------------------------------------------------------------------------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "pilocarpus"


###1.2. LOAD VCF FILES, GEOGRAPHICAL INFORMATIONS AND ADMIXTURE COEFFICIENTS: 
#A. Load neutral .vcf file with geographical information and admixture and :
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.

#B. Create a data frame with the information on metafile
df = as.data.frame(snps_neutral@meta)
head(df)
tail(df)

#C. Create vectors with samples' position by population in sNMF approach:
for (i in 1:length(unique(df$PopID_snmf))){
  pop = which(df$PopID_snmf == i)
  assign(paste0("pop_snmf_", i), pop)
}

#D. Create vectors with samples' position by population in DAPC approach:
for (i in 1:length(unique(df$PopID_DAPC))){
  pop = which(df$PopID_DAPC == i)
  assign(paste0("pop_DAPC_", i), pop)
}

#E. Create vectors with samples' position by population in TESS approach:
for (i in 1:length(unique(df$PopID_tess))){
  pop = which(df$PopID_tess == i)
  assign(paste0("pop_TESS_", i), pop)
}


###1.3. DEFINE A EXTENSION FOR THE MAPS
#A. Calculate the range for the geographical coordinates:
long_extension = c(min(df$longitude), max(df$longitude))
lat_extension = c(min(df$latitude), max(df$latitude))
#verify the range of coordenates
long_extension
lat_extension

#B. Plot points for samples 
plot(snps_neutral@meta$longitude, snps_neutral@meta$latitude)

#C. Choose a extension for maps that cover all points. In this case:
ext = extent(-50.7, -49.7, -6.5, -5.8)

#D. Define the axes from maps
axis.x = c(-50.7, -49.7)
axis.y = c(-6.5, -5.8)


###1.4. CREATING RASTER FOR THE BASE OF MAPS - HILLSHADE, SLOPE AND LANDCOVER
#A. Load DEM raster (elevation raster):
DEM = raster("maps/rasters/DEM_Carajas.tif")
#plot the map
plot(DEM, col = gray.colors(20, start = 0, end = 1))
#crop according to study area
DEM_base = crop(DEM, ext)
plot(DEM_base, col = gray.colors(20, start = 0, end = 1))

#B. Create slope, aspect, and hillshade for the study area:
slope = terrain(DEM_base, opt="slope", unit='radians')
aspect = terrain(DEM_base, opt="aspect", unit='radians')
hillshade = hillShade(slope, aspect, angle=45, direction=315)
#plot hillshade to verify
plot(hillshade, col = gray.colors(20, start = 0, end = 1))
#add smoothing to hillshade
hillshade2 = aggregate(hillshade , fact = 5 , method = "bilinear" )
hillshade2 = focal(hillshade2, w=matrix(1/9, nc=3, nr=3), mean)
plot(hillshade2, col = gray.colors(20, start = 0, end = 1))
#save as raster
writeRaster(hillshade2, "maps/rasters/hillshade_carajas.tif", driver='GTiff')

#C. Add a smoothing effect to slope and save it
slope2 = aggregate(slope , fact = 5 , method = "bilinear" )
slope2 = focal(slope2, w=matrix(1/9, nc=3, nr=3), mean)
#verify
plot(slope2, col = gray.colors(20, start = 0, end = 1))
#save as raster
writeRaster(slope2, "maps/rasters/slope_carajas.tif", driver='GTiff')

#D. If you do not have a landcover raster, download it from NaturalEarth site
landcover = ne_download(scale = 10, type = 'NE1_LR_LC', category="raster")
#save as raster
writeRaster(landcover, "maps/rasters/LC_NaturalEarth.tif", bandorder='BIL', driver='GTiff', overwrite=TRUE)


###1.5. LOAD RASTERS AND SHAPEFILES FOR MAPS AND CORRECT THE CRS
#A. Load, crop and verify the CRS for all rasters:
#landcover1
#landcover = stack("maps/rasters/LC_NaturalEarth.tif") #multi band raster as satelite images
#landcover = crop(landcover, ext)
#proj4string(landcover)

#landcover2
landcover2 = stack("maps/rasters/NoSea_CL.tif") #multi band raster as satelite images
landcover2 = crop(landcover2, ext)
proj4string(landcover2)

#landcover3
landcover3 = stack("maps/rasters/LC_Carajas_2018.tif")
landcover3 = crop(landcover3, ext)
proj4string(landcover3)

#hillshade
hillshade = stack("maps/rasters/hillshade_carajas.tif")
hillshade = crop(hillshade, ext)
proj4string(hillshade)

#slope
#slope = stack("maps/rasters/slope_carajas.tif")
#slope = crop(slope, ext)
#proj4string(slope)

#if you need to fix the CRS base on a model, also a raster:
raster_to_fix = spTransform(raster_to_fix, proj4string(raster_model))

#B. Load, crop and verify the CRS for all shapefiles:
#cangas
cangas = shapefile("maps/shapefiles/Cangas.shp")
#if your shapefile is smaller than study area you do not need to crop. If you do it, the object will be NULL
#cangas = crop(cangas, ext) 
crs(cangas)

#if you need to fix the CRS base on a model, also a raster:
cangas = spTransform(cangas, proj4string(landcover2))
crs(cangas)


###1.6. COMBINE THE ELEMENTS TO TEST CONSTRAST, BRIGHTNESS, AND OTHER PARAMETERS FOR THE MAPS
#A. test the plot for a new map. For Carajas area and maps. These are a sample:
plot.new()

plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

dev.off()


#------------------------------------------------------------------------------
#                           2. Mapping Genetic Clusters 
#------------------------------------------------------------------------------
###2.1. MAPS USING sNMF POPULATIONS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_sNMF_pops_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axis
plot.window(xlim= axis.x, ylim= axis.y)

#B. Plot maps for the base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#C. Plot long and lat axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos = -6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos = -50.7)

#D. Plot the points:
#set population names and colors. SAME THAN STEP 2!
col_snmf= c('blue', 'red', 'darkslategrey', "yellow")
legend_snmf = c("POP1", "POP2", "POP3", "POP4")

#plot points. One line for each populations. Add or remove lines according to your necessity.
points(df$longitude[pop_snmf_1], df$latitude[pop_snmf_1], col = 'black', bg = col_snmf[1], cex = 2.5, pch=21) #POP1
points(df$longitude[pop_snmf_2], df$latitude[pop_snmf_2], col = 'black', bg = col_snmf[2], cex = 2.5, pch=21) #POP2
points(df$longitude[pop_snmf_3], df$latitude[pop_snmf_3], col = 'black', bg = col_snmf[3], cex = 2.5, pch=21) #POP3
points(df$longitude[pop_snmf_4], df$latitude[pop_snmf_4], col = 'black', bg = col_snmf[4], cex = 2.5, pch=21) #POP4


#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.08), (lat_names_serra-0.025), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=0.9, r=0.1)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.65, r=0.08, pos=4)
}

#G. Plot the legend:
legend(x=-49.85, y=-5.83, legend = legend_snmf, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_snmf, pch= 21, pt.cex=1.5)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.47), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)

#I. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()


###2.2. MAPS USING DAPC POPULATIONS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_DAPC_pops_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axis
plot.window(xlim= axis.x, ylim= axis.y)

#B. Plot maps for the base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#C. Plot lat and long axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos=-6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos=-50.7)

#D. Plot the points:
#set population names and colors. SAME THAN STEP 2!
col_dapc= c('red', 'blue', 'yellow', "darkslategrey")
legend_dapc = c("POP1", "POP2", "POP3", "POP4")
#plot points. One line for each populations. Add or remove lines according to your necessity.
points(df$longitude[pop_DAPC_1], df$latitude[pop_DAPC_1], col = 'black', bg = col_dapc[1], cex = 2.5, pch=21) #POP1
points(df$longitude[pop_DAPC_2], df$latitude[pop_DAPC_2], col = 'black', bg = col_dapc[2], cex = 2.5, pch=21) #POP2
points(df$longitude[pop_DAPC_3], df$latitude[pop_DAPC_3], col = 'black', bg = col_dapc[3], cex = 2.5, pch=21) #POP3
points(df$longitude[pop_DAPC_4], df$latitude[pop_DAPC_4], col = 'black', bg = col_dapc[4], cex = 2.5, pch=21) #POP4

#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.08), (lat_names_serra-0.025), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=0.9, r=0.1)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.65, r=0.08, pos=4)
}

#G. Plot the legend:
legend(x=-49.85, y=-5.85, legend = legend_dapc, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_dapc, pch= 21, pt.cex=1.5)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.47), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)

#I. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()


###2.3. MAPS USING TESS POPULATIONS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_TESS_pops_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axis
plot.window(xlim= axis.x, ylim= axis.y)

#B. Plot maps for the base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#C. Plot lat and long axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos=-6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos=-50.7)

#D. Plot the points:
#set population names and colors. SAME THAN STEP 2!
col_tess= c('red', 'pink', "yellow", 'blue', 'darkslategrey')
legend_tess = c("POP1", "POP2", "POP3", "POP4", "POP5")

#plot points. One line for each populations. Add or remove lines according to your necessity.
points(df$longitude[pop_TESS_1], df$latitude[pop_TESS_1], col = 'black', bg = col_tess[1], cex = 2.5, pch=21) #POP1
points(df$longitude[pop_TESS_2], df$latitude[pop_TESS_2], col = 'black', bg = col_tess[2], cex = 2.5, pch=21) #POP2
points(df$longitude[pop_TESS_3], df$latitude[pop_TESS_3], col = 'black', bg = col_tess[3], cex = 2.5, pch=21) #POP3
points(df$longitude[pop_TESS_4], df$latitude[pop_TESS_4], col = 'black', bg = col_tess[4], cex = 2.5, pch=21) #POP4
points(df$longitude[pop_TESS_5], df$latitude[pop_TESS_5], col = 'black', bg = col_tess[5], cex = 2.5, pch=21) #POP5


#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.08), (lat_names_serra-0.025), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=0.9, r=0.1)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.65, r=0.08, pos=4)
}

#G. Plot the legend:
legend(x=-49.85, y=-5.85, legend = legend_tess, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_tess, pch= 21, pt.cex=1.5)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.47), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)

#I. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()


#------------------------------------------------------------------------------
#                     3. Mapping Ancestry Coefficients as Pie Chart
#------------------------------------------------------------------------------
###3.1. MAPS USING sNMF COEFFICIENTS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_sNMF_coeficients_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axes
plot.window(xlim= axis.x, ylim= axis.y)

#B. Plot maps for the base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#C. Plot long and lat axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos=-6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos=-50.7)

#D. Plot the pies charts:
#set population names and colors. SAME THAN STEP 2!
legend_snmf 
col_snmf
#plot the pies. Z must be a matrix! Radius is the size of pies!
names(df) #verify number of cols with ancestry information
draw.pie(z=as.matrix(df[8:11]), x=df$longitude, y=df$latitude, radius = 0.02, col=col_snmf, labels="")

#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.08), (lat_names_serra-0.025), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=0.9, r=0.1)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.65, r=0.08, pos=4)
}

#G. Plot the legend:
legend(x=-49.85, y=-5.83, legend = legend_snmf, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_snmf, pch= 21, pt.cex=1.5)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.47), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)

#I. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()


###3.2. MAPS USING TESS COEFFICIENTS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_TESS_coeficients_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axes
plot.window(xlim= axis.x, ylim= axis.y)

#B. Plot maps for the base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#C. Plot long and lat axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos=-6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos=-50.7)

#D. Plot the pies charts:
#set population names and colors. SAME THAN STEP 2!
legend_tess
col_tess
#plot the pies. Z must be a matrix! Radius is the size of pies!
names(df) #verify number of cols with ancestry information
draw.pie(z=as.matrix(df[18:22]), x=df$longitude, y=df$latitude, radius = 0.02, col=col_tess, labels="")

#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.08), (lat_names_serra-0.025), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=0.9, r=0.1)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.65, r=0.08, pos=4)
}
  
#G. Plot the legend:
legend(x=-49.85, y=-5.85, legend = legend_tess, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_tess, pch= 21, pt.cex=1.5)
  
#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.47), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)
  
#I. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)
  
#J. Turn off pdf:
dev.off()


#------------------------------------------------------------------------------
#                  4. Mapping Ancestry Coefficients as Interpolation
#------------------------------------------------------------------------------
###4.1. MAPS USING TESS COEFFICIENTS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_sNMF_coeficients_Interpolation_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()

#B. Create objects to plot maps
#Define the coordinates as matrix
names(df) #verify number of cols with ancestry information
coordinates = data.matrix(df[,4:5])
rownames(coordinates) = NULL
class(coordinates)

#Define colors and pop names:
legend_snmf
col_snmf
my.palette = CreatePalette(col_snmf, 9)

#C. Plot maps for the base
#Plot the Interpolate map. Select colunms with ancestry coefficient:
plot(as.qmatrix(df[,8:11]), coordinates, method = "map.max", interpol = FieldsKrigModel(10),
     main = NULL, xlab = NULL, ylab = NULL, ann=FALSE,
     resolution = c(700,700), cex = 1.5,
     col.palette = my.palette, add=T)

plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#D. Plot the points:
#set population names and colors. SAME THAN STEP 2!
points(df$longitude[pop_snmf_1], df$latitude[pop_snmf_1], col = 'black', bg = col_snmf[1], cex = 1.5, pch=21) #POP1
points(df$longitude[pop_snmf_2], df$latitude[pop_snmf_2], col = 'black', bg = col_snmf[2], cex = 1.5, pch=21) #POP2
points(df$longitude[pop_snmf_3], df$latitude[pop_snmf_3], col = 'black', bg = col_snmf[3], cex = 1.5, pch=21) #POP3
points(df$longitude[pop_snmf_4], df$latitude[pop_snmf_4], col = 'black', bg = col_snmf[4], cex = 1.5, pch=21) #POP4

#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.00), (lat_names_serra-0.0), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=1, r=0.07, pos=2)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.7, r=0.06, pos=4)
}

#G. Plot the legend:
legend(x=-50.1, y=-5.95, legend = legend_snmf, cex = 0.6, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_snmf, pch= 21, pt.cex=1.3)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.39), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.6)

#I. Plot a north arrow:
north.arrow(xb=-50, yb=-6.42, len=0.008, cex.lab=0.4, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()


###4.2. MAPS USING TESS COEFFICIENTS
#A. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_Maps/Map_TESS_coeficients_Interpolation_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()

#B. Create objects to plot maps
#Define the coordinates as matrix
names(df) #verify number of cols with ancestry information
coordinates = data.matrix(df[,4:5])
rownames(coordinates) = NULL
class(coordinates)

#Define colors and pop names:
legend_tess
col_tess
my.palette = CreatePalette(col_tess, 9)

#C. Plot maps for the base
#Plot the Interpolate map. Select colunms with ancestry coefficient:
plot(as.qmatrix(df[,18:22]), coordinates, method = "map.max", interpol = FieldsKrigModel(10),
     main = NULL, xlab = NULL, ylab = NULL, ann=FALSE,
     resolution = c(700,700), cex = 1.5,
     col.palette = my.palette, add=T)

plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =40, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")


#D. Plot the points:
#set population names and colors. SAME THAN STEP 2!
points(df$longitude[pop_TESS_1], df$latitude[pop_TESS_1], col = 'black', bg = col_tess[1], cex = 1.5, pch=21) #POP1
points(df$longitude[pop_TESS_2], df$latitude[pop_TESS_2], col = 'black', bg = col_tess[2], cex = 1.5, pch=21) #POP2
points(df$longitude[pop_TESS_3], df$latitude[pop_TESS_3], col = 'black', bg = col_tess[3], cex = 1.5, pch=21) #POP3
points(df$longitude[pop_TESS_4], df$latitude[pop_TESS_4], col = 'black', bg = col_tess[4], cex = 1.5, pch=21) #POP4
points(df$longitude[pop_TESS_5], df$latitude[pop_TESS_5], col = 'black', bg = col_tess[5], cex = 1.5, pch=21) #POP5

#E. Plot the names of Serras in Carajás or any area names according to shapefile information
posi_names_serra = which(!is.na(cangas$Serra2))
lon_names_serra = coordinates(cangas)[,1][posi_names_serra]
lat_names_serra = coordinates(cangas)[,2][posi_names_serra]
#plot names:
shadowtext((lon_names_serra-0.00), (lat_names_serra-0.0), cangas$Serra2[posi_names_serra], col="black", bg="white", cex=1, r=0.07, pos=2)

#F. Plot the names of localities according to vcf information in metafile
localities = df$local_ID
for (i in 1:length(localities)){
  #pick the position of the first sample from that localities
  posi_local = which(df$local_ID %in% localities[i] == TRUE)[1]
  shadowtext((df$longitude[posi_local]+0.0), (df$latitude[posi_local]+0.0), localities[i], col="black", bg="white", cex=0.7, r=0.06, pos=4)
}

#G. Plot the legend:
legend(x=-50.1, y=-5.95, legend = legend_tess, cex = 0.6, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_tess, pch= 21, pt.cex=1.3)

#H. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.39), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.6)

#I. Plot a north arrow:
north.arrow(xb=-50, yb=-6.42, len=0.008, cex.lab=0.4, lab="N", col='black', font.sub=1)

#J. Turn off pdf:
dev.off()

#END
