# Spatial Network Models ----
# 10 KM
# Joao Braga et al.
rm(list = ls())

library(raster) 
library(Matrix)
library(rgeos)
library(rgdal)
library(foreign)
library(NetIndices)
library(cheddar)
library(ade4)
library(ggplot2)
library(reshape)
library(vegan)

# * WD ----
loc <- "/Users/braga/Documents/FACULDADE/Universite Joseph Fourier/PhD"
setwd(paste0(loc,"/FoodWebs/data/BioticData/"))


source(paste0(loc,"/FoodWebs/R/FunctionsHDD.R"))

# *Loading species distributions @ 10K ----
# Format: rows are pixel, columns are species, first column is pixel ID (PAGENAME)
load("feeding/Spp_traits_habs/Species_presence/MASTER.bin10000_allhab_tresh0.Rdata")
# Making sure pixel ID are characterers
master$PAGENAME <- as.character(master$PAGENAME)
# Give pixel ID as rowname
rownames(master) <- master$PAGENAME

# * Loading 10k resolution grid with pixel IDs
template <- read.dbf(file = "./10k/grid/grid_10Km.dbf")

# *Loading species properties (NETWORKS10KALLHAB.R output) ----
# Format rows are pixel, 
load('10k/froggy_outputs/networkpropEUROPE_plus_basalDecomp.Rdata')
head(finalDATA)

# Network properties columns into numerical ----
interORD <- finalDATA
interORD <- as.data.frame(interORD)
interORD$Pagename <- as.character(interORD$Pagename) 

for(i in 2:ncol(interORD)){
  interORD[,i] <- as.numeric(as.character(interORD[,i]))
}

# Ordering rows by the original raster/dbf ----
interORD <- interORD[rownames(master),]

# Remove big lakes from dataset ----
lakes10k <- raster(x ="~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/data/envData/wise_large_lakes/10k/lakes10k.img")

interORD[lakes10k[which(!is.na(values(lakes10k)))],-1] <- NA

# Loading Environmental variables ----
# Worldclim
tempY10 <- read.table(file= paste0(loc,"/FoodWebs/data/envData/bio/10k/temp10k.txt"),header = TRUE,sep = ';')
tempSEASON10 <- read.table(file=paste0(loc,"/FoodWebs/data/envData/bio/10k/tempSeas10k.txt"),header = TRUE,sep = ';')
precY10 <- read.table(file=paste0(loc,"/FoodWebs/data/envData/bio/10k/PrecipY10k.txt"),header = TRUE,sep = ';')
amptem <-read.table(file=paste0(loc,"/FoodWebs/data/envData/bio/10k/tempAmplitude10k.txt"),header = TRUE,sep = ';')
precipCV <- read.table(file=paste0(loc,"/FoodWebs/data/envData/bio/10k/BIO15.txt"),header = TRUE,sep = ';')

# Habitat variables
habs <- read.dbf(file = "10k/tabulateGLC_10Km.dbf") # in area occupied 
head(habs)

Nhabs <- specnumber(habs[,-1],MARGIN = 1)
shannon_index <- diversity(x = habs[,-1],index = "shannon",MARGIN = 1)
evenness <- shannon_index/log(specnumber(habs[,-1],MARGIN = 1))
rownames(habs) <- as.character(habs$PAGENAME) # Assigning rownames
habs$Nhabs <- Nhabs
habs$shannon_index <- shannon_index
habs$simpson_index <- simpson_index
habs$evenness <- evenness

# Loading NPP
NPP <- read.table(file = paste0(loc,"/FoodWebs/data/envData/NPP/NPP.txt"),header = TRUE)
NPP <- NPP[-119179,] # last line is NAs
str(NPP)
NPP$PageName <- as.character(NPP$PageName)
rownames(NPP) <- as.character(NPP$PageName)

# Loading Human footprint
HF <- read.csv(file=  paste0(loc,"/FoodWebs/data/envData/Human_Footprint/hfp_global_geo_grid/10k/HF_10k.txt"),header = TRUE,sep = ";")
head(HF)
HF <- HF[,-1]
HF$PAGENAME <- as.character(HF$PAGENAME)
rownames(HF) <- as.character(HF$PAGENAME)

# Extracting coordinates (X and Y)
griid <- shapefile(x = "10k/grid/grid_10Km.shp")
coords <- griid@polygons[[1]]@Polygons[[1]]@coords
coords <- data.frame(PAGENAME= griid@data$PageName,X = coordinates(griid)[,1], Y = coordinates(griid)[,2],row.names = griid@data$PageName)

# Predictors data frame ----
EnvData <- data.frame(PAGENAME= as.character(tempY10$PAGENAME), 
                      TempMean = tempY10$MEAN, 
                      TempSeason = tempSEASON10$MEAN, 
                      Precip = precY10$MEAN, 
                      precipCV = precipCV$MEAN,
                      Amplitude = amptem$MEAN,
                      shannon_index =  habs[as.character(tempY10$PAGENAME),'shannon_index'],
                      NHabs = habs[as.character(tempY10$PAGENAME),'Nhabs'],
                      Evenness = habs[as.character(tempY10$PAGENAME),'evenness'],
                      X = coords[as.character(tempY10$PAGENAME),2],
                      Y = coords[as.character(tempY10$PAGENAME),3],
                      NPP = NPP[as.character(tempY10$PAGENAME),5],
                      HF = HF[as.character(tempY10$PAGENAME),5]
                      )

# Making sure that pixel ID is a character (not factor)
EnvData$PAGENAME <-as.character(EnvData$PAGENAME)
rownames(EnvData) <- EnvData$PAGENAME # assigning row names 

# CLEAN PREDICTOR DATA OF NAS
EnvDataClean <- EnvData[complete.cases(EnvData),]
dim(EnvDataClean)

# Predictors basic stats ----
summary(EnvDataClean)

#
# Local food web properties extent consistent with EnvData extent ----
networkprops <- interORD[rownames(EnvDataClean),] 
networkprops <- networkprops[complete.cases(networkprops),] # Only pixels with data
EnvDataClean <- EnvDataClean[rownames(networkprops),]

# Calculate link density 
networkprops$LinkDensity <- networkprops$L/networkprops$S

# Saving data as Rdata ----
save(networkprops, file = '../../BioticData/10k/Network_properties_with_basal_decomposition.Rdata')
save(EnvDataClean, file = '../../BioticData/10k/Environmental_table.Rdata')

#
# Follow up ----
# 1 Null models (cluster) script: NETWORKS10KRANDOM_TAXA.R
# 2 Infer statistical significance from null (with FDR): NullModel.R
# 3 Run PCA with significant variables (Model.R)
# 4 Run models                         (Model.R)
# 5 Run variable importance            (Variable Importance.R)

# 1 Null models Preparation ----
# Creating of data frame with combinations of number of species and proportion of local taxa
# Load local composition
load("../network10kcompositionALLstages.Rdata")
composition_10k <- finalDATA[rownames(networkprops),-1] ; rm(finalDATA) # First column is pixel ID (Pagename)

load('../Network_properties_with_basal_decomposition.Rdata') # food web metrics over study area (pixel by metric)

# Quantifying proportion of birds, mammals, reptiles and amphibian using sp codes
# e.g. A21 = amphibian species 21
taxa_type <- substr(colnames(composition_10k), 1, 1)
taxa_S_pix_numeric <- apply(X = sppXpix,MARGIN = 2, 
                            FUN = function(x) {
                              setNames(object = as.vector(table(taxa_type[as.logical(x)])),
                                       names(table(taxa_type[as.logical(x)])))
                            }
)


taxa_type_df <- data.frame(matrix(NA,ncol = 4,nrow = length(taxa_S_pix_numeric)),row.names = names(taxa_S_pix_numeric)) ; colnames(taxa_type_df) <- c("A","B","M","R")

for(x in names(taxa_S_pix_numeric)){
  taxa_type_df[x,names(taxa_S_pix_numeric[[x]])] <-  taxa_S_pix_numeric[[x]]
} 
taxa_type_df[is.na(taxa_type_df)] <- 0 # Empty pixels get zero for taxa

save(taxa_type_df,file = "10k/random networks/propTaxa.Rdata")

# end ----