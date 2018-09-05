#################################################
# Master tables for presence/absence 
# 10KM
# Joao Braga 22 oct 2015
#################################################
# Define habitat threshold (means higher than #, e.g. k = 0 than % higher than 0)
# This threshold defines the presence of a species
# If habitat proportion is equal or higher than k, species will be present.
k <- 0 # 10 or 20 

for(k in c(1:19,30,40,50)){
# Table name to save
master.file <- paste0("MASTER.bin10000_allhab_tresh",k,".Rdata")

# Data Working directory
setwd("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/data/")

# Requires
library(raster) 
library(Matrix)
library(rgeos)
library(rgdal)
library(foreign)

# Run this
# Modify only to specify if all or opt habitat was considered.
#####
source("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/R/froggyscripts/Functions.R")
#####
# Loading the mask and the dbf file
mask10kID <-read.dbf(file = "BioticData/10k/mask/reference_grid_10km.img.vat.dbf")
mask10k <- raster(x = "BioticData/10k/mask/reference_grid_10km.img")

# Loading the spatial distribution of species 
SPP.dist.files <- list.files(path = "BioticData/10k/COMPLETE/",full.names = T,pattern = ".RData$")
spp.dist <- fun.load.many.rdata(file.list = SPP.dist.files)

# Indexed species codes
sppCodes <- names(spp.dist)
sppCodes <- gsub(pattern = "^tab_m",replacement = "",x = sppCodes)
sppCodes <- gsub("(\\w)","\\U\\1",sppCodes,perl=TRUE)
names(spp.dist) <- sppCodes

# Creating a master matrix with CELL ID and presence/absence (given a certain % of primary habitat threshold) with all species in the study area
# Run me (less than 1 minute)
date()
master <- as.data.frame(matrix(data = NA,nrow = nrow(mask10kID),ncol = 1+length(spp.dist)))
names(master) <- c("PAGENAME",sppCodes)
master[,1] <- mask10kID$PageName
for(i in 1:length(spp.dist)){
  test <- fun.PRESENCE.ABSENCE(species.dist = spp.dist[[i]],threshold = k,opt.only = FALSE)
  cellID <- as.character(test$PAGENAME)
  master[master$PAGENAME %in% cellID,i+1] <- test[test$PAGENAME %in% cellID,"Present"] 
}
date()

save(master,file = paste0("BioticData/feeding/Spp_traits_habs/Species_presence/",master.file))
}
#####
# End