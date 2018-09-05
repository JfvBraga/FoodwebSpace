#
# Spatial Network Models ----
# 10 KM
# Joao Braga
#
rm(list = ls())

# Requires
library(raster)
library(vegan)
library(foreign)

# *WD ----
# in Mac
loc <- "/Users/braga/Documents/FACULDADE/Universite Joseph Fourier/PhD"

source(paste0(loc,"/FoodWebs/R/FunctionsHDD.R"))
setwd(paste0(loc,"/FoodWebs/data/BioticData/"))

#
# *Loading species distributions @ 10K ----
load("feeding/Spp_traits_habs/Species_presence/MASTER.bin10000_allhab_tresh20.Rdata")
dim(master)
master$PAGENAME <- as.character(master$PAGENAME)
rownames(master) <- master$PAGENAME

# *Loading species properties (NETWORKS10KALLHAB.R output) ----
# load('10k/networkpropEUROPE.Rdata')
load('10k/froggy_outputs/networkpropEUROPE_20HAB.Rdata')   # Same as Europe (as in the paper, but with basal decomposition into two categories)
head(finalDATA)

template <- read.dbf(file = "10k/grid/grid_10Km.dbf")

# Network properties columns into numerical ----
interORD <- finalDATA
interORD <- as.data.frame(interORD)
interORD$Pagename <- as.character(interORD$Pagename)

for(i in 2:ncol(interORD)){
  interORD[,i] <- as.numeric(as.character(interORD[,i]))
}

# Ordering rows by the original raster/dbf ----
interORD <- interORD[rownames(master),]

# Loading EnvData ----
load(file = '../BioticData/10k/Environmental_table.Rdata')

networkprops <- interORD[rownames(EnvDataClean),]
networkprops <- networkprops[complete.cases(networkprops),] # Same pixel selection ----
EnvDataClean <- EnvDataClean[rownames(networkprops),]

save(networkprops, file = '../BioticData/10k/Network_properties_Thr20.Rdata')

#