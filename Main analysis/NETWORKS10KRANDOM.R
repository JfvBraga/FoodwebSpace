#################################################
# Spatial network
# 10 KM
# Joao Braga 6 May 2016
#################################################
rm(list = ls())

# Requires
library(parallel)
library(NetIndices, lib.loc = "/home/bragaj/froggy_libs")
library(igraph, lib.loc = "/home/bragaj/froggy_libs")
library(cheddar, lib.loc = "/home/bragaj/froggy_libs")

library(foreign)
# Data Working directory
setwd("/scratch/bragaj/BioticData/")

results.dir <- "10k/Random webs/"
dir.create(results.dir, recursive = TRUE)

source("Functions.R")
source("FoodWebFunctions.R")
#################################################
# Loading the network
# Diet
load("BARMdiet_bin.RData")

#################################################
# Loading species dist at 10k for opt and secundary habitat of species
load('10k/MASTER.bin10000_allhab_tresh0.Rdata')

#################################################
# Loading habitat per species
BARM.HAB <- read.csv("BARM_allhabs.csv",header = TRUE,row.names = 'ID',sep = ';')

#################################################
# Loading habitat per pixel
#pix.hab <- read.dbf(file = '10k/tabulateGLC_10Km.dbf')
#pix.hab$PAGENAME <- as.character(pix.hab$PAGENAME)
#rownames(pix.hab) <- pix.hab$PAGENAME

#names(pix.hab)[-1] <- sub(pattern = "VALUE_",replacement = "X",x = names(pix.hab)[-1])

#pix.hab <- pix.hab[,-42] # remove the habitat category 230 (no habitat)

#pix.hab <- pix.hab[, c(1, which(names(pix.hab)%in% colnames(BARM.HAB)))]     # remove habitats that are not present in spp x habs table

#pix.hab[, -1] <- 1 # every habitat is present in a pixel



# subset habitats by spp and habitats present in pixel
dietcat <- head(row.names(BARMdiet.binary),12) # 13 in the case of eggscat or 12 withouth them

str(BARM.HAB)
BARM.HAB$SPPname <- as.character(BARM.HAB$SPPname)
BARM.HAB <- BARM.HAB[,-1]
BARM.HAB[is.na(BARM.HAB)] <- 0
sppXsppHAB <- as.matrix(BARM.HAB) %*%  t(as.matrix(BARM.HAB))
sppXsppHAB[sppXsppHAB >=1] <- 1


pix.hab <- BARM.HAB[1,]
pix.hab[,] <- 1
BARMdiet.binary <- BARMdiet.binary[c(dietcat,colnames(sppXsppHAB)),c(dietcat,colnames(sppXsppHAB))]

BARMdiet.binary[colnames(sppXsppHAB),colnames(sppXsppHAB)] <- BARMdiet.binary[colnames(sppXsppHAB),colnames(sppXsppHAB)] * sppXsppHAB[colnames(sppXsppHAB),colnames(sppXsppHAB)]


# Network properties
#################################################
newBARM.DIET <- as.matrix(BARMdiet.binary)

# In case of metaweb with dietcategory
# dietcat <- head(row.names(BARMdiet.binary),12) # 13 in the case of eggscat or 12 withouth them

valuess <- c(33,79,76,54,78,65,85,57,99,98,102,75,63,55,83,82,86,87,100,101,92,73,62,81,84,90,93,94,77,80,89,91,95,96,74,72,97,71,59,68,105,88
,67,61,69,70,104,106,56,58,60,103,107,110,111,108,119,64,109,123,114,115,124,126,128,118,116,66,122,131,117,125,127,130,113,129,133,112,120,135,136,134,121,137
,132,53,138,139,45,140,141,142,147,149,146,144,154,152,148,145,150,143,159,165,162,160,158,163,153,161,151,157,156,155,168,164,167,170,169,171,166,173,172,174,176,175
,177,178,180,179,181,182,184,183,185,186,187,188,197,194,190,191,195,189,193,192,196,198,200,202,207,199,203,201,206,210,225,217,218,214,213,215,212,211,205,216,209,204
,208,219,27,46,221,220,47,52,222,223,224,42,227,43,28,226,228,230,229,232,231,233,234,236,44,237,235,238,239,240,241,34,242,243,244,245,246,10,247,51,49,248
,23,50,251,250,249,253,252,254,12,16,48,255,256,257,262,259,258,37,25,39,22,261,260,21,264,263,266,265,267,41,268,269,13,36,35,32,38,270,272,279,274,273
,276,275,293,284,280,287,283,285,294,288,291,295,297,286,281,278,277,296,298,289,271,282,290,292,305,26)


# RUN
need.to.done <- rep(valuess,30000)
  
NCOREs <- 48

fun.to.par <- function(y){
    results.dir <- paste0("10k/Random webs/",y)
    dir.create(results.dir, recursive = TRUE)  
  
    dddd <- sample(colnames(BARMdiet.binary),y+y*0.05,replace = FALSE)
    
    while(length(intersect(x = dddd,dietcat))>= 1) {
      dddd <- sample(colnames(BARMdiet.binary),y+y*0.05,replace = FALSE)
      }
    
    dddd <- c(dietcat, dddd)                     # If spp are found in that pixel, then add the diet categories to the analysis
    
    hab <- pix.hab[1,]
 
    mmmm <- sub_web(metaweb = BARMdiet.binary,SPPCODE = dddd,PIX.HAB = hab,SPP.HAB = BARM.HAB,HELP = FALSE)             # Select those species from the network to a new network + Diet categories
    
    spp <- colnames(mmmm)[!(colnames(mmmm) %in% dietcat)]  
    # print(length(spp))
    while(length(spp) != y){
      
      dddd <- sample(colnames(BARMdiet.binary),y+y*0.05,replace = FALSE)
      
      while(length(intersect(x = dddd,dietcat))>= 1) {
        dddd <- sample(colnames(BARMdiet.binary),y+y*0.05,replace = FALSE)
      }
      
      dddd <- c(dietcat, dddd)                     # If spp are found in that pixel, then add the diet categories to the analysis
      
      hab <- pix.hab[1,]
      
      mmmm <- sub_web(metaweb = BARMdiet.binary,SPPCODE = dddd,PIX.HAB = hab,SPP.HAB = BARM.HAB,HELP = FALSE)             # Select those species from the network to a new network + Diet categories
      
      spp <- colnames(mmmm)[!(colnames(mmmm) %in% dietcat)]  
      # print(length(spp))
    }
    
    if(!is.matrix(mmmm)) mmmm <- as.matrix(mmmm)
    
    res <- fun.net.prop(t(mmmm),dietcat = dietcat)

  write.table(res, file = paste0(results.dir,"/Rand",y,"_",round(runif(n = 1,min = 0,max = 500000)),".txt"))
}

#date()
#mclapply(X =  need.to.done,FUN = fun.to.par,mc.cores = NCOREs)
#date()

################################################
# Assemble the table

valuess <- c(33,79,76,54,78,65,85,57,99,98,102,75,63,55,83,82,86,87,100,101,92,73,62,81,84,90,93,94,77,80,89,91,95,96,74,72,97,71,59,68,105,88
             ,67,61,69,70,104,106,56,58,60,103,107,110,111,108,119,64,109,123,114,115,124,126,128,118,116,66,122,131,117,125,127,130,113,129,133,112,120,135,136,134,121,137
             ,132,53,138,139,45,140,141,142,147,149,146,144,154,152,148,145,150,143,159,165,162,160,158,163,153,161,151,157,156,155,168,164,167,170,169,171,166,173,172,174,176,175
             ,177,178,180,179,181,182,184,183,185,186,187,188,197,194,190,191,195,189,193,192,196,198,200,202,207,199,203,201,206,210,225,217,218,214,213,215,212,211,205,216,209,204
             ,208,219,27,46,221,220,47,52,222,223,224,42,227,43,28,226,228,230,229,232,231,233,234,236,44,237,235,238,239,240,241,34,242,243,244,245,246,10,247,51,49,248
             ,23,50,251,250,249,253,252,254,12,16,48,255,256,257,262,259,258,37,25,39,22,261,260,21,264,263,266,265,267,41,268,269,13,36,35,32,38,270,272,279,274,273
             ,276,275,293,284,280,287,283,285,294,288,291,295,297,286,281,278,277,296,298,289,271,282,290,292,305,26)

#for(y in valuess){
#  fls <- list.files(path = paste0("10k/Random webs/",y,"/"), pattern = ".txt$",full.names = TRUE)
#  
#  if(length(fls) >= 1000) fls <- sample(x = fls,size = 1000,replace = FALSE)
#
# 
#  finalDATA <- matrix(NA,nrow =length(flsnames),ncol = 17)
#  row.names(finalDATA) <- flsnames
#  colnames(finalDATA) <- rownames(read.table(file = fls[1],header = TRUE))
#  #
#  for(i in 1:length(fls)) {
#    finalDATA[i,] <- as.character(read.table(file = fls[i],header = TRUE)[1:17,1])
#  }
#  
#  save(finalDATA,file = paste0("10k/Random webs/",y,".Rdata"))
#}


fun.to.par <- function(y){
  
  simsY <-  list.files(path = paste0("10k/Random webs/",y,"/"), pattern = ".txt$",full.names = TRUE)

  if(length(simsY) >= 1000) fls <- sample(x = simsY,size = 1000,replace = FALSE)

  SimTab <- matrix(NA,nrow =length(simsY),ncol = 16)

  colnames(SimTab) <- rownames(read.table(file = simsY[1],header = TRUE))[1:16]
  
  for(i in 1:length(simsY)) SimTab[i,] <- as.character(read.table(file =simsY[i],header = TRUE)[1:16,1])
  
  save(SimTab,file = paste0("10k/Random webs/",y,".Rdata"))

  }

date()
mclapply(X =  valuess,FUN = fun.to.par,mc.cores = NCOREs)
date()

q('no')
