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

results.path <- "10k/Random_webs_TAXA/"
dir.create(results.path, recursive = TRUE)

source("Functions.R")
#################################################
# Loading the metaweb
# Diet
load("BARMdiet_bin.RData")
#################################################
# Loading habitat per species
BARM.HAB <- read.csv("BARM_allhabs.csv",header = TRUE,row.names = 'ID',sep = ';')

#################################################
# Loading habitat per species (from Data_preparation.R)
load("propTaxa.Rdata")

# Species richness
S_richness <- rowSums(taxa_type_df)      # species richness across Europe
valuesS <- unique(S_richness)   # Interval of species richness across Europe

#################################################
# Subset habitats by spp and habitats present in pixel
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
dietcat <- head(row.names(BARMdiet.binary),12) # first 12 columns and rows are diet categories in metaweb

# Function to sample a number of species based on an empirical distirbution from the metaweb
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

sampling_spp <- function(y,                # number of species
                         BARMdiet.binary,  # metaweb
                         sampled_line = NULL){
  if(is.null(sampled_line)) {
    lines.to.sample <- which(S_richness == y)
    sampled_line <- taxa_type_df[resamp(lines.to.sample,size = 1),]
    }
   
spp_pools <- colnames(BARMdiet.binary)[!(colnames(BARMdiet.binary) %in% dietcat)]  # Species within metaweb that are not dietcats

# Sampling species
# Note the I multiple the number of taxa by a constant. This makes me sample always a bit more species than necessary. This is needed, because many species are disconnected after building the null food web. 
Asp <- if(sampled_line$A != 0)
  sample(x = grep(pattern = "A",spp_pools,value = TRUE), round(sampled_line$A*1.08),replace = F) else   
    NULL

Bsp <- if(sampled_line$B != 0 )
  sample(x = grep(pattern = "B",spp_pools,value = TRUE),round(sampled_line$B*1.02),replace = F) else
    NULL

Msp <-  if(sampled_line$M != 0 )
  sample(x = grep(pattern = "M",spp_pools,value = TRUE),round(sampled_line$M*1.08),replace = F) else 
    NULL

Rsp <- if(sampled_line$R != 0 )
  sample(x = grep(pattern = "R",spp_pools,value = TRUE),round(sampled_line$R*1.05),replace = F) else 
    NULL

# Output as list with two vectors
# species refers to the species sampled
# sampled_line refers to the proportion of taxa present in a particular local food web
spp_web <- list() ; spp_web$species <- c(Asp,Bsp,Msp,Rsp)
spp_web$Line_taxa <- sampled_line        # 
return(spp_web)
}

# RUN the generation of random food webs
replicates <- rep(valuesS,300000) # There are 300 unique species values times 1000 replicates.
  

fun.to.par <- function(y) {
    results.dir <- paste0(results.path,y)
    dir.create(results.dir, recursive = TRUE)  
  
# First run builds the random based on y (number of species)
    spp_web <- sampling_spp(y,BARMdiet.binary )    # Samples that y number of species
    spp_web_diet <- c(dietcat, spp_web$species)    # Adds diet categories to build food web
    
    hab <- pix.hab[1,] # needed for the code to run, although no selection by the habitat
    web <- sub_web(metaweb = BARMdiet.binary,SPPCODE = spp_web_diet,PIX.HAB = hab, SPP.HAB = BARM.HAB,HELP = FALSE)             # Select those species from the network to a new network + Diet categories
    
    spp <- colnames(web)[!(colnames(web) %in% dietcat)]  # Which species remained
    
    if(!is.matrix(web)) web <- as.matrix(web)
    
    same_tg <- length(spp_web$Line_taxa[spp_web$Line_taxa!= 0]) == length(table(substr(spp, 1, 1))) # calculating the proportion of taxa in the web
    dif_tg_sr <- if(same_tg) {sum(abs(spp_web$Line_taxa[spp_web$Line_taxa!= 0] - table(substr(spp, 1, 1)))) != 0 } else TRUE  # logic statement to match proportion of taxa between local and random
    
    # Rerun the sampling until logic statement is true (i.e. local and random proportion of taxa are equal)
    while(dif_tg_sr){   
      spp_web <- sampling_spp(y,BARMdiet.binary,sampled_line = spp_web$Line_taxa)
      spp_web_diet <- c(dietcat, spp_web$species)
      # If spp are found in that pixel, then add the diet categories to the analysis
      hab <- pix.hab[1,]
      web <- sub_web(metaweb = BARMdiet.binary,SPPCODE = spp_web_diet,PIX.HAB = hab,SPP.HAB = BARM.HAB,HELP = FALSE)             # Select those species from the network to a new network + Diet categories
      spp <- colnames(web)[!(colnames(web) %in% dietcat)]  
      if(!is.matrix(web)) web <- as.matrix(web)
    
      same_tg <- length(spp_web$Line_taxa[spp_web$Line_taxa != 0]) == length(table(substr(spp, 1, 1)))
      dif_tg_sr <- if(same_tg) {sum(abs(spp_web$Line_taxa[spp_web$Line_taxa!= 0] - table(substr(spp, 1, 1))))!=0 } else TRUE
    }
    
    # Calculate food web properties for random food web
    res <- fun.net.prop(t(web),dietcat = dietcat)

    # Save output based on y number of species and a random unique ID number
  write.table(res, file = paste0(results.dir,"/Rand",y,"_",round(runif(n = 1,min = 0,max = 500000)),".txt"))
}

NCOREs <- 48   # Cores used
date()
mclapply(X =  replicates,FUN = fun.to.par,mc.cores = NCOREs)
date()


################################################
# Assemble the table with results
compile_results <- function(files.loc.vector = fls,y){
  flsnames <- gsub(pattern = ".txt$",replacement = "",x = files.loc.vector)
  flsnames <- gsub(pattern = paste0(results.path,y,"//"),replacement = "",x = flsnames)

  print("Loading files into memory")

  PropsVector <- vector("list", length(files.loc.vector))
  for(i in 1:length(files.loc.vector)){
    PropsVector[[i]] <- if(is.vector(loadTxt(filePath = files.loc.vector[i]))) {eval(expr = parse(text=paste0("loadTxt(filePath = files.loc.vector[i])")))} else
       structure(1,names="NA")
  }
  print("Files are loaded. Bind output into matrix format.")
  PropDF <- do.call(foo, PropsVector)
  rownames(PropDF) <- flsnames
  return(PropDF)
}

fun.to.par <- function(y){
  simsY <-  list.files(path = paste0(results.path,y,"/"), pattern = ".txt$",full.names = TRUE)
  SimTab <- compile_results(files.loc.vector =simsY,y)
  save(SimTab,file = paste0(results.path,y,".Rdata"))
  }

date()
mclapply(X =  valuesS,FUN = fun.to.par,mc.cores = NCOREs)
date()

q('no')