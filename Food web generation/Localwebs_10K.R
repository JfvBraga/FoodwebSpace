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

source("Functions.R")
# Data Working directory
setwd("/scratch/bragaj/BioticData/")


results.dir <- "10k/NETWORKPROPSEUROPE_plus_basalDecomp/"
dir.create(results.dir, recursive = TRUE)
#################################################
# Cores for analysis
NCOREs <- 48

# Loading the metaweb (format: predators as rows, preys as columns)
#  with diet categories (12 rows and 12 columns added, one per diet category)
load("BARMdiet_bin.RData") 
# head(BARMdiet.binary)

#################################################
# Loading species dist at 10k for opt and secundary habitat of species
load('10k/MASTER.bin10000_allhab_tresh0.Rdata') # From MasterGen.R

# In case of different habitat habitat thresholds
#load('10k/MASTER.bin10000_allhab_tresh10.Rdata') 
#load('10k/MASTER.bin10000_allhab_tresh20.Rdata') 

#################################################
# Loading species habitat preferences 
BARM.HAB <- read.csv("BARM_allhabs.csv",header = TRUE,row.names = 'ID',sep = ';')

#################################################
# Loading habitat per pixel
pix.hab <- read.dbf(file = '10k/tabulateGLC_10Km.dbf')
pix.hab$PAGENAME <- as.character(pix.hab$PAGENAME)  # Pixel unique ID
rownames(pix.hab) <- pix.hab$PAGENAME               # Pixel unique ID as rownames

names(pix.hab)[-1] <- sub(pattern = "VALUE_",replacement = "X",x = names(pix.hab)[-1]) # Matching the names with BARM.HAB
pix.hab <- pix.hab[, c(1, which(names(pix.hab)%in% colnames(BARM.HAB)))]               # Removes habitats that are not present in species  habitat preferences table (BARM.HAB)

# Network properties ----
#################################################
master[is.na(master)] <- 0 # NAs become absences
# Species vector with their names
species <- colnames(master)[-1]  # First column are pixel names
master[,1] <- as.character(master[,1])  # Making sure the pixel IDs are characteres and not factors 

# * Store diet categories names as vectors ----
dietcat <- head(row.names(BARMdiet.binary),12) # 12 diet categories, first 12 rows (and columns) from BARMdiet.binary

# Function to parallelise ----  
#############################
fun.to.par <- function(y){

  # if a pixel in the master matrix is NA (outside europe or water), then NA                                                                       
  if(sum(master[y,2:ncol(master)]) == 0) { 
      # results for these pixel are NAs
      res  <- structure(c(master[y,1],rep(NA,18)), 
                        names=c("Pagename","S","L","C","clust.coef","char.path.length","propOmn","propB","propBNonH","propBH","propT","propI","GEN","GENSD","VUL","VULSD","MAXSIM","meanTL","maxTL")
                        )
    }     
  else{  # if not NA, then which species are present in that pixel
    
    composition <- species[which( master[y,2:length(master)]==1, arr.ind=FALSE)] # pixel y composition
    if(length(composition)!= 0) composition <- c(dietcat, composition)           # If spp are found in that pixel, then add the diet categories to the analysis                     
    hab <- pix.hab[y,] # Habitat composition for pixel y

    # Select those species + diet categories from the metaweb to form a local web for pixel y (adjacency matrix)
    local.web <- sub_web(metaweb = BARMdiet.binary,SPPCODE = composition,PIX.HAB = hab,SPP.HAB = BARM.HAB,HELP = FALSE)           
    # Making sure the output is a matrix
    if(!is.matrix(local.web)) local.web <- as.matrix(local.web)
    
    
    
    # local.web can be empty if no species remain (no resources; disconnected food web)
    if(sum(dim(local.web))!=0) res <- c(master[y,1],fun.net.prop(t(local.web),dietcat = dietcat))  # fun.net.prop calculates food web properties from a given local web. 
    else {
      # If local.web is empty, results for local food webs properties are NAs
      res  <- structure(c(master[y,1],rep(NA,18)), 
                        names=c("Pagename","S","L","C","clust.coef","char.path.length","propOmn","propB","propBNonH","propBH","propT","propI","GEN","GENSD","VUL","VULSD","MAXSIM","meanTL","maxTL")
                        )
    }
  }
  write.table(res, file = paste0(results.dir,master[y,1],".txt")) # Output is a text file with food web properies for a pixel y (based on y ID)
}
#############################

# To run
#############################
to.do <- 1:nrow(master)
date()
mclapply(X =  to.do,FUN = fun.to.par,mc.cores = NCOREs)
date()
#############################
# Assemble the table with results
# Compiling results into Rdata
fls <- list.files(path =results.dir, pattern = ".txt$",full.names = TRUE) # Listing output files

compile_results <- function(files.loc.vector = fls){
  # Clean and store the pixel names (IDs)
  flsnames <- gsub(pattern = ".txt$",replacement = "",x = flsnames)
  flsnames <- gsub(pattern = paste0(results.dir,"/"),replacement = "",x = flsnames)

  print("Loading files into memory")
  pb <- txtProgressBar(min = 0, max = length(files.loc.vector), style = 3)
  PropsVector <- vector("list", length(files.loc.vector))
  for(i in 1:length(files.loc.vector)){
    if(is.vector(loadTxt(filePath = files.loc.vector[i]))) PropsVector[[i]] <- eval(expr = parse(text=paste0("loadTxt(filePath = files.loc.vector[i])")))
    else PropsVector[[i]] <- structure(1,names="NA")
  }
  
  print("Files are loaded. Bind output into matrix format.")
  PropDF <- do.call(foo, PropsVector)
  rownames(PropDF) <- flsnames
  return(PropDF)
}

# Running the function to list and merge all output files
# Format: rows are pixels, columns are food web properties, first column is pixel ID
finalDATA <- compile_results(files.loc.vector = fls)

# Saving output table as a Rdata file
save(finalDATA,file = "networkpropEUROPE_plus_basalDecomp.Rdata")

q('no')