#
# Variables importance
#
#####################

rm(list = ls())

# Data Working directory
setwd("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/data/BioticData/10k/")

# Requires
library(raster) 
library(Matrix)
library(mgcv)
library(rgeos)
library(rgdal)
library(foreign)
library(NetIndices)
library(cheddar)
library(ade4)
library(ggplot2)
library(reshape)

# Variable importance
models <- list.files(path = "models/Revision/",pattern = "GAMmod.Rdata",full.names = TRUE)
RFmodels <- list.files(path = "models/Revision/",pattern = "^Random",full.names = TRUE)
modelsNAMES <- list.files(path = "models/Revision/",pattern = "GAMmod.Rdata",full.names = FALSE)
modelsNAMES <- gsub(pattern = "VARS",replacement = "",x = modelsNAMES)
modelsNAMES <- gsub(pattern = "models/Revision/",replacement = "",x = modelsNAMES)

# Load data
load("../10k/Network_properties_with_basal_decomposition.Rdata") # local food web properties
load("../10k/Environmental_table.Rdata")                         # Local environmental data

# keep same number of pixels and same order
EnvDataClean <- EnvDataClean[rownames(networkprops),]            
EnvDataClean <- EnvDataClean[complete.cases(EnvDataClean),]
networkprops <- networkprops[rownames(EnvDataClean),]

# Variable importance for Model 1
# Model axis 1
load(models[1]) #load model
load(RFmodels[1]) 
varnames <- c("HDNVAR","TempMean","TempSeason","Precip","precipCV","Evenness","shannon_index","NPP","HF")
data.M <- EnvDataClean
data.M$HDNVAR <- RF$y-RF$predicted
fitted_values <- fitted(gam.model)
data.M <- data.M[varnames]

# Function to parallel
fun.to.par.1 <- function(vari){

      res <- matrix(data = NA,nrow = 1000,ncol = 1)
      for(iii in 1:1000){

        fitted_values
      
        rand <-  data.M
        rand[,vari] <- sample(rand[,vari], length(rand[,vari]),replace = FALSE)
        preds <- predict.gam(object = gam.model,newdata = rand)
        res[iii,1]<- cor(fitted_values,preds,method = "spearman")
        
      }
      save(res,file = paste0("models/Revision/model1",vari,".Rdata"))
    }

# To run for multiple variables
mclapply(X = varnames,FUN = fun.to.par.1,mc.cores = 3)


# Model axis 2
load(models[2])
load(RFmodels[2])
varnames <- c("HDNVAR","TempMean","TempSeason","Precip","precipCV","Evenness","shannon_index","NPP","HF")
data.M <- EnvDataClean
data.M$HDNVAR <- RF$y-RF$predicted
fitted_values <- fitted(gam.model)
data.M <- data.M[,varnames]

# Function to parallel
fun.to.par.2 <- function(vari){
  
  res <- matrix(data = NA,nrow = 1000,ncol = 1)
  for(iii in 1:1000){
    rand <- data.M[,varnames]
    rand[,vari] <- sample(rand[,vari], length(rand[,vari]),replace = FALSE)
    preds <- predict.gam(object = gam.model,newdata = rand)
    res[iii,1]<- cor(fitted_values,preds,method = "spearman")
    
  }
  save(res,file = paste0("models/Revision/model2",vari,".Rdata"))
}

# To run for multiple variables
mclapply(X = varnames,FUN = fun.to.par.2,mc.cores = 3)


# Grab outputs & results ----
# Model Axis 1
var_import_outputs <- list.files("../10k/models/Revision/",pattern = "model1",full.names = TRUE)
varimport1 <- unlist(lapply(X = var_import_outputs, function(x) {
  load(x)
  return(1-mean(res))
  }))

names(varimport1) <- gsub(pattern = "../10k/models/Revision//model1",replacement = "" ,x = var_import_outputs) # assign the names for the var
varimport1
# Evenness.Rdata        HDNVAR.Rdata            HF.Rdata           NPP.Rdata        Precip.Rdata      precipCV.Rdata 
# 0.024971717         0.114091072         0.023601002         0.113166199         0.003003513         0.017350762 
# shannon_index.Rdata      TempMean.Rdata    TempSeason.Rdata 
# 0.078527904         0.408387250         0.457997179 

# Model Axis 2
var_import_outputs <- list.files("../10k/models/Revision/",pattern = "model2",full.names = TRUE)
varimport2 <-unlist(lapply(X = var_import_outputs, function(x) {
  load(x)
  return(1-mean(res))
}))

names(varimport2) <- gsub(pattern = "../10k/models/Revision//model2",replacement = "" ,x = var_import_outputs) # assign the names for the var
varimport2
# Evenness.Rdata        HDNVAR.Rdata            HF.Rdata           NPP.Rdata        Precip.Rdata      precipCV.Rdata 
# 0.0009274396        0.3668100368        0.2541902657        0.1110498345        0.0797488220        0.0258969112 
# shannon_index.Rdata      TempMean.Rdata    TempSeason.Rdata 
# 0.0751228458        0.3447961294        0.1096373395 

#End