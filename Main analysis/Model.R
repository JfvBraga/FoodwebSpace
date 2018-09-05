# Analysis

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
library(vegan)

source("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/R/FunctionsHDD.R")
#################################################
# Load data
load("../10k/Network_properties_with_basal_decomposition.Rdata") # local food web properties
load("../10k/Environmental_table.Rdata")                         # Local environmental data

# keep same number of pixels and same order
EnvDataClean <- EnvDataClean[rownames(networkprops),]            
EnvDataClean <- EnvDataClean[complete.cases(EnvDataClean),]
networkprops <- networkprops[rownames(EnvDataClean),]

# Results from the null models and FDR approach to test significance of network properties (from NullModel.R)

#         S                L*              C*       clust.coe* char.path.lengt*          propOm*           propB* 
# 0.0000000        0.9992266        0.9992266        0.9964628        0.9983391        0.9577559        0.9565642 
# propBNonH           propBH            propT            propI*             GEN            GENSD              VUL 
# 0.5430998        0.4476957        0.0000000        0.9678225        0.8646212        0.8522472        0.0000000 
# VULSD           MAXSIM               meanTL*            maxTL      LinkDensity* 
# 0.0000000        0.7729192        0.9641204        0.0000000        0.9992266 

# Spearman correlation between food web properties
corrs <- round(cor(networkprops[,-1],method = "spearman"),digits = 2)
corrs[lower.tri(corrs)] <- NA
corrs

# keep same number of pixels and same order
EnvDataClean <- EnvDataClean[rownames(networkprops),]            
EnvDataClean <- EnvDataClean[complete.cases(EnvDataClean),]
networkprops <- networkprops[rownames(EnvDataClean),]

# Extracting coordinates (X and Y) needed to the models
griid <- shapefile(x = "grid/grid_10Km.shp")
template <- read.dbf(file = "./grid/grid_10Km.dbf")
coords <- cbind(EnvDataClean$X, EnvDataClean$Y)           # Extracting spatial X and Y coordinates
geopoints <- SpatialPoints(coords = coords,proj4string = crs(griid)) # converting coordinates into spatial object (to be used by the autocov model)


##################################################
# Running the PCA on the significantly different (from null model) variables
NTWPROPNewRes <- dudi.pca(df = networkprops[,c("S","char.path.length","LinkDensity","propOmn","C","propI","clust.coef","meanTL" )],scannf = F) # Characteristic path length with infinite values (removed)
summary(NTWPROPNewRes)
s.arrow(NTWPROPNewRes$c1, lab = names(NTWPROPNewRes$tab),xax = 1,yax = 2)

# Model formula and function ----
#############
GAM.formula = {
  response ~  s(HDNVAR, bs = "cr", k = 3) + s(TempMean, bs = "cr", k = 3) + 
    s(TempSeason,bs = "cr", k = 3) + s(Precip, bs = "cr", k = 3) +  
    s(precipCV, bs = "cr", k = 3) + s(Evenness, bs = "cr", k = 3) + 
    s(shannon_index, bs = "cr", k = 3) + s(NPP, bs = "cr", k = 3) + 
    s(HF, bs = "cr", k = 3
    )
}

GAM.models.autoc <- function(res.dir,               # directory to save results
                             data.M,                # Data for the model (predictors)
                             NETVARS,               # Dataframe with network properties to test
                             networkprops,          # index of the properties in test of NETVARS
                             GAM.form,              # model formula
                             xypoints,              # pixels coordinates
                             nbs_dist,              # neighbours distance for autocovariate model
                             ncores = detectCores()){
  
  # Function to build my GAM models with random forest residuals from a autocovariate variables and landscape descriptors
  dir.create(res.dir)
  library(mgcv)
  if(parallel){
    NC <- ncores
    library(parallel)
    library(spdep)
    library(vegan)
    library(ade4)
    library(ape)
    library(randomForest)
    
    fun.to.par <- function(i){
      print(Sys.time(), paste0(" Starting autocovariate in ", i))
      auto <-  autocov_dist(z = networkprops[,i],xy = xypoints,longlat = FALSE,nbs = nbs_dist,style = "B")     # create the autocavariate variable
      save(auto, file =  paste0(res.dir,"autocovariate",names(networkprops)[i],"mod.Rdata"))
      print(Sys.time(), paste0(" Finished autocovariate in ", i))
      
      print(Sys.time(), paste0(" Starting RF Models in ", i))         
      # Random forest model (PCA AXIS ~ predictors)
      RF <- randomForest(auto ~  EnvDataClean$TempMean + EnvDataClean$TempSeason +                 
                           EnvDataClean$Precip + EnvDataClean$precipCV + 
                           EnvDataClean$NPP + EnvDataClean$HF + 
                           EnvDataClean$Evenness + EnvDataClean$shannon_index + 
                           (EnvDataClean$TempMean)^2 + (EnvDataClean$TempSeason)^2 +
                           (EnvDataClean$Precip)^2 + (EnvDataClean$precipCV)^2 +
                           (EnvDataClean$NPP)^2 + (EnvDataClean$HF)^2 + 
                           (EnvDataClean$shannon_index)^2 + 
                           (EnvDataClean$Evenness)^2
      )
      save(RF,file =  paste0(res.dir,"RandomForest",names(networkprops)[i],"mod.Rdata"))
      print(Sys.time(), paste0(" Finished RF Models in ", i))
      
      # Spatial residuals variable 
      RF.resid <- auto -  RF$predicted #   obs minus fitted
      data.M$response <- networkprops[,i]
      
      GAM.formula <- GAM.form
      data.M$HDNVAR <- RF.resid
      
      print(Sys.time(), paste0(" Starting GAM Models in ", i))
      # Running the gam.model PCA axis ~ predictors 
      gam.model <- gam(formula = GAM.formula ,data = data.M,family = "gaussian") 
      save(gam.model,file = paste0(res.dir,"VARS",names(networkprops)[i],"CV","GAMmod.Rdata"))
      print(Sys.time(), paste0(" Finished GAM Models in ", i))
      
      
      png(filename = paste0(res.dir,"VARS",names(networkprops)[i],"CV",".png"),width = 1750,height = 500)
      par(mfrow=c(1,9))
      plot(gam.model, residuals=F)
      dev.off()
      
      type <- "deviance"  
      resid <- residuals(gam.model, type = type)
      linpred <- napredict(gam.model$na.action, gam.model$linear.predictors)
      observed.y <- napredict(gam.model$na.action, gam.model$y)
      
      png(filename = paste0(res.dir,"VARS",names(networkprops)[i],"CV","RESID.png"),height=1000,width=1000)
      par(mfrow=c(2,2))
      qq.gam(gam.model, rep = 0, level = 0.9, type = type, rl.col = 2,
             rep.col = "gray80")
      hist(resid, xlab = "Residuals", main = "Histogram of residuals")
      plot(linpred, resid, main = "Resids vs. linear pred.",
           xlab = "linear predictor", ylab = "residuals")
      plot(fitted(gam.model), observed.y, xlab = "Fitted Values",
           ylab = "Response", main = "Response vs. Fitted Values")
      dev.off()
      
      resid.tab <- data.frame(PAGENAME = data.M$PAGENAME, Resid = residuals(gam.model),row.names =data.M$PAGENAME)
      dbbf <- data.frame(PAGENAME = template$PageName,VAL = NA, row.names = template$PageName)
      dbbf[rownames(resid.tab),2] <- resid.tab$Resid
      
      # create and save raster with GAM residuals 
      ResidRaster <- fun.dbf2raster(SPPPA = dbbf,mask.dir = "mask/")    
      writeRaster(x = ResidRaster,filename =  paste0(res.dir,"VARS",names(networkprops)[i],"CV","RESID.img"), overwrite=TRUE)
    }
    
    message(Sys.time()," Starting Models")
    mclapply(X = NETVARS ,FUN = fun.to.par,mc.cores = NC)
    message(Sys.time(), " Finished Models")
  }
}


# New results with new model ----
# PCA Axis 1
GAM.models.autoc(res.dir = "models/Revision/", # Output location
                 GAM.form = GAM.formula,
                 data.M = EnvDataClean,
                 NETVARS = 1,            # PCA Axis 1
                 networkprops = NTWPROPNewRes$li,
                 xypoints = geopoints,
                 nbs_dist = 110000,     # neighbours distance 110 000 km
                 ncores = 1)   # not parallel

GAM.models.autoc(res.dir = "models/Revision/", # Output location
                 GAM.form = GAM.formula,
                 data.M = EnvDataClean,
                 NETVARS = 2,            # PCA Axis 2
                 networkprops = NTWPROPNewRes$li,
                 xypoints = geopoints,
                 nbs_dist = 110000,     # neighbours distance 110 000 km
                 ncores = 1)          # not parallel




# Plots and summary: PCA axis 1
load("models/Revision/VARSAxis1CVGAMmod.Rdata") # Output
summary(gam.model)

png(filename = paste0("../../../R/Alpha Networks/MS/Review 1/models/","VARS","axis1","CV",".png"),width = 1750,height = 500)
par(mfrow=c(1,9))
plot(gam.model, residuals=F)
dev.off()


# Plots and summary: PCA axis 2
load("models/Revision/VARSAxis2CVGAMmod.Rdata") # Output
summary(gam.model)

png(filename = paste0("../../../R/Alpha Networks/MS/Review 1/models/","VARS","axis2","CV",".png"),width = 1750,height = 500)
par(mfrow=c(1,9))
plot(gam.model, residuals=F)
dev.off()

# End

