# Analysis for a selection of metrics (direct metrics instead of a pca)
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
# Model ----
#############
GAM.formula = {
  response ~  s(HDNVAR, bs = "cr", k = 3) + s(TempMean, bs = "cr", k = 3) + 
    s(TempSeason,bs = "cr", k = 3) + s(Precip, bs = "cr", k = 3) +  
    s(precipCV, bs = "cr", k = 3) + s(Evenness, bs = "cr", k = 3) + 
    s(shannon_index, bs = "cr", k = 3) + s(NPP, bs = "cr", k = 3) + 
    s(HF, bs = "cr", k = 3
    )
}

GAM.models.autoc <- function(res.dir,data.M,NETVARS, 
                             networkprops,plots =TRUE, 
                             spatial.resid=TRUE, 
                             GAM.form,xypoints, 
                             parallel= TRUE, 
                             nbs_dist, ncores = detectCores())
{
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
      # print(Sys.time(), paste0(" Starting autocovariate in ", i))
      auto <-  autocov_dist(z = networkprops[,i],xy = xypoints,longlat = FALSE,nbs = nbs_dist,style = "B")
      save(auto, file =  paste0(res.dir,"autocovariate",names(networkprops)[i],"mod.Rdata"))
      # print(Sys.time(), paste0(" Finished autocovariate in ", i))
      
      # print(Sys.time(), paste0(" Starting RF Models in ", i))
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
      # print(Sys.time(), paste0(" Finished RF Models in ", i))
      
      RF.resid <- auto -  RF$predicted #   obs minus fitted
      data.M$response <- networkprops[,i]
      
      GAM.formula <- GAM.form
      data.M$HDNVAR <- RF.resid
      
      # print(Sys.time(), paste0(" Starting GAM Models in ", i))
      gam.model <- gam(formula = GAM.formula ,data = data.M,family = "gaussian") 
      save(gam.model,file = paste0(res.dir,"VARS",names(networkprops)[i],"CV","GAMmod.Rdata"))
      # print(Sys.time(), paste0(" Finished GAM Models in ", i))
      
      
      png(filename = paste0(res.dir,"VARS",names(networkprops)[i],"CV",".png"),width = 1750,height = 500)
      par(mfrow=c(1,9))
      plot(gam.model, residuals=F)
      dev.off()
      
    }
    message(Sys.time()," Starting Models")
    mclapply(X = NETVARS ,FUN = fun.to.par,mc.cores = NC)
    message(Sys.time(), " Finished Models")
  }
}


# Load data
load("../10k/Network_properties_with_basal_decomposition.Rdata") # local food web properties
load("../10k/Environmental_table.Rdata")                         # Local environmental data

# keep same number of pixels and same order
EnvDataClean <- EnvDataClean[rownames(networkprops),]            
EnvDataClean <- EnvDataClean[complete.cases(EnvDataClean),]
networkprops <- networkprops[rownames(EnvDataClean),]

# Extracting coordinates (X and Y) needed to the models
griid <- shapefile(x = "grid/grid_10Km.shp")
template <- read.dbf(file = "./grid/grid_10Km.dbf")
coords <- cbind(EnvDataClean$X, EnvDataClean$Y)           # Extracting spatial X and Y coordinates
geopoints <- SpatialPoints(coords = coords,proj4string = crs(griid)) # converting coordinates into spatial object (to be used by the autocov model)


# Running the PCA on the significantly different (from null model) variables
NTWPROPNewRes <- dudi.pca(df = networkprops[,c("S","char.path.length","LinkDensity","propOmn","C","propI","clust.coef","meanTL" )],scannf = F) # Characteristic path length with infinite values (removed)
summary(NTWPROPNewRes)
s.arrow(NTWPROPNewRes$c1, lab = names(NTWPROPNewRes$tab),xax = 1,yax = 2)

# Vaiable selection (from PCA)
# C
# S
# Mean Trophic level
# Link density
# Char. path lengh

var.for.model <- c("C","S", "meanTL","LinkDensity","char.path.length","LinkDensity")

GAM.models.autoc(res.dir = "models/SingleMetricModels/", # Output location
                 GAM.form = GAM.formula,
                 data.M = EnvDataClean,
                 NETVARS = which(colnames(networkprops) %in% var.for.model),            # Axis
                 networkprops = networkprops,
                 xypoints = geopoints,
                 nbs_dist = 110000,     # neighbours distance 110 000 km
                 ncores = 1)

#####
# Variable importance
models <- list.files(path = "models/SingleMetricModels/",pattern = "GAMmod.Rdata",full.names = T)
RFmodels <- list.files(path = "models/SingleMetricModels/",pattern = "^Random",full.names = T)[c(2,1,3:5)]
modelsNAMES <- list.files(path = "models/SingleMetricModels/",pattern = "GAMmod.Rdata",full.names = FALSE)
modelsNAMES <- gsub(pattern = "VARS",replacement = "",x = modelsNAMES)
modelsNAMES <- gsub(pattern = "models/SingleMetricModels/",replacement = "",x = modelsNAMES)
modelsNAMES <- gsub(pattern = "CVGAMmod.Rdata",replacement = "",x = modelsNAMES)


varnames <- c("HDNVAR","TempMean","TempSeason","Precip","precipCV","Evenness","shannon_index","NPP","HF")

for(i in 1:5){
  print(paste("Loading",models[i]))  
  load(models[i])
  print(paste("Loading",RFmodels[i]))  
  load(RFmodels[i])
  
  data.M <- EnvDataClean
  data.M$HDNVAR <- RF$y-RF$predicted
  fitted_values <- fitted(gam.model)
  data.M <- data.M[varnames]
  
  fun.to.par <- function(vari){
    for(iii in 1:1000){
      
      fitted_values
      
      rand <-  data.M
      rand[,vari] <- sample(rand[,vari], length(rand[,vari]),replace = FALSE)
      preds <- predict.gam(object = gam.model,newdata = rand)
      res <- cor(fitted_values,preds,method = "spearman")
      write(res,file = paste0("models/SingleMetricModels/variable_importance/",modelsNAMES[i],"_",vari,"_",iii,".txt"))
    }
    
  }
  message(Sys.time()," Starting randomization for ",modelsNAMES[i])
  mclapply(X = varnames,FUN = fun.to.par,mc.cores = 3)
  message(Sys.time()," Finished randomization for ",modelsNAMES[i])
}

# Results
outputs <-   function(model,var) 1-mean(as.numeric(unlist(lapply(X = list.files(path = "./models/SingleMetricModels/variable_importance/",pattern=paste0(model,'_',var,'_'),full.names = T), FUN = function(x) readLines(x)))))
model_results <- function(model, varnames) {
  outp <- lapply(varnames,FUN = function(x) outputs(model,var = x))
  names(outp) <- varnames
  return(outp)
}
res <- lapply(X = modelsNAMES,FUN = function(x) model_results(model=x,varnames=varnames))
names(res) <- modelsNAMES

varImportance <- do.call(rbind, lapply(res, data.frame, stringsAsFactors=F))
#                       HDNVAR  TempMean TempSeason      Precip    precipCV    Evenness shannon_index       NPP          HF
# C                0.080544150 0.1674344 0.09684329 0.064831686 0.061373413 0.018465934    0.13935971 0.1151695 0.393738227
# char.path.length 0.001514634 0.1598480 0.05808640 0.196368423 0.019472041 0.006026897    0.14713597 0.1304480 0.341016307
# LinkDensity      0.115167426 0.4633469 0.52985955 0.009406719 0.007160971 0.037878657    0.10384647 0.2526503 0.013354312
# meanTL           0.024912207 0.5207503 0.39799796 0.001168098 0.021661556 0.038447005    0.07333707 0.2243299 0.005990797
# S                0.073121169 0.2995137 0.53267375 0.021258391 0.013853662 0.045092091    0.14672798 0.1396241 0.020238299


# Fig smoother effects
# for food web metrics
lbls <- c("Spatial residuals",
          "Annual average temp.",
          "Temperature season.",
          "Precipitation",
          "Coeff. of var. of precip.",
          "Evenness index",
          "Shannon-Weiver index",
          "Net primary prod.",
          "Human footprint")

xlabs <- c("Spatial residuals",
           "°C x 100",
           "°C x 100",
           expression(paste("mm of precipitation x year"^-1)),
           "Coefficient of Variation",
           "J'",
           "H'",
           expression(paste("g of Carbon x year"^-1)),
           "Human footprint index")

lbls_mod <- c(
          "Connectance",
          "Characteristic path length",
          "Link density",
          "Mean trophic level",
          "Species richness")

lblsy <- do.call("expression", lapply(1:9, function(i) substitute("F"[X] , list(X = lbls[i]))))

pdf(file = "models/SingleMetricModels/PartialEffects12.pdf",width = 15.75,height = 2.5*2)
par(mfrow=c(2,9))
ii <- 1
for(mod in models[1:2]){
  load(mod)
  for(i in c(2:9,1)) {
    plot.gam(gam.model,select = i, xlab = xlabs[i],
             ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
    if(i == 2) mtext(lbls_mod[ii], side=3, adj=0, line=3, cex=1, font=2) 
    
    
  }
  ii <- ii+1
}
dev.off()

pdf(file = "models/SingleMetricModels/PartialEffects22.pdf",width = 15.75,height = 2.5*3)
par(mfrow=c(3,9))
ii <- 3
for(mod in models[3:5]){
  load(mod)
  for(i in c(2:9,1)) {
    plot.gam(gam.model,select = i, xlab = xlabs[i],
             ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
    if(i == 2) mtext(lbls_mod[ii], side=3, adj=0, line=3, cex=1, font=2) 
    
    
  }
  ii <- ii+1
}
dev.off()
