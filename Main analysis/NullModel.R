############################
#
# Null model
# Author: Joao Braga
############################
rm(list=ls())
# Data Working directory
setwd("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/data/BioticData/10k/random_networks_taxa/")

# Requires
library(raster) 
library(Matrix)
library(rgeos)
library(rgdal)
library(foreign)
library(igraph)
library(mapview)
library(NetIndices)
library(cheddar)

source("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/R/FunctionsHDD.R")


# Testing (with FDR) ----
#################################################
# Gathering the resuts from the NULL models (From randomization scripts)
netpvalues <- matrix(data = NA,nrow = nrow(networkprops),ncol = 19) # Store matrix for pvalues


for(i in sort(unique(networkprops$S))){  # For each value of species richness
  rowsS <- which(networkprops$S==i)    # identify the rows
  obsS <- networkprops[rowsS,]         # extract them and store the observed values
  
  # loading the correct null model simulation for i
  load(file = paste0("./random_networks_taxa/",i,".Rdata"))
  SimTab <- cbind(SimTab, as.numeric(SimTab[,2])/as.numeric(SimTab[,1])) # Link density = S / L
  colnames(SimTab)[19] <- "LinkDensity"
  class(SimTab) <- "numeric"

  # p-value 
  NullMeans <- colMeans(SimTab)
  sim0 <- abs(sweep(SimTab,2,NullMeans))
  obs0 <- abs(sweep(obsS[,-1],2,NullMeans))
  
  for(rw in 1:nrow(obsS)){
    for(ii in 1:19){
      netpvalues[rowsS[rw],ii] <- (sum(sim0[,ii]  >=  obs0[rw,ii]) + 1)/(length(sim0[,ii]) + 1)   # two-sided test
      
    }
  }
  
  print(paste(i,"is done!"))
}

# apply(X = netpvalues,MARGIN = 2,FUN = function(x) sum(x < 0.05)/length(x))

# Names
colnames(netpvalues) <- colnames(SimTab)
rownames(netpvalues) <- rownames(networkprops)

# Running the FDR algorithm ----
apply(X = netpvalues,MARGIN = 2,FUN = function(x) sum(p.adjust(x,method = "fdr") <= 0.05)/length(x))
#         S                L*              C*       clust.coe* char.path.lengt*          propOm*           propB* # For reference. This metrics has been decomposed in non-herbivore and herbivore.
# 0.0000000        0.9992266        0.9992266        0.9964628        0.9983391        0.9577559        0.9565642 
# propBNonH           propBH            propT            propI*             GEN            GENSD              VUL 
# 0.5430998        0.4476957        0.0000000        0.9678225        0.8646212        0.8522472        0.0000000 
# VULSD           MAXSIM               meanTL*            maxTL      LinkDensity* 
# 0.0000000        0.7729192        0.9641204        0.0000000        0.9992266 


# Plots ----
#################################################
# Quantifying the upper and lower quantiles of null distribution
# Loading the null distribution for each value of species richness
sims <- list.files(pattern = ".Rdata$",full.names = TRUE)

AverProp  <- matrix(data = NA,nrow = length(sims),ncol = 19)
UpperProp <- matrix(data = NA,nrow = length(sims),ncol = 19)
LowerProp <- matrix(data = NA,nrow = length(sims),ncol = 19)
SD <-        matrix(data = NA,nrow = length(sims),ncol = 19)

for(i in 1:length(sims)){
  load(file = sims[i])
  SimTab <- cbind(SimTab, as.numeric(SimTab[,2])/as.numeric(SimTab[,1]))     # Link density

  for(ii in 1:19) {
    AverProp[i,ii]  <- median(as.numeric(SimTab[,ii]),na.rm=TRUE)
    UpperProp[i,ii]  <- quantile(as.numeric(SimTab[,ii]), probs = 0.975, na.rm=TRUE)
    LowerProp[i,ii]  <- quantile(as.numeric(SimTab[,ii]), probs = 0.025, na.rm=TRUE)
    SD[i,ii]  <- sd(SimTab[,ii],na.rm = TRUE)
    SD[i,1] <- SimTab[1,1]
  }
}  

AverProp <- as.data.frame(AverProp)
AverProp <- AverProp[with(AverProp, order(V1)),]

LowerProp <- as.data.frame(LowerProp)
LowerProp <- LowerProp[with(LowerProp, order(V1)),]

UpperProp <- as.data.frame(UpperProp)
UpperProp <- UpperProp[with(UpperProp, order(V1)),]

SD <- as.data.frame(SD)
SD <- SD[with(SD, order(V1)),]

# 
lbls <- c("Species richness",
          "Number of links",
          "Connectance",
          "Cluster coefficient",
          "Characteristic path length",
          "Proportion of Omnivore spp",
          "Proportion of basal spp",
          "Proportion of non-herbivore basal spp",
          "Proportion of herbivore basal spp",
          "Proportion of top spp",
          "Proportion of intermediate spp",
          "Generality",
          "SD of generality",
          "Vulnerability",
          "SD of vulnerability",
          "Maximum trophic similarity",
          "Mean trophic level",
          "Maximum trophic level",
          "Link density")

png(filename = "../../../../R/Alpha Networks/MS/Review 1/Figs/S51.png",width = 5,height = 6.66666,units = "in",res = 300)
par(mfrow=c(2,2))
for(i in c(3,19,4,5)) {
  # Lower limit
 ll <- c(min(networkprops[,colnames(AverProp)[i]],na.rm = TRUE),min(LowerProp[,i],na.rm = TRUE))
  # Upper limit
  uu  <-  c(max(networkprops[,colnames(AverProp)[i]],na.rm = TRUE), max(UpperProp[,i],na.rm = TRUE))
  limm <- c(min(ll), max(uu))
  
  
  plot(x = AverProp[,1], AverProp[,i], main=lbls[i], type="n", ylim= limm, xlab="Number of species", ylab=lbls[i])
  points(x = networkprops$S,y = networkprops[,colnames(AverProp)[i]],pch=".", cex = 1)
  
  polygon(c(LowerProp[,1], rev(UpperProp[,1])), c(LowerProp[,i], rev(UpperProp[,i])), border=NA, col= '#d78c3933')
 lines(x = AverProp[,1], y = AverProp[,i], main=colnames(AverProp)[i], col = "red") 
}
dev.off()

png(filename ="../../../../R/Alpha Networks/MS/Review 1/Figs/S52.png",width = 5,height = 10,units = "in",res = 300)
par(mfrow=c(3,2))
for(i in c(6,8:12)) {
  # Lower limit
  ll <- c(min(networkprops[,colnames(AverProp)[i]]),min(LowerProp[,i]))
  # Upper limit
  uu  <-  c(max(networkprops[,colnames(AverProp)[i]]), max(UpperProp[,i]))
  limm <- c(min(ll), max(uu))
  
  
  plot(x = AverProp[,1], AverProp[,i], main=lbls[i], type="n", ylim= limm, xlab="Number of species", ylab=lbls[i])
  points(x = networkprops$S,y = networkprops[,colnames(AverProp)[i]],pch=".", cex = 1)
  
  polygon(c(LowerProp[,1], rev(UpperProp[,1])), c(LowerProp[,i], rev(UpperProp[,i])), border=NA, col= '#d78c3933')
  lines(x = AverProp[,1], y = AverProp[,i], main=colnames(AverProp)[i], col = "red") 
}
dev.off()

png(filename = "../../../../R/Alpha Networks/MS/Review 1/Figs/S53.png",width = 5,height = 10,units = "in",res = 300)
par(mfrow=c(3,2))
for(i in 13:18) {
  # Lower limit
  ll <- c(min(networkprops[,colnames(AverProp)[i]]),min(LowerProp[,i]))
  # Upper limit
  uu  <-  c(max(networkprops[,colnames(AverProp)[i]]), max(UpperProp[,i]))
  limm <- c(min(ll), max(uu))
  
  if(i == 16) limm <- c(0, 1) 
  
  plot(x = AverProp[,1], AverProp[,i], main=lbls[i], type="n", ylim= limm, xlab="Number of species", ylab=lbls[i])
  points(x = networkprops$S,y = networkprops[,colnames(AverProp)[i]],pch=".", cex = 1)
  
  polygon(c(LowerProp[,1], rev(UpperProp[,1])), c(LowerProp[,i], rev(UpperProp[,i])), border=NA, col= '#d78c3933')
  lines(x = AverProp[,1], y = AverProp[,i], main=colnames(AverProp)[i], col = "red") 
}
dev.off()
# end 