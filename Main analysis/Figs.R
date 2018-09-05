# Analysis
rm(list = ls())
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
load("../10k/Environmental_table.Rdata")  
# keep same number of pixels and same order
EnvDataClean <- EnvDataClean[rownames(networkprops),]            
EnvDataClean <- EnvDataClean[complete.cases(EnvDataClean),]
networkprops <- networkprops[rownames(EnvDataClean),]

# Data Working directory
setwd("~/Documents/FACULDADE/Universite Joseph Fourier/PhD/FoodWebs/data/BioticData/10k/")

# Running the PCA on the significantly different (from null model) variables

NTWPROPNewRes <- dudi.pca(df = networkprops[,c("S","char.path.length","LinkDensity","propOmn","C","propI","clust.coef","meanTL" )],scannf = F) # Characteristic path length with infinite values (removed)
summary(NTWPROPNewRes)
s.arrow(NTWPROPNewRes$c1, lab = names(NTWPROPNewRes$tab),xax = 1,yax = 2)


# Figs:
#Figure: Partial effects ----
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

lblsy <- do.call("expression", lapply(1:9, function(i) substitute("F"[X] , list(X = lbls[i]))))

pdf(file = "../../../R/Alpha Networks/MS/Review 1/models/partialEffects.pdf",width = 15.75,height = 5)
# X11()
par(mfrow=c(2,9)
    # ,mar=c(5, 0, 5, 0)
)
load("models/Revision/VARSAxis1CVGAMmod.Rdata")
for(i in c(2:9,1)) {
  # plot(x=1,y=1, xlim=c(0,100), ylim=c(-7,4),type="n", xlab="", ylab="")
  # par(new=TRUE)
  
  plot.gam(gam.model,select = i, xlab = xlabs[i],
           ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
  if(i == 2) mtext("(a)", side=3, adj=0, line=3, cex=1, font=2) 
  # axis(side = 1)
  # axis(side = 2)
  
}
load("models/Revision/VARSAxis2CVGAMmod.Rdata")
for(i in c(2:9,1)) {
  # plot(x=1,y=1, xlim=c(0,100), ylim=c(-7,4),type="n", xlab="", ylab="")
  # par(new=TRUE)
  
  plot.gam(gam.model,select = i, xlab = xlabs[i],
           ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
  if(i == 2) mtext("(b)", side=3, adj=0, line=3, cex=1, font=2) 
  # axis(side = 1)
  # axis(side = 2)
  
}

dev.off()


# PCA plot ----
LBLS <- c("Sp rich.","Char. path. length","Link dens.","Prop. omn.","Connect.","Prop. int.","Cluster coef.","Mean TL" )
png(filename = "models/Revision/Fig2PCA.png",width = 6,height = 6,units = "in",res = 300)
par(mar=c(5, 5, 5, 5) + 0.1)

plot(x = c(-15,15),y = c(-15,15),type="n",xlab="PC 1, 49.17%", ylab="PC 2, 27.47%",axes=F)
axis(1, xlim=c(-15,15),lwd=2,line=1)
axis(2, xlim=c(-15,15),lwd=2,line=1)
abline(h = 0,v = 0)

points(x = NTWPROPNewRes$li[,1],y = NTWPROPNewRes$li[,2], pch=".",cex=2,col= rgb(red = 0,green = 0,blue = 0,alpha = 0.2))

par(new=TRUE)
plot(x = c(-0.6,0.6),y = c(-0.6,0.6),type="n", xlab="", ylab="",axes =F)
axis(3, xlim=c(-0.6,0.6),lwd=2,line=1,col = "red")
# mtext(text  = "Correlation with PC 1",side = 3,line = 3, col="red")
axis(side = 4, ylim=c(-0.6,0.6),lwd=2,line=1, col= "red")
# mtext(text  = "Correlation with PC 2",side = 4,line = 3, col="red")
marrows(NTWPROPNewRes$c1,xax = 1,yax = 2, lab = LBLS,grid = FALSE,addaxes = FALSE,origin = c(0,0),boxes = TRUE,add.plot = TRUE,xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), clabel = 0.8)
dev.off()
# Selection of 4 networks (1 from each quadrant of PCA)
# + +
POSPOS <- NTWPROP$li[NTWPROP$li$Axis1>=0 & NTWPROP$li$Axis2>=0,1:2]
mid <- c(mean(POSPOS$Axis1),mean(POSPOS$Axis2))

diisst <- dist(rbind(mid,POSPOS),method = "euclidean")
diisstM <- as.matrix(diisst)
which(diisstM[,1]==min(diisstM[-1,1]))
# FW404 
# 1254 

# + -
POSNEG <- NTWPROP$li[NTWPROP$li$Axis1>=0 & NTWPROP$li$Axis2 < 0,1:2]
mid <- c(mean(POSNEG$Axis1),mean(POSNEG$Axis2))
diisst <- dist(rbind(mid,POSNEG),method = "euclidean")
diisstM <- as.matrix(diisst)
which(diisstM[,1]==min(diisstM[-1,1]))
# QP518 
# 9291 

# - -
NEGNEG <- NTWPROP$li[NTWPROP$li$Axis1 < 0 & NTWPROP$li$Axis2 < 0,1:2]
mid <- c(mean(NEGNEG$Axis1),mean(NEGNEG$Axis2))
diisst <- dist(rbind(mid,NEGNEG),method = "euclidean")
diisstM <- as.matrix(diisst)
which(diisstM[,1]==min(diisstM[-1,1]))
# QZ198 
# 16730 

# - +
NEGPOS <- NTWPROP$li[NTWPROP$li$Axis1 < 0 & NTWPROP$li$Axis2 >= 0,1:2]
mid <- c(mean(NEGPOS$Axis1),mean(NEGPOS$Axis2))
diisst <- dist(rbind(mid,NEGPOS),method = "euclidean")
diisstM <- as.matrix(diisst)
which(diisstM[,1]==min(diisstM[-1,1]))
# KT530 
# 7809 

# Plot food webs from mid pixel ----
to.plot <- c("FW404","QP518","QZ198","KT530")

# Loading species presence
load("../../BioticData/feeding/Spp_traits_habs/Species_presence/MASTER.bin10000_allhab_tresh0.Rdata")
master$PAGENAME <- as.character(master$PAGENAME)
rownames(master) <- as.character(master$PAGENAME)
xx <- colnames(master)[-1]

pix.hab <-  read.dbf(file = "./tabulateGLC_10Km.dbf")
pix.hab$PAGENAME <- as.character(pix.hab$PAGENAME)
rownames(pix.hab) <- as.character(pix.hab$PAGENAME)

colnames(pix.hab) <- gsub(pattern = "VALUE_",replacement = "X",x = colnames(pix.hab))
# Loading habitat per species
BARM.HAB <- read.csv("../../BioticData/feeding/Spp_traits_habs/Habitats/BARM_allhabs.csv",header = TRUE,row.names = 'ID',sep = ';')
load("../../BioticData/feeding/Spp_traits_habs/Spp_interactions/All/BARMdiet_bin.RData")
dietcat <- colnames(BARMdiet.binary)[1:12]

#
# To plot 
for(y in to.plot){
  # Generate the web
  if(sum(master[y,2:ncol(master)]) == 0) stop() # if a pixel in the master matrix is NA, then NA in the interaction map
  else{
    dddd <- xx[which( master[y,2:length(master)]==1, arr.ind=FALSE)] # if not NA, then which species are present in that pixel
    if(length(dddd)!= 0) dddd <- c(dietcat, dddd)                     # If spp are found in that pixel, then add the diet categories to the analysis
    hab <- pix.hab[y,]
    web <- sub_web(metaweb = BARMdiet.binary,SPPCODE = dddd,dietcat = dietcat,PIX.HAB = hab,SPP.HAB = BARM.HAB,HELP = FALSE)             # Select those species from the network to a new network + Diet categories
  }
  save(web,file = paste0("../10k/models/Revision/",y,".Rdata"))
  
  # Calulate trophic level for y coordinates
  community <- Community(data.frame(node = colnames(t(as.matrix(web)))),
                         trophic.links = PredationMatrixToLinks(t(as.matrix(web))),
                         properties = list(title = "TL")
  )
  community <- RemoveCannibalisticLinks(community, title='community');
  # Get all spp trophic levels
  temp <- PreyAveragedTrophicLevel(community)
  
  # Using Igraph to plot by torphic level
  t <- graph.adjacency(adjmatrix = t(as.matrix(web)),mode = "directed")
  dt <- rownames(web)[rownames(web) %in% dietcat]
  ff <- rownames(web)[rownames(web)!=(dt)]
  temp <- temp[ff]
  t <- delete_vertices(graph = t,v = dt)
  rk <- rank(degree(graph = t))
  x <- quantile(rk)
  lay<-matrix(nrow=length(temp),ncol=2) # create a matrix with one column as runif, the other as trophic level
  lay[,1]<- runif(n = length(rk))
  lay[,2]<- temp
  par(mar=c(.1,.1,.1,.1))
  V(graph = t)$color[grep(pattern = "^A",x =  names(V(t)))] <- "black"
  V(graph = t)$color[grep(pattern = "^B",x =  names(V(t)))] <- "black"
  V(graph = t)$color[grep(pattern = "^R",x =  names(V(t)))] <- "black"
  V(graph = t)$color[grep(pattern = "^M",x =  names(V(t)))] <- "black"
  png(filename = paste0("../10k/models/Revision/",y,".png"))
  plot.igraph(t,layout=lay,vertex.size=2,vertex.label=NA, edge.arrow.size=.5,edge.width=.5,edge.color=rgb(0.1,0.1,0.1,0.1), vertex.color=V(t)$color)
  dev.off()
}

# Plot individual metrics rasters - Figure S1 ----
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


for(i in 2:length(networkprops)) {
  rtr <- fun.dbf2raster(SPPPA = data.frame(PAGENAME = rownames(networkprops), EXPL = networkprops[,i]),mask.dir = "./mask/")
  
  rtrDT  <- as.data.frame(as(rtr, "SpatialPixelsDataFrame"))
  
  brks <- seq(min(rtrDT$reference_grid_10km), max(rtrDT$reference_grid_10km), length.out = 4)
  brks <- floor(brks * 100)/100
  if(i == 2) brks <- round(brks)
  
  propPLOT <-   ggplot() + 
    geom_tile(data = rtrDT, aes(x = x, y = y, fill= reference_grid_10km)) +  
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_fill_distiller(breaks = c(round(min(rtrDT$reference_grid_10km), digits = 2),brks),
                         palette = "Spectral") +
    theme_bw() + theme(text = element_text(size = 14), legend.text = element_text(size = 14),
                       axis.ticks = element_blank(),
                       legend.key.width = unit(x = 0.06, units = "npc"), legend.position = "bottom",plot.title = element_text(hjust = 0.5)) +
    labs(x = "", y = "", title = lbls[i-1], fill = "")
  
  assign(names(networkprops)[i], propPLOT)
  rm(rtr)
  rm(rtrDT)
  rm(propPLOT)
}

pdf(file = "models/Revision/FIG1.pdf",width =12,height = 24)
grid.arrange(S,
             C,
             LinkDensity,
             char.path.length,
             clust.coef,
             propI,
             propOmn,
             ncol=2)

dev.off()




##### Figure 3 Smoother ----
lbls <- c("Spatial residuals",
          "Annual average temp.",
          "Temperature season.",
          "Precipitation",
          "Coeff. of var. of precip.",
          "Habitat evennesss index",
          "Habitat Shannon index",
          "Net primary productivity",
          "Human footprint")

xlabs <- c("Auto covariate var.",
           "°C x 100",
           "°C x 100",
           expression(paste("mm of precipitation x year"^-1)),
           "Coefficient of Variation",
           "J'",
           "H'",
           expression(paste("g of Carbon x year"^-1)),
           "Human footprint index")

lblsy <- do.call("expression", lapply(1:9, function(i) substitute("F"[X] , list(X = lbls[i]))))

# Plot 
pdf(file = "models/Revision/Fig3.pdf",width = 14,height = 5)
par(mfrow=c(2,8)
    # ,mar=c(5, 0, 5, 0)
)
load("models/Revision/VARSAxis1CVGAMmod.Rdata")
for(i in c(2:9,1)) {
  # plot(x=1,y=1, xlim=c(0,100), ylim=c(-7,4),type="n", xlab="", ylab="")
  # par(new=TRUE)
  plot.gam(gam.model,select = i, xlab = xlabs[i],
           ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
  if(i == 2) mtext("a)", side=3, adj=0, line=3, cex=1, font=2) 
  # axis(side = 1)
  # axis(side = 2)
  
}

load("models/Revision/VARSAxis2CVGAMmod.Rdata")
for(i in c(2:9,1)) {
  # plot(x=1,y=1, xlim=c(0,100), ylim=c(-7,4),type="n", xlab="", ylab="")
  # par(new=TRUE)
  plot.gam(gam.model,select = i, xlab = xlabs[i],
           ylab =lblsy[i],rug = FALSE,shade = TRUE, main=lbls[i])
  if(i == 2) mtext("b)", side=3, adj=0, line=3, cex=1, font=2) 
  # axis(side = 1)
  # axis(side = 2)
  
}
dev.off()

# FIG S4.1 Individual food web metrics rasters ----
library(gridExtra)
names(networkprops)
pdf(file = "models/Revision/S4_1.pdf",width =12,height = 18)
grid.arrange(S,
             C,
             LinkDensity,
             clust.coef,
             char.path.length,
             propOmn,
             ncol=2)
dev.off()

pdf(file = "models/Revision/S4_2.pdf",width =12,height = 18)
grid.arrange(propBH,
             propBNonH,
             propI,
             propT,
             MAXSIM,
             ncol=2)

dev.off()

pdf(file = "models/Revision/S4_3.pdf",width =12,height = 18)
grid.arrange(GEN,
             GENSD,
             VUL,
             VULSD,
             meanTL,
             maxTL,
             ncol=2)

dev.off()

# Plot predictors distribution ----
# Supplementarty material
# Supplementarty material
load('../models/Revision/VARSAxis1CVGAMmod.Rdata')
data.M <- gam.model$model

lbls <- c("Spatial residuals",
          "Annual average temp.",
          "Temperature season.",
          "Precipitation",
          "Coeff. of var. of precip.",
          "Habitat evennesss index",
          "Habitat Shannon index",
          "Net primary productivity",
          "Human footprint")


for(i in 2:length(data.M)) {
  rtr <- fun.dbf2raster(SPPPA = data.frame(PAGENAME = rownames(data.M), EXPL = data.M[,i]),mask.dir = "../../10k/mask/")
  rtrDT  <- as.data.frame(as(rtr, "SpatialPixelsDataFrame"))
  
  propPLOT <- 
    ggplot() + 
    geom_tile(data = rtrDT, aes(x = x, y = y, fill= reference_grid_10km)) +    
    scale_fill_distiller(palette = "Spectral") +
    theme_bw() + theme(text = element_text(size = 14), legend.text = element_text(size = 14),
                       legend.key.width = unit(x = 0.06, units = "npc"), legend.position = "bottom",plot.title = element_text(hjust = 0.5)) +
    labs(x = "", y = "", title = lbls[i-1], fill = "") +
    theme(legend.position="bottom") 
  
  assign(names(data.M)[i], propPLOT)
  rm(rtr)
  rm(rtrDT)
  rm(propPLOT)
}

pdf(file = "../../../../R/Alpha Networks/MS/Review 1/Supplementary material/S2_1.pdf",width =12,height = 18)
grid.arrange(TempMean,
             TempSeason, 
             Precip,
             precipCV,
             shannon_index,
             Evenness,
             ncol=2)

dev.off()

pdf(file = "../../../../R/Alpha Networks/MS/Review 1/Supplementary material/S2_2.pdf",width =12,height = 6)
grid.arrange(HF,
             NPP,
             ncol=2)
dev.off()

# end