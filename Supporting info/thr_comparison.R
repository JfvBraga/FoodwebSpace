#
# Spatial Network ----
# 10 KM
# Joao Braga
# Comparison between three different habitat threshold (0,10 and 20%)

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

# As table
bioregions <- shapefile(x = "../envData/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp")  # shape file containing metadata
regionsCode <- bioregions@data[-10,1:4]

lbls <- regionsCode$code
names(lbls) <- c(1:9,11,12)

Bioreg <- read.table(file="../envData/bio/10k/bioregion.txt",header = TRUE,sep = ';')   # Bioregions in pixel by region format, @ 10km resolution (from ArcGis)
rownames(Bioreg) <- as.character(Bioreg$PAGENAME)  # Pixel ID


# Loading properties ----
Thr20 <- get(load(file = '../BioticData/10k/Network_properties_Thr20.Rdata'))
Thr10 <- get(load(file = '../BioticData/10k/Network_properties_Thr10.Rdata'))
Thr0 <- get(load(file = '../BioticData/10k/Network_properties_with_basal_decomposition.Rdata'))

# Basic stats ----
# Same pixel selection ----
Thr10 <- Thr10[rownames(Thr0),]
Thr20 <- Thr20[rownames(Thr0),]

summary(Thr0[,-1])
summary(Thr10[,-1])
summary(Thr20[,-1])

# Boxplots ----
boxplot(Thr0$S, Thr10$S, Thr20$S)
boxplot(Thr0$L, Thr10$L, Thr20$L)
boxplot(Thr0$C, Thr10$C, Thr20$C)
boxplot(Thr0$propI, Thr10$propI, Thr20$propI)


# Analysing differences
S_sensitive <- data.frame(PageName=Thr0$Pagename, S0 = Thr0$S, S10 = Thr10$S, S20 = Thr20$S,BioRegion= Bioreg[Thr0$Pagename,'MAX'],row.names = Thr0$Pagename)
head(S_sensitive)

S_sensitive <- reshape2::melt(S_sensitive, id.vars = c("PageName", "BioRegion"), variable.name = "HabThresh", value.name = "S")
names(S_sensitive)[grepl("variable", names(S_sensitive))] <- "HabThresh"
names(S_sensitive)[grepl("value", names(S_sensitive))] <- "S"               
S_sensitive$HabThresh <- sub("S", "", S_sensitive$HabThresh)

x11()
library(Hmisc)
S_sen_plot <- ggplot(na.omit(S_sensitive), aes(x = HabThresh, y = S, colour = as.factor(BioRegion), group = BioRegion)) +
  stat_summary(fun.y = "mean", geom = "line") +
  stat_summary(fun.data = "mean_sdl") +
  theme_bw()+
  scale_colour_manual(name = "Biogeographic regions",values = c("#a6cee3",
                                 "#1f78b4",
                                 "#b2df8a",
                                 "#33a02c",
                                 "#fb9a99",
                                 "#e31a1c",
                                 "#fdbf6f",
                                 "#ff7f00",
                                 "#cab2d6",
                                 "#6a3d9a",
                                 "#b15928"),
                      labels = lbls) + 
  theme(text = element_text(size = 10), legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.key.width = unit(x = 0.06, units = "npc"), legend.position = "",plot.title = element_text(hjust = 0.5))+
  labs(x = "Habitat threshold (%)", y = "Species richness", fill = "")

# Connectance ----
C_sensitive <- data.frame(PageName=Thr0$Pagename, C0 = Thr0$C, C10 = Thr10$C, C20 = Thr20$C,BioRegion= Bioreg[Thr0$Pagename,'MAX'],row.names = Thr0$Pagename)
head(C_sensitive)

C_sensitive <- reshape2::melt(C_sensitive, id.vars = c("PageName", "BioRegion"), variable.name = "HabThresh", value.name = "C")
names(C_sensitive)[grepl("variable", names(C_sensitive))] <- "HabThresh"
names(C_sensitive)[grepl("value", names(C_sensitive))] <- "C"               
C_sensitive$HabThresh <- sub("C", "", C_sensitive$HabThresh)

x11()
library(Hmisc)
C_sen_plot<- ggplot(na.omit(C_sensitive), aes(x = HabThresh, y = C, colour = as.factor(BioRegion), group = BioRegion)) +
  stat_summary(fun.y = "mean", geom = "line") +
  stat_summary(fun.data = "mean_sdl") +
  theme_bw()+
  scale_colour_manual(name = "Biogeographic regions",values = c("#a6cee3",
                                                                "#1f78b4",
                                                                "#b2df8a",
                                                                "#33a02c",
                                                                "#fb9a99",
                                                                "#e31a1c",
                                                                "#fdbf6f",
                                                                "#ff7f00",
                                                                "#cab2d6",
                                                                "#6a3d9a",
                                                                "#b15928"),
                      labels = lbls) + 
  theme(text = element_text(size = 10), legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.key.width = unit(x = 0.06, units = "npc"), legend.position = "none",plot.title = element_text(hjust = 0.5))+
  labs(x = "Habitat threshold (%)", y = "Connectance", fill = "")

library(gridExtra)
library(ggpubr)
X11()
leg <- as_ggplot(get_legend(C_sen_plot))
pdf(file = '../../R/Alpha Networks/MS/Review 1/Sensitive.pdf',height = 4,width = 9)
grid.arrange(S_sen_plot,C_sen_plot,leg,ncol=3)
dev.off()
