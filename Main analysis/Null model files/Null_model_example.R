#################################################
# Spatial network null model and examples
# Joao Braga 9 May 2019
# 
# I strongly advise to parallelize the whole process, since the number of iterations can be very high
#################################################
rm(list = ls())


# Setting WD
setwd("(YOUR directory with data here)")

# Requires
library(igraph)    # to generate a random network
library(plyr)
library(NetIndices)
library(parallel)


source('Functions_for_dray.R')


results.path <- "./"


# E.g. of data
######

# For this example I will generate some data, but I will conserve the original names of the European metaweb
# Let BARMdiet.binary be an adjacency matrix of 300 species, with a connetance 0.05
# This metaweb is composed by 50% of birds, 25% reptiles, 20% of mammals, and 5% of amphibian
Spp_names <- c(paste0("B",seq(1,300*0.50,length.out = 300*0.50)),
               paste0("R",seq(1,300*0.25,length.out = 300*0.25)),
               paste0("M",seq(1,300*0.20,length.out = 300*0.20)),
               paste0("A",seq(1,300*0.05,length.out = 300*0.05))
               )

BARMdiet.binary <-matrix(rbinom(300*300,size = 1,prob = 0.05),nrow = 300,ncol = 300,
                         dimnames = list(Spp_names,Spp_names)
                           )

# Let us add some diet categories, they will help up us define who are the basal species later.
# In the original metaweb, aprox. 70% of species eat diet categories, so lets define the categories as follows:
dietcat_names <- paste0("diet",1:12)  # naming the diet categories (12 like in the original data)

dietcat_adj   <- matrix(rbinom(300*12,size = 1,prob = 0.70),nrow = 12,ncol = 300,dimnames = list(dietcat_names,Spp_names))  # defining who eats diet categories with 70% of prob.

BARMdiet.binary <- rbind(dietcat_adj, BARMdiet.binary) # adding the dietcategories to the original adjacency matrix

BARMdiet.binary <- cbind(matrix(0,nrow = nrow(BARMdiet.binary),ncol = length(dietcat_names),dimnames = list(rownames(BARMdiet.binary),dietcat_names)),BARMdiet.binary) # adding 12 columns of zero that correspond to the dietcategories as predators (zeros because they do not eat)

dim(BARMdiet.binary) # now it is a square matrix


plot.graph.trophic.level(BARMdiet.binary)

# image(BARMdiet.binary) # there's a stripe of zero, that correspond to the diet categories "prey".

# Now, lets generate some presence and absence for species
# Let say that we have 500 pixels
# Let us also generate some gradient within our pixels, e.g. elevation. The higher it is, the less likely it is
# to find species over the top. Lets define the probability of a species being present as:
master <- t(apply(X = matrix(seq(0.95,0.0,length.out = 500),nrow = 500,ncol = 1),MARGIN = 1,FUN = function(x) rbinom(300,size = 1,prob = (0.05 + x))))
colnames(master) <- Spp_names
rownames(master) <- paste0("pix",1:500)
master <- as.data.frame(master)
head(master)

# In my GEB paper, I used different constraints in my null model:
# 1) Species need to have at least one prey (being diet category or vertebrate prey);
# 2) No disconnected nodes or portions of network;
# 3) Conserve observed proportion of the different taxonomical groups.

# Note that in the GEB I had more 70,000 pixels, so instead of doing it by pixel, I ran my null model by species richness values.
# E.g. In 70,000 pixels, I had only 300 different values of species richness. So, I ran my null model 1000 times for each value of richness, so 300,000 times
# instead of > 70,000,000 times. Then, I compared the network property value of pixel i with species richness j, with the correspondent null distribution of species 
# richness j.

# To meet constrain 1) I made sure that everytime I generated a random network all predators (except diet categories) add prey. If sampled species had no prey,
# I continued resampling until 1) was met.
# To meet constrain 2), certain metrics such as  characteristic path length would be Inf or NA for certain species, this would be mean that those species are disconnected (even if they meet 1)
# To meet constrain 3), I generated a distribution of "observed" proportion of taxa from all pixels of species richness j. Then, I sampled from this distribution 1000 times to obtained
# networks with real observed proportions of taxa. 

# For Everytime any of these constrain were not meet, I re-did them.


# How I calculated the observed Proportion of taxa
taxa_type <- substr(colnames(master), 1, 1)
taxa_S_pix_numeric <- apply(X = master,MARGIN = 1, 
                            FUN = function(x) {
                              setNames(object = as.vector(table(taxa_type[as.logical(x)])),
                                       names(table(taxa_type[as.logical(x)])))
                            }
)

taxa_type_df <- data.frame(matrix(NA,ncol = 4,nrow = length(taxa_S_pix_numeric)),row.names = names(taxa_S_pix_numeric)) ; colnames(taxa_type_df) <- c("A","B","M","R")

for(x in names(taxa_S_pix_numeric)){
  taxa_type_df[x,names(taxa_S_pix_numeric[[x]])] <-  taxa_S_pix_numeric[[x]]
} 
taxa_type_df[is.na(taxa_type_df)] <- 0 ; rm(taxa_S_pix_numeric)

head(taxa_type_df,1)
#head(master,1)


# NULL model
#####################################
# Values needed
richness_pix <- rowSums(taxa_type_df)          # Vector with richness values of all pixels
length(richness_pix)

S_richness    <- unique(richness_pix)          # unique values of species richness
length(S_richness)

dietcat <- head(row.names(BARMdiet.binary),12)  # Vector containing the diet categories

# Total iterations
replicates <- rep(S_richness,1000)

# E.g. for RUN replicates[75]
# Function to generate the random network
# The output should be an adjacency matrix of a network that meets all criteria
randWEB <- random_network(dietcat = dietcat,
               metaweb = BARMdiet.binary,
               y = replicates[75],
               S_richness =  richness_pix,
               prop.taxa = taxa_type_df)

iweb <- graph.adjacency(randWEB, mode="directed")  
char.path.length <- mean(shortest.paths(iweb))   # Are all nodes should be somehow reachable, if NA, then re-run

                                                                                                                          
plot.graph.trophic.level(randWEB) # Ploting the web

# E.g. for RUN replicates[242]  --> Smaller webs take more time to converge!!
# Function to generate the random network
# The output should be an adjacency matrix of a network that meets all criteria
randWEB <- random_network(dietcat = dietcat,
                          metaweb = BARMdiet.binary,
                          y = replicates[242],
                          S_richness =  richness_pix,
                          prop.taxa = taxa_type_df)

iweb <- graph.adjacency(randWEB, mode="directed")  
char.path.length <- mean(shortest.paths(iweb))   # Are all nodes should be somehow reachable, if NA, then re-run


plot.graph.trophic.level(randWEB) # Ploting the web


# Then, to create a thousand of these random webs, just use some embarrassingly parallel computations. In my case I did as follow:
# RUN
replicates <- rep(S_richness,1000)

NCOREs <- 48

fun.to.par <- function(i) {
  results.dir <- paste0(results.path,i)       # a folder for each value of richness to put the results inside
  dir.create(results.dir, recursive = TRUE)  
  
  randWEB <- random_network(dietcat = dietcat,
                            metaweb = BARMdiet.binary,
                            y = i,
                            S_richness =  richness_pix,
                            prop.taxa = taxa_type_df)
  
  res <- fun.net.prop(randWEB,dietcat = dietcat)  # A custom function that I used to compute all metrics in GEB paper
  
  write.table(res, file = paste0(results.dir,"/Rand",i,"_",round(runif(n = 1,min = 0,max = 500000)),".txt")) # save a file with metrics per network
  
}

date()
mclapply(X =  replicates,FUN = fun.to.par,mc.cores = NCOREs)    # from parallel package.
date()

# Then compile results for each value of richness into a table
fun.to.par <- function(y){
  simsY <-  list.files(path = paste0(results.path,y,"/"), pattern = ".txt$",full.names = TRUE)
  SimTab <- compile_results(files.loc.vector =simsY,y)
  save(SimTab,file = paste0(results.path,y,".Rdata"))
}

date()
mclapply(X =  S_richness,FUN = fun.to.par,mc.cores = NCOREs)    # from parallel package.
date()

#####
# End
