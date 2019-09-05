#################################################
# Functions 
# Joao Braga
# 22/10/2015
#################################################
# Function to Identify a spp by the code
whois <- function(SPPCODE = NULL, SPPNAME = NULL) {
  # Function to Identify a spp by the code
  
  if(is.null(SPPCODE) & is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  if(!is.null(SPPCODE) & !is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  
  SppID <- read.table(file = "/feeding/Spp_traits_habs/SppID.txt", header = TRUE)
  
  if(length(SPPCODE) > 1){
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if(length(SPPNAME) > 1){
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }
  
  if(!is.null(SPPCODE))    who <- SppID[grepl(pattern = SPPCODE,x = SppID$ID),]
  if(!is.null(SPPNAME))    who <- SppID[grepl(pattern = SPPNAME,x = SppID$Spp),]
  
  return(who)
}


#################################################
# Function to build local networks
sub_web <- function(metaweb = NULL, SPPCODE = NULL,
                    
                    SPP.HAB = NULL, PIX.HAB = NULL,
                    
                    HELP=TRUE){
  
  # Function to subset the meta-web square matrix according
  
  # a subgroup of species into a sub web square matrix 
  
  if(HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   
                   \n    SPP.HAB: spp x habitat matrix
                   
                   \n    PIX.HAB: pix x habitat matrix
                   
                   \n    SPP.EXT: vector of spp to be primarily extinct
                   
                   \n To avoid this message insert the following argument: HELP=FALSE")
  
  if(is.null(metaweb)) stop("Must provide a meta-web")
  
  if(is.null(SPPCODE)) stop("Must provide species codes")
  
  if(!is.null(SPP.HAB) & is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  
  if(is.null(SPP.HAB) & !is.null(PIX.HAB)) stop("Must provide both spp x habitat and pixel x habitat matrices")
  
  # -------------------------------------------------------
  
  # CALCULATE ORIGINAL LOCAL WEB
  
  # -------------------------------------------------------
  # calculate pixel potential web
  subweb <- metaweb[SPPCODE, SPPCODE]
  
  
  
  # FILTER 1. SPP REMOVAL BASED ON HABITAT OVERLAP
  
  # remove forbidden links, based on habitat overlap
  
  if(!is.null(SPP.HAB)){
    sppHAB <- as.matrix(SPP.HAB[SPPCODE[!SPPCODE %in% dietcat], names(PIX.HAB[, -1])[PIX.HAB[, -1] > 0], drop = FALSE])    # subset habitats by spp and habitats present in pixel
    sppHAB[is.na(sppHAB)] <- 0                           # NAs are considered 0s
    sppHAB[sppHAB > 1] <- 1                              # secondary and optimal habitats treated equally
    
    # adding diet categories with presences in all habitats (makes their addition to the final web easier)
    sppHAB <- rbind(matrix(1, nrow = length(dietcat), ncol = ncol(sppHAB), dimnames = list(dietcat, colnames(sppHAB))),
                            sppHAB)
    # creating a matrix of co-occurrences based on habitats
    sppXspp <- sppHAB %*% t(sppHAB)                      
    sppXspp[sppXspp > 1] <- 1                            # transforming to binary
    
    # multiplying each cell of the co-occurence matrix by the subweb will remove spp interactions that cannot occur in the pixel
    subweb <- sppXspp * subweb
  }
  # FILTER 2. REMOVAL BASED ON DISTRIBUTIONS
  # remove SPP with no resources from metaweb (unconnected nodes + pure cannibals) iteratively
  while(any(rowSums(subweb[!row.names(subweb) %in% dietcat,, drop = FALSE], na.rm = TRUE) == 0)){   
    
    subweb <- rbind(subweb[row.names(subweb) %in% dietcat,],                             # add diet categories predators back to keep the matrix square
                    subweb[rowSums(subweb, na.rm = TRUE) > 0,, drop = FALSE])
    
    
    
    # update codes and subweb, since some species were removed from predator list and have to be removed as prey
    
    SPPCODE <- row.names(subweb)
    subweb <- subweb[SPPCODE, SPPCODE]
    
    # remove pure cannibals
    pu.cannib <- colnames(subweb)[rowSums(subweb) == diag(subweb)]       # species (and diet categories) with sum of rows equal to diag represent pure cannibals
    pu.cannib2 <- pu.cannib[!pu.cannib %in% dietcat]                     # excluding diets from cannibal list
    
    subweb <- subweb[!rownames(subweb) %in% pu.cannib2, !colnames(subweb) %in% pu.cannib2] # remove pure cannibals
  }
  # remove unconnected nodes
  subweb <- subweb[which(!(rowSums(subweb) == 0 & colSums(subweb) == 0)),which(!(rowSums(subweb) == 0 & colSums(subweb) == 0)), drop = FALSE]
  return(subweb)
    
}


#################################################
fun.load.many.rdata <- function(file.list=NULL){
  # function to read and assign rdata to a list
  if(is.null(file.list)) stop("Must provide list of fully path files") 
  e1 = new.env()
  invisible(lapply(file.list, load, envir = e1))
  my_list = as.list(e1)
  
  return(my_list)
}
#################################################
fun.PRESENCE.ABSENCE <- function(species.dist=NULL, threshold=NULL, opt.only = NULL){
  if(is.null(species.dist)) stop("Must specify a species distribution dbf file")
  if(is.null(threshold)) stop("Must specify a cell habitat % threshold to consider a species present")
  if(is.null(opt.only)) stop("Must specify whether to use only optimal habitats, or both optimal and secondary habitats")
  
  if(is.null(species.dist$VALUE_1_prct)) warning(paste("Species", i,"has no secondary habitat in the study area"))
  if(is.null(species.dist$VALUE_2_prct)) warning(paste("Species", i,"has no primary habitat in the study area"))
  
  species.dist$Present <- 0
  species.dist[which(species.dist$VALUE_2_prct > threshold),"Present"] <- 1
  
  if(!opt.only){
    species.dist[which(species.dist$VALUE_1_prct > threshold),"Present"] <- 1   # if using secondary habitats, add these as presences
  }
  
  return(species.dist[,c(1,length(species.dist))])
}

#################################################
# Function to transform spp distribution database files into rasters
fun.dbf2raster <- function(SPPPA, mask.dir = NULL){
  # SPPPA must be in this format - first colmun with CELL ID and second column with the value (to plot)
  if(is.null(mask.dir)) stop("Must specify the mask file directory")
  library(raster)
  
  maskID <- read.dbf(list.files(path = mask.dir, full.names = TRUE, pattern = ".img.vat.dbf$"))
  maskk <- raster(x = list.files(path = mask.dir, full.names = TRUE, pattern = ".img$"))
  
  spp <- maskID
  spp$val <- NA
  spp$PageName <- as.character(spp$PageName)
  row.names(spp) <- spp$PageName
  
  SPPPA$PAGENAME <- as.character(SPPPA$PAGENAME)
  SPPPA[,2] <- as.numeric(as.character(SPPPA[,2]))
  row.names(SPPPA) <- SPPPA$PAGENAME
  
  cellID <- as.character(SPPPA$PAGENAME)
  if( nrow(spp[cellID,]) != nrow(SPPPA[cellID,])) stop("Cell IDs do not match")
  spp <- spp[cellID,] 
  spp$val <- SPPPA[,2]
  
  xx <- values(maskk)
  
  if( length(xx[!is.na(xx)]) != nrow(spp)) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val[xx[!is.na(xx)]]
  
  values(maskk) <- xx
  return(maskk)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

loadTxt <- function(filePath){
  #loads an text df file, and returns it
  props  <- t(read.table(file = filePath,header = TRUE))
  colnames(props)[1] <- 'Pagename' # First element is always Pagename
  props <- props[1,]            # As vector, keeping the names of elements
  get(ls()[ls() != "filePath"])
}

foo <- function (...) 
{
  dargs <- list(...)
  if (!all(vapply(dargs, is.vector, TRUE))) 
    stop("all inputs must be vectors")
  if (!all(vapply(dargs, function(x) !is.null(names(x)), TRUE))) 
    stop("all input vectors must be named.")
  all.names <- unique(names(unlist(dargs)))
  out <- do.call(rbind, lapply(dargs, `[`, all.names))
  colnames(out) <- all.names
  out
}

###### Network properties #######################
#################################################
# Function to calculate network prop

fun.net.prop <- function(web, dietcat=NULL){
  #   Metrics calculated:
  #   1: Chain length (Trophic level)
  #   2: Taxa  
  #   3: complexity
  #   4: strategy

  #NEEDS
  library(igraph)
  
  # Trophic level
  community <- Community(data.frame(node = colnames(web)), 
                         trophic.links = PredationMatrixToLinks(web),
                         properties = list(title = "Test2"))
  
  community <- RemoveCannibalisticLinks(community, title='community');
  TL <- PreyAveragedTrophicLevel(community)
  
  if(!is.null(dietcat)) webless <- web[!row.names(web) %in% dietcat,!colnames(web) %in% dietcat, drop=FALSE] # without diet category
  if(!is.null(dim(web))){
    
    Chain.res <- c(meanTL = mean(TL), maxTL= max(TL))
       
    # TAXA
    nCarnivorous_cat <- c("Invertebrates","Fish","DomesticAnimals","Carrion")   # diet categories that are not plant based. Used to defined which basal are herbivore basal and non-herbivore basal species.
    
    nCarnivorous_cat_in_web <- nCarnivorous_cat[nCarnivorous_cat %in% row.names(web)]
    Basalspp <- colnames(webless)[which(colSums(webless)==0)]
    propB = sum(colSums(webless)==0)/ nrow(webless)
    
    # The subset function keeps matrix class of the object, which is required for the colSums in the case of only 1 basal spp is present
    if(length(nCarnivorous_cat_in_web) !=0) propBNonH <- sum(colSums(subset(x = web[nCarnivorous_cat_in_web,,drop=FALSE],select =Basalspp))!=0)/nrow(webless)
      else propBNonH <- 0
    
    propBH <- abs(propBNonH-propB)
    propT = sum(rowSums(webless)==0 & colSums(webless)!=0 )/ nrow(webless)
    sum(colSums(web) == 0) / ncol(web)
    
    Taxa.res <- c(propOmn = OmnivoryCUS(web,dietcat = dietcat),
                  propB = propB,
                  propBNonH = propBNonH,
                  propBH = propBH,
                  propT = propT,
                  propI = 1- (propB + propT))
    
    # Complexity
    iweb <- graph.adjacency(web, mode="directed")
    Complex.res <- c(S = nrow(webless),
                     L = sum(webless),
                     C = sum(webless)/length(webless[,1])^2,
                     clust.coef = transitivity(iweb),
                     char.path.length = mean(shortest.paths(iweb)))
    
    #Strategy
    strat.res <- c(GEN =MeanGenerality(webless),
                   # GEN =MeanGenerality_norm(webless), # Used in previous revision.
                   GENSD = SDGenerality_norm(webless),
                   VUL = MeanVulnerability(webless),
                   # VUL = MeanVulnerability_norm(webless), # Used in previous revision.
                   VULSD = SDVulnerability_norm(webless),
                   MAXSIM = Maxsim(webless))
    #Compilation of results
    res <- c(Complex.res, Taxa.res, strat.res, Chain.res) 
  } 
  
  return(res)
}

# Centrality metrics
Centrality <- function(web, dietcat){
 
  community <- Community(data.frame(node = colnames(web)), 
                         trophic.links = PredationMatrixToLinks(web),
                         properties = list(title = "Test2"))
  
  # Get all spp trophic levels
  TL <- PreyAveragedTrophicLevel(community)
  
  if(!is.null(dietcat)) {
    TL <-  TL[!(names(TL) %in% (dietcat))] # to remove diet cat from trophic level
    web <- web[!row.names(web) %in% dietcat,!colnames(web) %in% dietcat, drop=FALSE] # remove diet category from centrality measures
  }
 
# Now I'll keep the diet categories
  iweb <- graph.adjacency(web, mode="directed")
  res <- list()
  res$iweb        <- iweb  
  res$TL          <- TL
  res$DegreeA     <- degree(graph = iweb,mode = "total")
  res$DegreeIN    <- degree(graph = iweb,mode = "in")
  res$DegreeOUT   <- degree(graph = iweb,mode = "out")
  res$Unconnected <- names(which(res$DegreeA == 0)) # registered nodes that after diet catogery removal are unconnected (centrality = 0 or null)
  res$BC          <- betweenness(graph = iweb,directed = TRUE)
  res$Eig.Cent    <- eigen_centrality(graph = iweb,directed = FALSE,scale = TRUE)
  res$CC          <- closeness(iweb, vids = V(iweb), mode = c("in"), weights = NULL, normalized = FALSE)

  return(res)
}

InDegree <- TrophicGenerality <- NumberOfResources <- function(M){
  return(colSums(M));
}

OutDegree <- TrophicVulnerability <- NumberOfCosumers <- function(M){
  return(rowSums(M));
}

Degree <- function(M){
  return(InDegree(M)+OutDegree(M));
}

NormalisedGenerality <- function(M){
  return(TrophicGenerality(M)/(sum(M)/dim(M)[1]));
}

NormalisedVulnerability <- function(M){
  return(TrophicVulnerability(M)/(sum(M)/dim(M)[1]));
}

# Function used by fun.net.prop   
# food web metrics  ----
# Calculate if a spp is an omnivore based on the number of different trophic levels predated by a spp
IsOmnivoreCUS <- function(M, level = PreyAveragedTrophicLevel){
  #Use PredationMatrixToLinks() to create a Cheddar community from a predation
  community <- Community(data.frame(node = colnames(M)), trophic.links = PredationMatrixToLinks(M),
                         properties = list(title = "Test2"))
  community <- RemoveCannibalisticLinks(community, title='community');
  # get resource spp for each predator
  resource.spp <- ResourcesByNode(community)
  #get no. resources
  n.resources <- sapply(resource.spp, length)
  # get all spp trophic levels
  tl <- level(community)
  # get the number of different trophic levels predated upon
  resource.tl <- sapply(resource.spp, FUN = function(x){
    return(length(unique(tl[x])))
  })
  # a sp is an omnivore if it predates 2 or more spp of different trophic levels
  return(n.resources >= 2 & resource.tl >= 2)
}

OmnivoryCUS <- function(M, dietcat = NULL, level = PreyAveragedTrophicLevel){
  omnivs <- IsOmnivoreCUS(M, level = level)
  if(!is.null(dietcat)){
    omnivs <- omnivs[!names(omnivs) %in% dietcat]
  }
  return(sum(omnivs) / length(omnivs)) 
  
}

MeanFoodChainLength <- function(M){
  # Use PredationMatrixToLinks() to create a Cheddar community from a predation
  # matrix
  node <- 1:dim(M)[1];
  for(n in 1:length(node)){
    node[n] <- paste(node[n],'-');
  }
  
  pm <- matrix(M, ncol=dim(M)[2], dimnames=list(node, node), byrow=TRUE);
  community <- Community(nodes=data.frame(node=node), trophic.links=PredationMatrixToLinks(pm), properties=list(title='Community'));
  community <- RemoveCannibalisticLinks(community, title='community');
  chain.stats <- TrophicChainsStats(community)
  ch_lens <- (chain.stats$chain.lengths + 1)
  
  return(sum(ch_lens)/length(ch_lens));
}


MeanGenerality <- function(M){
  return(sum(colSums(M))/sum((colSums(M)!=0)));
}

MeanVulnerability <- function(M){
  return(sum(rowSums(M))/sum((rowSums(M)!=0)));
}

MeanGenerality_norm <- function(M){
  norm_g <- NormalisedGenerality(M)
  return(mean(norm_g[norm_g!=0]))
}

SDGenerality_norm <- function(M){
  norm_g <- NormalisedGenerality(M)
  return(sd(norm_g[norm_g!=0]))
}

SDGenerality <- function(M){
  return(sd(InDegree(M)[InDegree(M)!=0]));
}

MeanVulnerability_norm <- function(M){
  norm_v <- NormalisedVulnerability(M)
  return(mean(norm_v[norm_v!=0]))
}

SDVulnerability_norm <- function(M){
  norm_v <- NormalisedVulnerability(M)
  return(sd(norm_v[norm_v!=0]))
}

SDVulnerability <- function(M){
  return(sd(OutDegree(M)[OutDegree(M)!=0]));
}

# Maximum trophic similiarity (Function from https://github.com/opetchey/dumping_ground/blob/master/random_cascade_niche/FoodWebFunctions.r)
Maxsim <- function(web){
  sims <- matrix(0, length(web[,1]), length(web[,1]))
  for(i in 1:length(web[,1]))
    for(j in 1:length(web[,1]))
      sims[i,j] <- T.sim.ij(web, i, j)
    diag(sims) <- NA
    mean(apply(sims, 1, function(x) max(x[!is.na(x)])))
}

T.sim.ij <- function(web, i, j){
  same <- sum(web[i,] & web[j,]) + sum(web[,i] & web[,j])
  total <- sum(web[i,] | web[j,]) + sum(web[,i] | web[,j])
  same / total
}
