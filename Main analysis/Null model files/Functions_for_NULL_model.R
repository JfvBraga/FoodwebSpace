
# Plot food web by TL
#################################################
plot.graph.trophic.level <- function(web,...){
  
  TL <- TrophInd(as.matrix(web))$TL
  
  iweb <- graph.adjacency(web, mode="directed")  
  
  lay<-matrix(nrow=length(V(iweb)),ncol=2) # create a matrix with one column as runif, the other as trophic level
  lay[,1]<- runif(n = length(V(iweb)))
  lay[,2]<- TL
  
  plot.igraph(iweb,layout=lay)
}

################################################


#################################################
# Function to build local networks
sub_web <- function(metaweb = NULL, SPPCODE = NULL,
                    
                    SPP.HAB = NULL, PIX.HAB = NULL,
                    
                    HELP=TRUE, dietcat){
  
  # Function to subset the meta-web square matrix according
  
  # a subgroup of species into a sub web square matrix 
  
  if(HELP) warning("Description of HAB arguments: if spp x hab and pixel x hab matrices are provided habitat overlap will be used to filter spp links.
                   
                   \n    SPP.HAB: spp x habitat matrix
                   
                   \n    PIX.HAB: pix x habitat matrix
                   
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
    
    subweb <- metaweb[SPPCODE, SPPCODE]
    
    
    
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

# Functions to sample an exact number of species based on empirical distributions
#################################################
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

sampling_spp <- function(y,
                         S_richness, 
                         metaweb, 
                         prop.taxa, 
                         sampled_line_ref = NULL){
  # browser()
  if(is.null(sampled_line_ref)) {
    lines.to.sample <- which(S_richness == y)
    sampled_line_ref <- prop.taxa[resamp(lines.to.sample,size = 1),]
  }
  
  
  spp_pools <- colnames(BARMdiet.binary)[!(colnames(BARMdiet.binary) %in% dietcat)]  # Species within metaweb, without dietcat
  
  Asp <- if(sampled_line_ref$A != 0)
    sample(x = grep(pattern = "A",spp_pools,value = TRUE), round(sampled_line_ref$A),replace = F) else
      NULL
  
  Bsp <- if(sampled_line_ref$B != 0 )
    sample(x = grep(pattern = "B",spp_pools,value = TRUE),round(sampled_line_ref$B),replace = F) else
      NULL
  
  Msp <-  if(sampled_line_ref$M != 0 )
    sample(x = grep(pattern = "M",spp_pools,value = TRUE),round(sampled_line_ref$M),replace = F) else 
      NULL
  
  Rsp <- if(sampled_line_ref$R != 0 )
    sample(x = grep(pattern = "R",spp_pools,value = TRUE),round(sampled_line_ref$R),replace = F) else 
      NULL
  
  spp_web <- list() ; spp_web$species <- c(Asp,Bsp,Msp,Rsp)
  spp_web$sampled_line_ref <- sampled_line_ref
  return(spp_web)
}

#################################################


#################################################
# Null model function
random_network <- function(y,
                           metaweb = NULL,
                           S_richness = NULL,
                           prop.taxa = NULL,
                           ...) {
  # Arguments:
  # y: a species richness value being randomized
  # metaweb: adjacency matrix of the metaweb
  # S_richness: vector with all species richness across pixel
  # prop.taxa: matrix site x taxon. group containing the amount of species for each group at each pixel
  
  # browser()
  
  # First run 
  # ------------------------------
  # Sampling the species 
  spp_web <- sampling_spp(y,
                          S_richness,
                          metaweb = metaweb,
                          prop.taxa)  
  
  
  spp_web_diet <- c(dietcat, spp_web$species)                                                                            # vector with sampled species and diet categories
  
  # If spp are found in that pixel, then add the diet categories to the analysis
  web <- sub_web(...,
                 metaweb = metaweb,
                 SPPCODE = spp_web_diet,
                 SPP.HAB = NULL,
                 PIX.HAB = NULL,
                 HELP = FALSE)             # Select those species from the metaweb to a new network + Diet categories

  if(!is.matrix(web)) web <- as.matrix(web)             # Ensure output is a matrix
  
  # Evaluation process
  # ------------------------------  
  spp <- colnames(web)[!(colnames(web) %in% dietcat)]   # vector with number of nodes that are not diet catgories. This value may differ from y, since sub_web function removes unconnected nodes
  
  # Assessing if the spp vector contains the same proportion of taxa as the target
  same_tg <- length(spp_web$sampled_line_ref[spp_web$sampled_line_ref!= 0]) == length(table(substr(spp, 1, 1))) 

  dif_tg_sr <- if(same_tg) {sum(abs(spp_web$sampled_line_ref[spp_web$sampled_line_ref!= 0] - table(substr(spp, 1, 1)))) != 0 } else TRUE   
  # if dif_tg_sr is TRUE, then it mean that the random web does not have the same amount of taxa as the target. So, re-run

  # Second or more runs 
  # ------------------------------
  while(dif_tg_sr){
    spp_web <- sampling_spp(y,
                            S_richness,
                            metaweb = metaweb,
                            prop.taxa) 
    spp_web_diet <- c(dietcat, spp_web$species)
    
    web <- sub_web(...,
                   SPPCODE = spp_web_diet,
                   metaweb = metaweb,
                   PIX.HAB = NULL,
                   SPP.HAB = NULL,
                   HELP = FALSE)             # Select those species from the network to a new network + Diet categories
    
    if(!is.matrix(web)) web <- as.matrix(web)
    
    # Evaluation process
    # ------------------------------  
    spp <- colnames(web)[!(colnames(web) %in% dietcat)]  
    
  
    # Assessing if the spp vector contains the same proportion of taxa as the target
    same_tg <- length(spp_web$sampled_line_ref[spp_web$sampled_line_ref!= 0]) == length(table(substr(spp, 1, 1))) 
    
    dif_tg_sr <- if(same_tg) {sum(abs(spp_web$sampled_line_ref[spp_web$sampled_line_ref!= 0] - table(substr(spp, 1, 1)))) != 0 } else TRUE   
    # if dif_tg_sr is TRUE, then it means that the random web does not have the same amount of taxa as the target. So, re-run
  }
  
 
  return(web)
}

#################################################

# Functions to read and compile results
#################################################
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

compile_results <- function(files.loc.vector = fls,y){
  flsnames <- gsub(pattern = paste0(results.path,y,"//"),replacement = "",x = files.loc.vector)
  flsnames <- gsub(pattern = ".txt$",replacement = "",x = flsnames)
  
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

##################################################