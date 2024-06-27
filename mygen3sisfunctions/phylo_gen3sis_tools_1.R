# Required packages.
require(ape)
require(phytools)
require(phangorn)

################################################################################
# getAncestorPath
################################################################################
# Get the sequence of ancestors of a node or tip of a phylogeny.
# `tree`: a phylogenetic tree in phylo format.
# `node`: a node/tip label or node/tip index.
# return: a vector with a sequence of heights of ancestor nodes of the focal 
# node. The names of elements of this vector corresponds to the indexes of the
# ancestor nodes.
getAncestorPath <- function(tree, node){
  
  # Verify if the focal node informed is in an appropriate format.
  if(is.character(node)){
    node_index <- which(tree$tip.label==node)
  }
  else if(is.numeric(node)){
    node_index <- node
  }
  else{
    print("Error: node format unrecognized. Please, use node label or node index.")
  }
  
  # Get the indexes of the ancestor nodes of the focal node.
  anc <- nodepath(tree, from=node_index, to=(tree$Nnode+2))
  
  # Create a vector to get the heights of the ancestor nodes.
  node_heights <- rep(NA, length(anc))
  # Set the indexes of ancestor nodes as names of the vector elements.
  names(node_heights) <- as.character(anc)
  
  # Get the heights of the ancestor nodes and sum it with the height of the 
  # root.
  for(i in 1:length(anc)){
    node_heights[i] <- tree$root.edge + phytools::nodeheight(tree, node=anc[i])
  }
  
  # Return the vector with the indexes of ancestor nodes an its respective 
  # heights.
  return(node_heights)
}

################################################################################
# mapAncestorSpecies
################################################################################
# This function map the ancestor species index in their correspondent species 
# names at that time.
# This function work based in the fact that species names are assigned in the 
# same order that species emerge.
# `tree`: phylogenetic tree in phylo format.
# `ancestor_indexes`: a vector with a sequence of ancestor indexes.
# return: a vector with the ids (gen3sis species ids) of the ancestors.
mapAncestorSpecies <- function(tree, ancestor_indexes){

  # Number of ancestor nodes.
  n_nodes <- length(ancestor_indexes)
  # Coerce the ancestor indexes to a numeric format.
  ancestor_indexes <- as.numeric(ancestor_indexes)
  
  # Create a vector to get the ids (gen3sis species ids) associated to each 
  # ancestor.
  ancestor_ids <- rep(NA, n_nodes)

  # Get the id for each ancestor of the sequence.
  for(i in 1:n_nodes){
    # Index (in the phylogenetic tree) of the ancestor (node) i.
    node <- ancestor_indexes[i]
    # Get the list of indexes of the descendants (tips) of the node i.
    descendant_indexes <- phangorn::Descendants(tree, node)[[1]]
    # Get the labels (tip names) of the descendants of node i.
    descendant_names <- tree$tip.label[descendant_indexes]
    # Get the gen3sis ids of the descendants from the tip names.
    descendant_ids <- as.numeric(gsub("species", "", descendant_names))
    # Use the pattern of nomenclature of species in the gen3sis package to get 
    # the id of the ancestor node.
    # The lesser id among the descendants corresponds to the id of the founder 
    # of the clade, i. e. the id of the ancestor node.
    ancestor_ids[i] <- min(descendant_ids)
  }
  
  return(ancestor_ids)
}

################################################################################
# listAncestorSpecies
################################################################################  
# Function to create a temporal list of ancestor species from a species ancestor
# path (a data.frame within species ancestor id and the time step in which it
# speciate.).
# This list can be used to get information of a lineage along the phylogeny.
# This data can be used in plots, or to calculate metrics for this lineage, 
# beyond tracking the evolution of its parameters.
# `species_ancestor_path`: a data.frame with the heights of the nodes of the 
# focal lineage and the ids of the nodes (ancestors). 
# `Ntimesteps`: number of time-steps of the gen3sis simulation.
# return: a vector with a temporal series of the gen3sis id of the ancestor of 
# focal species at each time.
listAncestorSpecies <- function(species_ancestor_path, Ntimesteps){
  
  # Create a vector to get the id of the ancestors of focal species at each 
  # time-step. 
  ancestor_species_id <- rep(NA, Ntimesteps)
  
  # Number of different ids of ancestors.
  k <- length(species_ancestor_path$ancestor_height) + 1
  # Create a vector with the times when occur speciations in the lineage.
  time_threshold <- c(Ntimesteps, species_ancestor_path$ancestor_height)
  # Create a sequence with the ids of the ancestor of the focal species.
  ancestor_ids <- c(species_ancestor_path$ancestor_id[1], species_ancestor_path$ancestor_id)
  
  # For each time-step get the id of the ancestor of the focal species.
  for(t in 1:Ntimesteps){
    
    # Set the ancestor id at time-step `t`.
    ancestor_species_id[t] <- ancestor_ids[k]
    
    # Verify if the time exceed the next `time_threshold`.
    # If TRUE (exceed), then change the id to set in the next time-steps.
    if(t>=time_threshold[k]){
      k <- k-1
    }
  }
  
  return(ancestor_species_id)
}


# Function to rename tips and nodes of the phylogeny.
# `tree`: a phylogenetic tree in phylo format.
renamePhylo <- function(tree){

  # Number of nodes + tips.
  nTipNodes <- 2*tree$Nnode + 1
  # Create a data.frame to get information of node heights.
  nameIndexGuide <- data.frame(id=rep(NA, nTipNodes), 
                           height=rep(NA, nTipNodes))
  
  # For each node/tip, it gets the id and the height of the node and use these 
  # information to make a new name for the nodes.
  for(i in 1:nTipNodes){
    #i=1130
    #  
    species_ancestors <- getAncestorPath(tree, node=i)
    ancestor_ids <- mapAncestorSpecies(tree, names(species_ancestors))
    ancestor_ids <- rev(ancestor_ids)
    
    nAncestors <- length(ancestor_ids)
    new_name <- paste("sp", ancestor_ids[1], sep="_")
    id <- ancestor_ids[1]
    
    if(nAncestors>1){
      for(id in ancestor_ids[2:nAncestors]){
        new_name <- paste(new_name, id, sep=".")
      }
    }
    
    rownames(nameIndexGuide)[i] <- new_name
    nameIndexGuide$id[i] <- id
    nameIndexGuide$height[i] <- species_ancestors[1]

  }
  
  # Rename the tips and the nodes of the phylogeny.
  tree$tip.label <- rownames(nameIndexGuide)[1:(tree$Nnode + 1)]
  tree$node.label <- rownames(nameIndexGuide)[(tree$Nnode + 2):nTipNodes]
  
  new_tree <- list(tree=tree, nameIndexGuide=nameIndexGuide)
  
  return(new_tree)
}


################################################################################
# Function to plot species abundance.
plotSpeciesAbundance <- function(metacommunity,
                                 landscape_basis,
                                 species_id=1,
                                 colorscale=c(0, 1, 3, 10, 30, 100, 300),
                                 colors=c("#ffffff", "#ffff80", "#ffff00", "#ff8000", "#ff4000", "#ff0000", "#800000")){
  
  
  # Creating a a basic data.frame to species abundance data.
  species_data <- data.frame(
    "long"=landscape_basis$long, 
    "lat"=landscape_basis$lat,
    "abundance"=landscape_basis$hab)
  
  rownames(species_data) <- rownames(landscape_basis)
  
  # Site names of the sites where the species occur.
  select_site_names <- names(metacommunity[[species_id]]$abundance)
  species_data[select_site_names, "abundance"] <- 
    metacommunity[[species_id]]$abundance[select_site_names]
  
  # Defining color scale and breaks position.
  colorbreaks <- truncatedLog(colorscale)
  colorscale_labels <- as.character(colorscale[-length(colorscale)])
  
  # Plot.
  species_abundance_plot <- ggplot(data=species_data) +
    geom_raster(mapping=aes(x=long, y=lat, fill=truncatedLog(abundance))) +
    scale_fill_gradientn(
      colors=colors, 
      name='Abundance',
      space="Lab", # I don't know what this is for.
      limits=c(min(colorbreaks), max(colorbreaks)),
      values=scales::rescale(colorbreaks),
      breaks=colorbreaks[-length(colorbreaks)],
      labels=colorscale_labels) +
    coord_equal() +
    theme_minimal()
  
  return(species_abundance_plot)  
}

