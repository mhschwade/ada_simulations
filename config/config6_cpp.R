
######################################
###            METADATA            ###
######################################
# gen3sis configuration
#
# Version: 1.0
#
# Author: Marcelo Henrique Schwade
#
# Date: 04.08.2024 (MM.DD.YYYY)
#
# Landscape: landscape
#
# Publications: these data will are used in a publication in colaboration with 
# Leandro Duarte, Gabriel Nakamura, Renan Maestri, Fabricio Villalobos, Carolina 
# Prauchner and Juliene costa.
#
# Description: this is a first simulation to generate data to test models of
# species ancestral range reconstruction.
# 
#
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation.
random_seed = 42

# set the starting time step or leave NA to use the earliest/highest time-step.
start_time = 8000

# set the end time step or leave as NA to use the latest/lowest time-step (0).
end_time = 0

# maximum total number of species in the simulation before it is aborted.
max_number_of_species = 100000

# maximum number of species within one cell before the simulation is aborted.
max_number_of_coexisting_species = 10000

# a list of traits to include with each species
# a "dispersal" trait is implicitly added in any case
trait_names = c("cond1", "cond2")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list("cond1"=NA, "cond2"=NA) 

# Parameters to change.
dispersal_rate <- 0.5
#niche_breath <- 0.05 # sd=Proportion of the gradient.
niche_var2 <- 2*(0.05^2) #niche_breath
trait_variability <- 0.01 # Proportion of the gradient.
# Delta - Rate of mixing of the traits within the cluster.
delta <- 1 # Populations totally mixed into each cluster (original proposal of Hagen et al. 2021)
# Initial distribution of the ancestor species.
range_initial <- c(-5, 0, -5, 0) # longmin, longmax, latmin, latmax.

# Control parameters of timestep_observer.
observation_interval <- 10

# Import/compile Rcpp functions (to ecology dynamics).
# Function: lotka_volterra_comp()
# Function: lotka_volterra_comp_neutral()
Rcpp::sourceCpp("../mygen3sisfunctions/lotka_volterra_comp_neutral2.cpp")
lotka_volterra_comp_neutral2 <- lotka_volterra_comp_neutral2

######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation and collect additional information if desired.
end_of_timestep_observer = function(data, vars, config){
  # the list of all species can be found in data$all_species
  # the current landscape can be found in data$landscape
  
  # Verify speciation.
  speciation_event <- vars$n_new_sp_ti > 0

  ##############################################################################
  # Save the species data at each 10 timesteps and when speciations happen.
  if((vars$ti%%config$user$observation_interval==0) | speciation_event){
    
    # Save data of the species.
    config <- dynGet("config")
    data <- dynGet("data")
    vars <-  dynGet("vars")
    
    # Save abundances of the species.
    dir.create(file.path(config$directories$output, "abundance"), showWarnings=FALSE, recursive=TRUE)
    abundance <- lapply(data$all_species, function(x){return(x[["abundance"]])})
    names(abundance) <- sapply(data$all_species, function(x){x$id})
    saveRDS(object = abundance,
            file = file.path(config$directories$output, "abundance", paste0("abundance", "_t_", vars$ti, ".rds")))

    # Save traits of the species.
    dir.create(file.path(config$directories$output, "traits"), showWarnings=FALSE, recursive=TRUE)
    traits <- lapply(data$all_species, function(x){return(x[["traits"]])})
    names(traits) <- sapply(data$all_species, function(x){x$id})
    saveRDS(object = traits,
            file = file.path(config$directories$output, "traits", paste0("traits", "_t_", vars$ti, ".rds")))
    

    # Plot data at each time step.
    par(mfrow=c(2,1)) # To plot a mosaic of two graphs.
    
    # Plot richness.
    plot_richness(data$all_species, data$landscape)
    # Plot species 1 distribution.
    plot_species_presence(data$all_species[[1]], data$landscape)
    
  }
  
  mtext("STATUS", 1)
}


######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when 
# colonizing a cell during the dispersal.
initial_abundance = 1

# place species in the landscape:
create_ancestor_species <- function(landscape, config) {
  # stop("create the initial species here")
 
  # Define the distribution range of the ancestral species.
  # longmin, longmax, latmin, latmax.
  range <- config$user$range_initial
  # Using an initial range fully distributed.
  # The speciation starts only after the first ancestor colonizing all the sites.
  # Then, it's more practical and simple to start with the ancestor fully dispersed.
  
  # Select the sites within distribution range.
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] & co[, "x"] <= range[2] & 
    co[, "y"] >= range[3] & co[, "y"] <= range[4]
  
  initial_cells <- rownames(co)[selection]
  
  # Create the initial species in this sites.
  new_species <- create_species(initial_cells, config)
  
  # Set the species traits.
  new_species$traits[, "cond1"] <- landscape$environment[initial_cells, "cond1"]
  new_species$traits[, "cond2"] <- landscape$environment[initial_cells, "cond2"]
  # Here, we set the initial trait of each population species as the 
  # environmental condition on that local.
  
  # Return the list with species information.
  return(list(new_species))
  
}


######################################
###             Dispersal          ###
######################################

# the maximum range to consider when calculating the distances from local distance inputs.
max_dispersal <- 6.907755*dispersal_rate # Cutting radius.
# Distance where the accumulated probability of dispersal beyond is lesser than 0.001.
# This is computed to an exponential dispersal kernel.

# returns n dispersal values.
get_dispersal_values <- function(n, species, landscape, config) {
  # Generate values of dispersal.
  # Probabilistic dispersal using a dispersal kernel exponential.
  values <- rexp(n, rate=(1/config$user$dispersal_rate)) # is this a proper rate?
  
  return(values)
}


######################################
###          Speciation            ###
######################################

# threshold for genetic distance after which a speciation event takes place.
divergence_threshold = 500

# Define the factor by which the divergence is increased between geographically isolated population.
# can also be a matrix between the different population clusters.
# `species`:
# `cluster_indices`:
# `landscape`:
# `config`:
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  # stop("calculate divergence factor here")
  # Return 1 for each time step at which is accumulated genetic incompatibilities (?).
  return(1)
}


######################################
###            Evolution           ###
######################################

# mutate the traits of a species and return the new traits matrix.
apply_evolution <- function(species, cluster_indices, landscape, config) {
  # stop("mutate species traits here")
  
  # Define the variability of trait in a population in a time step.
  trait_variability <- config$user$trait_variability
  
  # Get the species traits.
  traits <- species[["traits"]]
  
  # Get the cells.
  cells <- rownames(traits)
  
  ##############################################################################
  # Homogenize the traits within a population cluster based on abundance.
  # So, it assume a fully mixed cluster and the absence of limitations to gene 
  # flow within the cluster, while between clusters the gene flow is zero.
  
  # Iterate about all the cluster
  for(cluster_index in unique(cluster_indices)){
    
    # Get the cells belonging to the cluster.
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    
    # Calculate the total abundance of the cluster.
    total_abundance <- sum(species$abundance[cells_cluster])
    
    # Calculate the abundance weighted average to the "cond1" and "cond2" traits.
    cluster_mean_cond1 <- sum(traits[cells_cluster, "cond1"]*species$abundance[cells_cluster])/total_abundance
    cluster_mean_cond2 <- sum(traits[cells_cluster, "cond2"]*species$abundance[cells_cluster])/total_abundance
    
    traits[cells_cluster, "cond1"] <- config$user$delta*cluster_mean_cond1 + (1 - config$user$delta)*traits[cells_cluster, "cond1"]
    traits[cells_cluster, "cond2"] <- config$user$delta*cluster_mean_cond2 + (1 - config$user$delta)*traits[cells_cluster, "cond2"]
  }
  
  # Generate a mutation/drift variation in the trait.
  # The variation is gaussian with mean 0 and standard deviation equal to `trait_variability`.
  traits[, "cond1"] <- traits[, "cond1"] + rnorm(length(traits[, "cond1"]), mean=0, sd=trait_variability)
  traits[, "cond2"] <- traits[, "cond2"] + rnorm(length(traits[, "cond2"]), mean=0, sd=trait_variability)
  
  # Return the modified traits.
  return(traits)
}


######################################
###             Ecology            ###
######################################

##############################################################################
# Define the parameters of Lotka-Volterra competition model (simplified).

bmax <- 3 # Maximum birth intrinsic rate.
# rmax is the maximum fundamental increment rate.
# rmax = bmax - 1 = 5 - 1 = 4 (all individuals dead after a generation time, 
# time-step)

Kr <- 500 # Species-local carrying capacity per recruitment. # Modification from config2_cpp.
# The maximum carrying capacity Kmax is
# Kmax = rmax*Kr
# Here, Kmax = (bmax-1)*Kr =(3-1)*1000 = 2*1000 = 2000

alpha <- 1 # Interespecific competition rate,
# Alpha equal to 1 represent an identical intra and inter specific competition.

# Environmental tolerance of species.
sigma2_cond1 <- config$user$niche_var2 # 2*sigma^2
sigma2_cond2 <- config$user$niche_var2 # 2*sigma^2

# Define abundance threshold.
abundance_threshold <- 1 # When abundance is lesser than 1 the population is extinct.

# Put all the parameters in a vector - please, the order should be respected.
lvc_pars <- c(bmax, Kr, alpha)


# called for every cell with all occurring species, this function calculates abundances and/or 
# who survives for each sites.
# returns a vector of abundances.
# set the abundance to 0 for every species supposed to die.
# `abundance`: 
# `traits`:
# `environment`:
# `config`:
apply_ecology <- function(abundance, traits, environment, config) {
  # stop("calculate species abundances and deaths here")
  
  # In this model, we consider the species from a neutral point of view in terms of 
  # theirs biological interactions. There are not difference in carrying 
  # capacities, basic increment rates or competitive ability between species.
  # Yet individuals from different species compete in the same way that 
  # individuals from different species. 
  # The only difference between species is its affinity with the environment in
  # each site, which result in an effective increment rate less than the basic 
  # rate.
  
  ##############################################################################
  # Extinguish the species with abundance below the threshold.
  
  # Extinguish the populations with abundance below the threshold.
  abundance[abundance < config$user$abundance_threshold] <- 0
  
  # Apply a Lotka-Volterra Competition Model (LVC Model, modified to include niche).
  # Here, in this config file (config6_cpp), we use a neutral version of LVC model. So, 
  # the environment and the environmental tolerances (sigmas) aren't passed as 
  # argument.
  # The function to integrate LVC model was coding in C++ using Rcpp.
  abundance <- config$user$lotka_volterra_comp_neutral2(abundance, config$user$lvc_pars, dt=0.01, times=1)
  
  # Apply drift using Poisson distribution with mean equal to deterministic abundance.
  # abundance <- rpois(n=length(abundance), lambda=abundance)
  
  # Apply abundance threshold to the updated populations.
  abundance[abundance < config$user$abundance_threshold] <- 0
  
  # Return the final abundance.
  return(abundance)
}


