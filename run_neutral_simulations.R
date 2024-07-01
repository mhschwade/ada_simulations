# Run gen3sis simulation.
# File name: 

# Author: 
# Marcelo Henrique Schwade

# Date: 
# June 20, 2024.

# Simulation path: 
# /home/mhs/Desktop/ada_simulations

# Output path:
# /home/mhs/Desktop/ada_simulations/neutral_simulations

# Config file: 
# config/config6_cpp.R

# Landscape:
# landscapes_20231114/landscapes/simple_landscapes_20231114.rds

# Description:
# Simulation to generate data to test species ancestral range reconstruction.
# This simulations use a neutral ecological model to the population dynamics of 
# the species.


################################################################################
# Import packages.
require(gen3sis)
require(raster)
require(parallel)
require(doParallel)
require(Rcpp)


################################################################################
# Get simulation path.
simulation_path <- "/Users/marceloschwade/Desktop/ada_simulations"

# Define the config file of the simulation.
config_file <- file.path(simulation_path, "config", "config6_cpp.R")

# Defining the basic paths to landscape and to save the results of the simulation.
landscape_path <- file.path(simulation_path, "landscapes_20231114", "landscapes", "distances")
output_directory <- file.path(simulation_path, "neutral_simulations")

# Number cores to use to parallelize the simulations.
no_cores <- 2

################################################################################
# Import config file and check its integrity.
config_object <- create_input_config(config_file)

# Before to start a simulation, we can verify if the config object works.
verify_config(config_object)
# If this verification returns TRUE it is good to go.

################################################################################
# Function to run simulations using a set of parameters.
# Function to run gen3sis changing the random seed, the initial range and other 
# basic parameters.
runGen3sisSim <- function(config_file, landscape_path, output_directory, par, save_state=NA, verbose=T){
  
  config <- create_input_config(config_file)
  
  config$gen3sis$general$random_seed <- par$rseed # Random seed of the simulation.
  config$user$trait_variability <- par$trait_variability # Trait variability.
  config$user$dispersal_rate <- par$dispersal_rate # Dispersal rate.
  config$user$delta <- 1 # Rate of mixture of the traits of the populations into a cluster.
  # For delta=1 the trait values of the populations are totally mixed into each cluster.
  config$user$niche_var2 <- 2*(par$niche_breath^2) # Environmental tolerance of the species. 
  config$user$range_initial <- par$range_initial # Initial range of distribution of the ancestor species.
  
  # Print the parameters of the simulation in the log file.
  #textlog <- paste(Sys.time(), par$niche_breath, par$dispersal_rate, par$trait_variability, par$rseed, par$range_initial[1], par$range_initial[2], par$range_initial[3], par$range_initial[4], sep='\t')
  #, config_file, landscape_path
  #time, niche_breath, dispersal_rate, trait_variability, random_seed, init_cond, config_file, landscape
  #write(textlog, file=file.path(output_directory, "simulations.log"), append=T)
  
  time_start <- Sys.time()
  simulation_results <- run_simulation(config=config, 
                                       landscape=landscape_path,
                                       output_directory=output_directory, save_state=save_state, verbose=verbose)
  time_end <- Sys.time()
  
  print(paste(paste("Simulation with seed =", pars$rseed), "finished."))
  print(time_end - time_start)
  return(simulation_results)
}

################################################################################
# Defining the parameters for the simulations.

# Random seeds, we are using 5 different random seeds.
rseeds <- c(344062, 953111, 329969, 575612, 364862)

# Initial range of the ancestor species, we are using 5 different initial ranges.
range_initial <- list(c(-5, 0, -5, 0), 
                      c(8, 13, -10, -5),
                      c(-3, 2, -17, -12),
                      c(-12, -7, -24, -19),
                      c(10, 15, -30, -25))
names(range_initial) <- c("A", "B", "C", "D", "E")

# Dispersal rate (dispersal), niche breath (ecology, it doesn't have effect in 
# neutral simulations) and trait variability (evolution).
dispersal_rates <- c(0.6)
niche_breaths <- c(0.05) # This parameter is not in use for neutral dynamics.
trait_variability <- 0.01

# The number of simulations is given by the multiplicity of the parameters.
# How only `rseeds` and `range_initial` have more than 1 'value'/state, then the 
# number of simulations is merely 5x5 = 25 simulations.

# Put all set of parameters into a list. Each set of parameters (element of the 
# list) is used to run a simulation.
pars <- list()
nsim <- 0
for(rseed in rseeds){
  for(range_i in 1:5){
    for(niche_breath in niche_breaths){
      for(dispersal_rate in dispersal_rates){
        pars <- c(pars, 
                  list(list(rseed=rseed,
                            dispersal_rate=dispersal_rate, 
                            niche_breath=niche_breath, 
                            trait_variability=trait_variability,
                            range_initial=range_initial[[range_i]])))
        nsim <- nsim + 1
      }
    }
  }
}


# Run the simulations!!!
################################################################################
# Finally, we can run a simulation with the landscape and the configuration that
# we define. 
####################################

# Run one simulation to test if the script is working.
sim.results <- runGen3sisSim(config_file, landscape_path, output_directory, pars[[1]], save_state=NA)

# Run all the missing simulations using parallelization.
# Create the cluster to run the simulations in parallel using no_cores.
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Run the simulations in parallel. Each simulation run in a core.
foreach(i=2:nsim, .packages=c('gen3sis', 'raster', 'Rcpp')) %dopar% {
  sim.results.par <- runGen3sisSim(config_file, landscape_path, output_directory, pars[[i]], verbose=F)
}

stopCluster(cl)

# Task finished.
cat("This script has finished its execution.\n")
cat(paste("The simulation results can be found into the path: ", output_directory, "\n"))
