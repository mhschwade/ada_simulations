# Script with utilitary functions related to usage of gen3sis package.
################################################################################
# Author: Marcelo H. Schwade
# Last modification: Apr 4, 2024.
# 
################################################################################

# Basic packages.
require(gen3sis)

################################################################################
# Function to resume status of the simulations.
################################################################################
resume_simulations <- function(simulation_main_dir, simulations_list){
  
  simulation_infos <- data.frame(simulation_name=simulations_list)
  nsim <- length(simulations_list)
  {
    simulation_infos$status <- rep(NA, nsim)
    simulation_infos$Nsites <- rep(NA, nsim)
    simulation_infos$Ntimesteps <- rep(NA, nsim)
    simulation_infos$trait_variability <- rep(NA, nsim)
    simulation_infos$dispersal_rate <- rep(NA, nsim)
    simulation_infos$delta <- rep(NA, nsim)
    simulation_infos$niche_breath <- rep(NA, nsim)
    simulation_infos$random_seed <- rep(NA, nsim)
    simulation_infos$nSpeciesAlive <- rep(NA, nsim)
    simulation_infos$nSpeciesTotal <- rep(NA, nsim)
    simulation_infos$maxFinalRichness <- rep(NA, nsim)
    simulation_infos$matrixOcupanceFinal <- rep(NA, nsim)
  }
  
  for(i in 1:nsim){
    
    # Simulation name.
    sim_path <- simulations_list[i]
    
    # Loop status.
    print(paste(paste("Get data from sim ", sim_path), ".", sep=""))
    
    simulation_file <- file.path(simulation_main_dir, "default_config", sim_path, "sgen3sis.rds")
    
    if(file.exists(simulation_file)){
      
      sim_result <- readRDS(simulation_file)
      
      simulation_infos$status[i] <- "OK"
      
      # Get the parameters of the simulation.
      Ntimes <- length(sim_result$summary$phylo_summary[, "total"])
      simulation_infos$Nsites[i] <- nrow(sim_result$summary$`richness-final`)
      simulation_infos$Ntimesteps[i] <- Ntimes - 1
      simulation_infos$trait_variability[i] <- sim_result$parameters$user$trait_variability
      simulation_infos$dispersal_rate[i] <- sim_result$parameters$user$dispersal_rate
      simulation_infos$delta[i] <- sim_result$parameters$user$delta
      simulation_infos$niche_breath[i] <- sqrt(sim_result$parameters$user$niche_var2/2)
      simulation_infos$random_seed[i] <- sim_result$parameters$gen3sis$general$random_seed
      simulation_infos$nSpeciesAlive[i] <- sim_result$summary$phylo_summary[Ntimes, "alive"]
      simulation_infos$nSpeciesTotal[i] <- sim_result$summary$phylo_summary[Ntimes, "total"]
      simulation_infos$maxFinalRichness[i] <- max(sim_result$summary$`richness-final`[,3], na.rm=T)
      simulation_infos$matrixOcupanceFinal[i] <- round(100*sum(ifelse(na.omit(sim_result$summary$`richness-final`[,3])>0, yes=1, no=0))/length(na.omit(sim_result$summary$`richness-final`[,3])), digits=2)
      
    }
    else{
      simulation_infos$status[i] <- "FAIL"
    }
  }

  write.table(simulation_infos, file=file.path(simulation_main_dir, "simulation_infos.txt"), sep="\t")
  
  return(simulation_infos)
}


