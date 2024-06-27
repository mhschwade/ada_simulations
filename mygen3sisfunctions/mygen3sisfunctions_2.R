# Script with utilitary functions related to usage of gen3sis package.
################################################################################
# Author: Marcelo H. Schwade
# Last modification: Apr 4, 2024.
# 
################################################################################

# Basic packages.
require(gen3sis)
require(ape)
require(tidyr)
source("phylo_gen3sis_tools_1.R")



################################################################################
# Built Community Matrix at specific time.
################################################################################
################################################################################

built_community_matrix <- function(sim, times=c(0), simulation_path="./", output_path="./"){
  # Set the path with the simulation results.
  setwd(file.path(simulation_path, sim))

  # Simulation resume.
  gen3sis_simulation <- readRDS("sgen3sis.rds")
  
  # Basic parameters of the simulation.
  Ntimesteps <- length(gen3sis_simulation$summary$phylo_summary[,"total"]) - 2
  Nsites <- length(gen3sis_simulation$summary$`richness-final`[,1])
  
  Com_Matrix_list <- list()
  
  for(t in times){
    # Load the gen3sis simulation data.
    nSpeciesTotal <- gen3sis_simulation$summary$phylo_summary[(Ntimesteps - t + 2), "total"]
    nSpeciesAlive <- gen3sis_simulation$summary$phylo_summary[(Ntimesteps - t + 2), "alive"]
    
    # Species data at present time.
    metacommunity_at_t <- readRDS(file.path("species", paste("species_t_", t, ".rds", sep="")))

    # This suits only to verify the richness compatibility.
    ##############################################################################
    
    species_presence <- list()
    
    for(species_id in 1:nSpeciesTotal){
      species_id_presence <- list(id=species_id, 
                                  sites_with_presence=as.numeric(names(metacommunity_at_t[[species_id]]$abundance)),
                                  abundance=metacommunity_at_t[[species_id]]$abundance)
      species_presence <- c(species_presence, list(species_id_presence))
      names(species_presence)[length(species_presence)] <- as.character(species_id)
    }
    
    # Verify the names of the list (names or the ids of species)
    names(species_presence)

    # Create the community matrix.
    community_matrix <- matrix(0, nrow=Nsites, ncol=nSpeciesTotal)
    colnames(community_matrix) <- names(species_presence)
    
    for(species_id in names(species_presence)){
      community_matrix[species_presence[[species_id]]$sites_with_presence, species_id] <- species_presence[[species_id]]$abundance
    }

    Com_Matrix_list <- c(Com_Matrix_list, list(community_matrix))    
    
    # Verify if the richness in each site matches with the final richness of the simulation.
    # richness_final <- replace_na(gen3sis_simulation$summary$`richness-final`[,3], 0)
    # richness <- rowSums(ifelse(community_matrix>0, yes=1, no=0))
    # if(all(richness==richness_final)){
    #   print("Richness verification is OK!")
    # }
    # else{
    #   print("An incompatibility was found in the richness verification. The process was aborted.")
    #   break
    # }
    
  }
  
  names(Com_Matrix_list) <- as.character(times)
  
  return(Com_Matrix_list)
}


sim_path <- "~/Desktop/simple_simulation_20231114/simulation_outputs5/default_config"
sim <- "202311301511-00304"
Coms <- built_community_matrix(sim, times=c(0, 10), simulation_path=sim_path, output_path="~/Desktop")
View(Coms[["10"]])
require(vegan)
install.packages("arkhe")
require(arkhe)
dim(Coms[[1]])
com0<-arkhe::remove_zero(Coms[[1]], margin=2, all=T)
com1<-arkhe::remove_zero(Coms[[2]], margin=2, all=T)
class(com0)
dim(com0)
dim(com1)
View(com0)
cor(Coms[[1]], Coms[[2]])



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


