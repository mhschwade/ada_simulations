################################################################################
# Create a simple landscape to simulate the eco-evolutionary dynamics of a 
# lineage.

################################################################################
# Defining work directory.
simulation_path <- "/home/mhs/Desktop/simple_simulation_20231114"
setwd(simulation_path)

################################################################################
# Import packages.
require(gen3sis)
require(raster)
require(parallel)
#require(plyr)

################################################################################
# Defining geographical landscape and coordinates.
################################################################################

print("Defining parameters and functions. 1/6")

# Fixing geographical limits.
latmin <- -30
latmax <- 0
longmin <- -15
longmax <- 15

# Defining the size of our landscape.
Ncol <- 30
Nrow <- 30
Nsites <- Ncol*Nrow
# N sites Ncol*Nrow = 100.

# Defining the temporal variability of the landscape.
Ntimesteps <- 10000

# Mapping site-center coordinates: lat , long.
coordinates <- data.frame(lat = rep(seq(latmax, latmin, length.out=Nrow), each=Ncol), long = rep(seq(longmin, longmax, length.out=Ncol), times=Nrow))
# It's needed ordering latitude from maximum to minimum for that the matrix/raster be properly built.

# Checking coordinates.
#coordinates$lat
#coordinates$long

# Set random seed.
set.seed(1859)

################################################################################
## Geographical landscape.
################################################################################

# Define habitability: plain landscape: all the sites are habitable.
# probh <- 0.9
# habitability <- rbinom(Nsites, size=1, prob=probh)
# Import the cost of landscape.
landscape_cost_matrix <- read.table(file="landscape_cost.txt", sep=" ")
landscape_cost <- data.frame(long=coordinates$long, 
                             lat=coordinates$lat, 
                             cost=rep(NA, Nsites))

for(i in 1:Nrow){
  for(j in 1:Ncol){
    landscape_cost$cost[(i - 1)*Nrow + j] <- landscape_cost_matrix[i, j]
  }
}

landscape_cost_raster <- rasterFromXYZ(landscape_cost)

habitability <- data.frame(long=coordinates$long, 
                           lat=coordinates$lat, 
                           hab=ifelse(landscape_cost$cost==1, yes=1, no=NA))


habitability_raster <- rasterFromXYZ(habitability)

# Viewing geographical matrix.
par(mfrow=c(1,2))
plot(landscape_cost_raster, main="landscape cost", xlab="long", ylab="lat", col=c("yellow", "black"), asp=3)
plot(habitability_raster, main="Habitability", xlab="long", ylab="lat", col=c("gray", "green"), asp=3)

################################################################################
# Environmental conditions.
################################################################################
# 3 environmental conditions.

####################################################################################
# Function to built the landscape (using sapply in parallel for time).
####################################################################################
create_landscape_environment_2Dpar <- function(habitability, landscape_cost, 
                                               dims, ncores=2, sd=0){
  
  # To counting time of processing.
  time_start <- proc.time()
  
  # Spatio-temporal dimensions.
  Nrow <- dims[1]
  Ncol <- dims[2]
  Ntimesteps <- dims[3]
  
  ################################################################################
  # Creating arrays of environmental variables. Each array has three dimensions.
  envir_cond <- list()
  
  # Creating environmental data.frames.
  envir_cond1_df <- as.data.frame(matrix(NA, nrow=(Nrow*Ncol), 
                                         ncol=(length(c(Ncol,Nrow))+Ntimesteps)))
  envir_cond2_df <- as.data.frame(matrix(NA, nrow=(Nrow*Ncol), 
                                         ncol=(length(c(Ncol,Nrow))+Ntimesteps)))
  
  #envir_cond3_df <- as.data.frame(matrix(NA, nrow=(Nrow*Ncol), 
  #                                       ncol=(length(c(Ncol,Nrow))+Ntimesteps)))
  
  ##############################################################################
  # Functions to determine the environmental conditions on space.
  
  # Environmental condition 1. Range: ~ 0 : 30. Latitudinal.
  evaluate_cond1 <- function(site_habitability, sd=0){
    
    if(!(is.na(site_habitability["hab"]))){
      cond1 <- 100 - 3.33*abs(site_habitability["lat"]) + sd*rnorm(1)
        #80 - 1.2*abs(site_habitability["lat"]) + sd*rnorm(1) #envir_mean
    }
    else{
      cond1 <- NA
    }
    
    return(cond1)
  }
  
  # Environmental condition 2. Range: ~ 10 : 90. Centered.
  evaluate_cond2 <- function(site_habitability, sd=0){
    
    if(!(is.na(site_habitability["hab"]))){
      cond2 <- 3.33*(site_habitability["long"] + 15) + sd*rnorm(1)
        #20 + 1.2*(site_habitability["long"] + 25) + sd*rnorm(1)
    }
    else{
      cond2 <- NA
    }
    
    return(cond2)
  }
  
  evaluate_landscape <- function(timestep_name, habitability, evaluate_cond, sd=0){
    
    # Loops for the spatial dimensions.
    envir_cond <- apply(habitability, 1, evaluate_cond, sd)
    
    # Time-step name.
    names(envir_cond) <- timestep_name
    
    return(envir_cond)
  }
  
  # Naming the timesteps.
  timesteps_names <- paste("X", 1:Ntimesteps, "", sep="")
  
  # Start cluster to run parallel.
  cl_time <- makeCluster(min(ncores, detectCores()))
  
  # For tests.
  #envir_cond1_df[, 3] <- evaluate_landscape(timesteps_names[1], habitability, evaluate_cond1, sd=0)
  #envir_cond1_df[, 1] <- coordinates$long
  #envir_cond1_df[, 2] <- coordinates$lat
  #envir_cond1_raster <- rasterFromXYZ(envir_cond1_df[,c(1,2,3)])
  #plot(envir_cond1_raster)
  
  # Use parallel Sapply for temporal loop - time uncorrelated.
  envir_cond1_df[,3:(Ntimesteps+2)] <- parSapply(cl_time, timesteps_names, evaluate_landscape, 
                                                 habitability, evaluate_cond1, sd=sd)
  envir_cond2_df[,3:(Ntimesteps+2)] <- parSapply(cl_time, timesteps_names, evaluate_landscape, 
                                                 habitability, evaluate_cond2, sd=sd)
  # Close the cluster.
  stopCluster(cl_time)
  
  #for(t in 1:Ntimesteps){
  #  envir_cond3_df[, t+2] <- landscape_cost$cost
  #}
  
  # Associating coordinates to dataframes.
  envir_cond1_df[, 1] <- habitability$long
  envir_cond1_df[, 2] <- habitability$lat
  envir_cond2_df[, 1] <- habitability$long
  envir_cond2_df[, 2] <- habitability$lat
  #envir_cond3_df[, 1] <- habitability$long
  #envir_cond3_df[, 2] <- habitability$lat
  
  # Add dataframes to list.
  envir_cond$cond1 <- envir_cond1_df
  envir_cond$cond2 <- envir_cond2_df
  #envir_cond$landscape_cost <- envir_cond3_df
  
  # Calculate processing time.
  envir_cond$processing_time <- proc.time() - time_start
  print(envir_cond$processing_time)
  
  return(envir_cond)
}

####################################################################################

####################################################################################
# Creating landscape.
#####################
print("Evaluating environmental conditions. 2/6")

envir_cond <- create_landscape_environment_2Dpar(habitability=habitability,
                                                 landscape_cost=landscape_cost,
                                                 dims=c(Nrow, Ncol, Ntimesteps), 
                                                 ncores=8)

# Plot environmental conditions.
#################################
print("Ploting environmental conditions. 3/6")

envir_raster_1 <- rasterFromXYZ(envir_cond$cond1[, c(1, 2, 3)])
envir_raster_2 <- rasterFromXYZ(envir_cond$cond2[, c(1, 2, 3)])
#envir_raster_3 <- rasterFromXYZ(envir_cond$landscape_cost[, c(1, 2, 3)])

# Parameters to plot the landscape.
# cond1min <- min(values(envir_raster), na.rm=T)
# cond1max <- max(values(envir_raster), na.rm=T)
# breaks_1 <- c(0, seq(cond1min, cond1max, length.out=50))
# colorscale_1 <- c('#4040ff', heat.colors(50, rev=T)) 
# 
# cond2min <- min(values(envir_raster_2), na.rm=T)
# cond2max <- max(values(envir_raster_2), na.rm=T)
# breaks_2 <- c(0, seq(cond2min, cond2max, length.out=50))
# colorscale_2 <- c('#4040ff', terrain.colors(50, rev=T)) 

# Plotting landscape.
par(mfrow=c(1,2))
pdf("simple_landscape2.pdf")
plot(envir_raster_1, col=heat.colors(50, rev=T), xlab="Long", ylab="Lat", main="Condition 1")
plot(envir_raster_2, col=terrain.colors(50, rev=T), xlab="Long", ylab="Lat", main="Condition 2")
dev.off()

# Saving landscape.
################################################################################

print("Saving landscape. 4/6")

# Define the work directory.



# Create a list of rasters with landscape information.
landscapes_list <- list()
for (timestep in 1:Ntimesteps){
  # Create a raster from the landscape associated to each time step.
  envir_cond1_raster <- rasterFromXYZ(envir_cond$cond1[, c(1, 2, timestep+2)])
  envir_cond2_raster <- rasterFromXYZ(envir_cond$cond2[, c(1, 2, timestep+2)])
  #envir_cond3_raster <- rasterFromXYZ(envir_cond$landscape_cost[, c(1, 2, timestep+2)])
  
  # Add the raster of each time step to the list.
  landscapes_list$cond1 <- c(landscapes_list$cond1, envir_cond1_raster)
  landscapes_list$cond2 <- c(landscapes_list$cond2, envir_cond2_raster)
  #landscapes_list$landscape_cost <- c(landscapes_list$landscape_cost, envir_cond3_raster)
}

# Save the landscapes list as a rds file.
saveRDS(landscapes_list, file.path(simulation_path, "landscapes", "simple_landscapes_20231114.rds"))


################################################################################
# Defining a cost function to calculate distances between sites.

print("Evaluating geographic distances. 5/6")

# HABITABILITY (NA) FUNCTION COST: the non-habitable sites (NA) has a different weight to 
# dispersion.
# In this case it has a double weight. 
# source: the index of source site (origin).
# habitable_src: habitability of source site.
# dest: the index of destiny site.
# habitable_dest: habitability of destiny site.
cost_function_hab <- function(source, habitable_src, dest, habitable_dest){
  # If any of sites (origin or destiny) are non-habitable return 2 (double 
  # weight), if both sites are habitable return 1.
  if(!all(habitable_src, habitable_dest)){
    return(2)
  }
  else{
    return(1)
  }
}



#cost_function_landscape <- function(source, cost_src, dest, cost_dest){
#  cat(cost_src, cost_dest, "\n")
#  # Compute the total cost to dispersal from source to dest.
#  total_cost <- cost_src + cost_dest
#  return(total_cost)
#}

# NULL FUNCTION COST: all the sites are equivalents.
# source: the index of source site (origin).
# habitable_src: habitability of source site.
# dest: the index of destiny site.
# habitable_dest: habitability of destiny site.
cost_function_null <- function(source, habitable_src, dest, habitable_dest){
  # Return 1 for any pair of sites.
  return(1)
}



################################################################################
# Calculating the distances matrices between sites using the cost function.

# Importing data.
landscapes <- readRDS(file.path(simulation_path, "landscapes", "simple_landscapes_20231114.rds"))
landscapes$cond1[[1]]
# Using the rds format the data is imported as a list, then it's not needed modify it.


################################################################################
# Creating the distances to input landscape in the simulation.
{
  time_start <- Sys.time()
  create_input_landscape(landscapes = landscapes,
                         # The cost function to calculate the distances in the landscape.
                         cost_function = cost_function_hab,
                         # Directory name to save the files in.
                         output_directory = file.path(simulation_path, "landscapes", "distances"),
                         # All surrounding sites from a focal site: 4, 8 or 16 neighboring sites.
                         directions = 8,
                         calculate_full_distance_matrices = T
                         # Calculate full distance matrices or only surrounding sites distances.
                         # In the second case the distance between other sites is calculated
                         # at each time step of simulation. Then:
                         # TRUE: large storage and higher speed computing.
                         # FALSE: lower storage and lower speed computing.
  )
  time_end <- Sys.time()
  print(time_end - time_start) # 31.217 mins # 29.286 mins
}

print("Landscape created. Work finished. 6/6")
################################################################################
# Finally, it's really important documenting the landscape you created. It is 
# doing writing the METADATA.txt file within the landscape folder. It's 
# automatically created by the function `create_input_landscape`.

