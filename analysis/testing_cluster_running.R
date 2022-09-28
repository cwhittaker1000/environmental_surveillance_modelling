###################################################################################################
##                                                                                               ##
##                      Prepping and setting up the cluster connection                           ##
##                                                                                               ##
###################################################################################################

# Load required libraries
library(tidyverse); library(tictoc); library(parallel); library(patchwork); library(odin)

# Source helper functions
source("functions/helper_functions.R")
source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R")

### Epidemiological
stochastic_sim <- 150
fixed_params <- list(# Misc Params
                     num_reads = 5,
                     dt = 0.2, 
                     seed = rpois(stochastic_sim, 200) * rpois(stochastic_sim, 200) * rpois(stochastic_sim, 20),
                     stochastic_sim = stochastic_sim, 
                     start_infections = 10, 
                     
                     # Epi Params
                     gamma = 0.25,
                     sigma = 0.2,
                     population_size = 10^6, 
                     
                     # Flight Params
                     proportion_AB = 0.1, # proportion flights going A -> B
                     vol_flight_ww = 800, 
                     sample_flight_ww = 1, 
                     capacity_per_flight = 200, 
                     
                     # Shedding Params
                     non_virus_shed = 2000000, 

                     # Sequencing Params
                     met_bias = 1)

# Epi Params - Beta
R0_sens <- seq(1.5, 3, 0.25)        
beta_sens <- R0 * fixed_params$sigma 
fixed_params$end <- data.frame(end = c(250, 200, 150, 125, 125, 100, 100), beta = beta_sens) # computational trick - run model for less time when R0 is higher

# Flight Params - Num Flights (Prop A->B Assumed)
num_flights_sens <- seq(50, 200, 20)
num_flightsAB_sens <- num_flights_sens * fixed_params$proportion_AB
prop_flying_sens <- num_flights_sens * fixed_params$capacity_per_flight / fixed_params$population_size
prop_flyingAB <- prop_flying_sens * fixed_params$proportion_AB

# Shedding Params - Frequency & Ratio Virus:Everything Else
shedding_sens <- seq(0.4, 2, 0.4)
ratio_sens <- lseq(1e-07, 1e-04, 10)   

# Sequencing Params - Total Sequencing Done
seq_tot_sens <- round(lseq(10^7, 10^9, 10))

# Create Overall Set of Parameter Values
variable_params  <- expand.grid(beta_sens = beta_sens, 
                                num_flights_sens = num_flights_sens, 
                                shedding_sens = shedding_sens, 
                                ratio_sens = ratio_sens, 
                                seq_tot_sens = seq_tot_sens)


cluster_model_running <- function(fixed_params, variable_params, generator) {
  
  # Setting Up Storage for Model Outputs
  variable_params$num_reads <- fixed_params$num_reads
  variable_params$ttd_lower <- NA
  variable_params$ttd_mean <- NA
  variable_params$ttd_upper <- NA
  variable_params$cuminf_lower <- NA
  variable_params$cuminf_mean <- NA
  variable_params$cuminf_upper <- NA
  variable_params$avg_inf <- NA
  variable_params$lower_inf <- NA
  variable_params$upper_inf <- NA
  variable_params$avg_flight_prev <- NA
  variable_params$lower_flight_pre <- NA
  variable_params$upper_flight_pre <- NA
  variable_params$perc_success <- NA
  
  # Setting Up the Model
  for (i in 1:dim(variable_params)[1]) {
    
    # Loading in set of parameters
    beta <- variable_params$beta_sens[i]
    shedding_freq <- variable_params$shedding_sens[i]
    ratio_virus_to_non_virus <- variable_params$ratio_sens[i]
    num_flights <- variable_params$num_flights_sens[i]
    seq_tot <- variable_params$seq_tot_sens[i]
    
    ##### RESTART WORKING FROM HERE ######
    
    # Generating the model instance
    mod <- generator$new(
      
      # Epidemiological Parameters
      beta = beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
      population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
      
      # Flight Parameters
      capacity_per_flight = fixed_params$capacity_per_flight, 
      num_flights = num_flights, num_flightsAB = num_flights * fixed_params$proportion_AB, 
      samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * num_flightsAB),
      
      # Sequencing Parameters
      shedding_freq = shedding_freq, virus_shed = fixed_params$non_virus_shed * ratio_virus_to_non_virus, 
      non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
      seq_tot = seq_tot, 
      
      # Miscellaenous Parameters
      dt = dt)
    
    # Setting up temporary storage for stochastic realisation results
    temp_storage <- matrix(data = NA, nrow = fixed_params$stochastic_sim, ncol = 5)
    colnames(temp_storage) <- c("ttd", "num_infs", "flight_prev", "cuminf", "success")
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Setting Seed
      set.seed(fixed_params$seed[j])
      
      # Running the Model
      end <- fixed_params$end$end[which(fixed_params$end$beta == beta)]
      output <- mod$run(1:(end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Extracting Outputs
      temp_storage[j, "ttd"] <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      temp_storage[j, "cuminf"] <- if(!is.na(temp_storage[j, "ttd"])) sum(output2$n_SE_Output[1:(temp_storage[j, "ttd"]/dt)])/fixed_params$population_size else NA
      temp_storage[j, "success"] <- if(!is.na(temp_storage[j, "cuminf"])) 1 else 0
    }
    
    variable_params$ttd_lower[i] <- min(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$ttd_mean[i] <- mean(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$ttd_upper[i] <- max(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$cuminf_lower[i] <- min(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$cuminf_mean[i] <- mean(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$cuminf_upper[i] <- max(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$perc_success[i] <- 100 * sum(temp_storage[, "success"]) / stochastic_sim
    
    rm(temp_storage)
    print(i) 
  }
  
  # Saving and Writing the Output for the MCMC Chains
  saveRDS(variable_params, file = "cluster_test_output.rds")
  
  # Returning the Output
  return(variable_params)
  
}





# Sets Working Directory to a Network Share
setwd("Q:/")

# Initial configuration of the cluster
didehpc::didehpc_config(home = "Q:/")

# Functions to be Run on the Cluster (Functions that the Wrapper Function calls)
cluster_function <- "Q:/nao_modelling/functions/helper_functions.R"
sources <- c(cluster_function)

# Create A Context Describing the Content You Want the Cluster to Utilise
#       Given a combination of R packages and function definitions that you provide
#       this sets up a copy of your work environment on the cluster.
general_context_directory <- "Q:/"
specific_context_directory <- "nao_modelling_test_26thSeptember2022" # change data as appropriate
context_directory <- paste0(general_context_directory, specific_context_directory)
ctx <- context::context_save(context_directory, # Directory that everything will be saved in
                             sources = sources, # Source files for functions called by the wrapper function
                             packages = c("odin")) # Packages required to run wrapper function (compiled code has to be in the form of a package, hence EPILOAModelling) 

# Configure the queue 
config <- didehpc::didehpc_config(rtools = TRUE, home = "Q:/", cluster = "dideclusthn") #  "fi--dideclusthn" is the small cluster. "fi--didemrchnb" is the big cluster- use when the small one is swamped/busy

# Create the Queue- represents the interface to the actual cluster queue should install all the packages required for the context
run <- didehpc::queue_didehpc(ctx, config = config) 

# Summary of all available clusters and checking various tasks
run$cluster_load(nodes = FALSE)
run$task_list()
run$task_times()

# Test Run
test_run <- run$enqueue(cluster_model_running(fixed_params, variable_params[100:102, ], stoch_seir_dust))
test_run$status()

test_run$context_id()
test_run$times() 
test_run$log()
bloop <- test_run$result()
saveRDS(bloop, "test_cluster_running_result.rds")
