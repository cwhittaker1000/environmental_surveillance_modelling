###################################################################################################
##                                                                                               ##
##                      Prepping and setting up the cluster connection                           ##
##                                                                                               ##
###################################################################################################

# Load required libraries
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel); library(patchwork)

# Source helper functions
source("functions/helper_functions.R")
source("models/temp_model_file_cluster.R")

### Epidemiological
R0 <- seq(1.5, 3, 0.25)         # Range 1.50 to 3, increments of 0.25, 7 values
fixed_params <- list(num_reads = 5,
                     dt = 0.2, 
                     seed = rpois(100, 200) * rpois(100, 200) * rpois(100, 20),
                     stochastic_sim = 50, 
                     vol_flight_ww = 800, 
                     sample_flight_ww = 1, 
                     population_size = 10^6, 
                     capacity_per_flight = 250, 
                     non_virus_shed = 2000000, 
                     start_infections = 10, 
                     met_bias = 1,
                     gamma = 1/4,
                     sigma = 1/5)
beta <- R0 * fixed_params$sigma 
fixed_params$end = data.frame(end = c(250, 200, 150, 125, 125, 100, 100), beta = beta)

## Parameters To Vary

### Shedding 
shedding_freq <- seq(0.2, 1.8, 0.4)                  # Range 0.2 to 1.8, increments of 0.4 - 5 values
ratio_virus_to_non_virus <- 1/round(lseq(10^6, 10^8, 7))   # Range 1/10^6 to 1/10^8 - 8 values

### Flight
proportion_flying <- seq(0.002, 0.010, 0.002)       # Range 0.002 to 0.010 (i.e. 0.2% to 1% population flying daily), increments of 0.002 - 5 values
num_flights <- round(proportion_flying * fixed_params$population_size / fixed_params$capacity_per_flight, 0)          
proportion_AB <- seq(0.005, 0.025, 0.005)           # Range 0.005 to 0.025 (i.e. 1 in every 200 flights to 1 in every 40 flights), increment of 0.005 - 5 values
num_flightsAB <- num_flights * proportion_AB
function_pop <- proportion_flying * proportion_AB

### Sequencing
seq_tot <- round(lseq(10^7, 10^9, 8)) # Range TBD but 8 values

# Create Overall Set of Parameter Values
variable_params  <- expand.grid(beta = beta, 
                                shedding_freq = shedding_freq, 
                                ratio_virus_to_non_virus = ratio_virus_to_non_virus, 
                                num_flights = num_flights, 
                                num_flightsAB = num_flightsAB, 
                                seq_tot = seq_tot)


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
