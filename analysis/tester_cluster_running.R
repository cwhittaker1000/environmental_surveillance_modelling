###################################################################################################
##                                                                                               ##
##                      Prepping and setting up the cluster connection                           ##
##                                                                                               ##
###################################################################################################

# Sets Working Directory to a Network Share
setwd("Q:/")

# Initial configuration of the cluster
didehpc::didehpc_config(home = "Q:/")

# Functions to be Run on the Cluster (Functions that the Wrapper Function calls)
cluster_function <- "Q:/environmental_surveillance_modelling/functions/helper_functions.R"
model_call <- "Q:/environmental_surveillance_modelling/models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"
sources <- c(cluster_function, model_call)

# Create A Context Describing the Content You Want the Cluster to Utilise
#       Given a combination of R packages and function definitions that you provide
#       this sets up a copy of your work environment on the cluster.
general_context_directory <- "Q:/"
specific_context_directory <- "nao_modelling_test_1stOctober2022" # change date as appropriate
context_directory <- paste0(general_context_directory, specific_context_directory)
pkg_src <- conan::conan_sources("mrc-ide/odin")
ctx <- context::context_save(context_directory, # Directory that everything will be saved in
                             sources = sources, # Source files for functions called by the wrapper function
                             packages = c("odin", "parallel", "dplyr", "scales", "pkgbuild", "pkgload", "dde"),
                             package_sources = pkg_src) 

# Load required libraries
library(tidyverse); library(odin); library(parallel); library(scales)

# Source helper functions
setwd("environmental_surveillance_modelling/")
source("functions/helper_functions.R")
options(dplyr.summarise.inform = FALSE)

# Fixed Params for Model Running
stochastic_sim <- 200
fixed_params <- list(# Misc Params
  num_reads = 5,
  dt = 0.2, 
  seed = rpois(stochastic_sim, 200) * rpois(stochastic_sim, 10000) * rpois(stochastic_sim, 20),
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
R0_sens <- seq(1.5, 3.5, 0.25)        
beta_sens <- R0_sens * fixed_params$sigma 
fixed_params$end <- data.frame(end = c(250, 200, 150, 125, 125, 100, 100, 100, 100), beta = beta_sens) # computational trick - run model for less time when R0 is higher

# Flight Params - Num Flights (Prop A->B Assumed)
num_flights_sens <- seq(50, 250, 20)
num_flightsAB_sens <- num_flights_sens * fixed_params$proportion_AB
prop_flying_sens <- num_flights_sens * fixed_params$capacity_per_flight / fixed_params$population_size
prop_flyingAB <- prop_flying_sens * fixed_params$proportion_AB

# Shedding Params - Frequency & Ratio Virus:Everything Else
shedding_sens <- seq(0.4, 2, 0.2)
ratio_sens <- lseq(1e-07, 1e-04, 12)   

# Sequencing Params - Total Sequencing Done
seq_tot_sens <- round(lseq(10^7, 10^9, 12))

# Create Overall Set of Parameter Values
variable_params  <- expand.grid(beta_sens = beta_sens, 
                                num_flights_sens = num_flights_sens, 
                                shedding_sens = shedding_sens, 
                                ratio_sens = ratio_sens, 
                                seq_tot_sens = seq_tot_sens)

variable_params_list <- vector(mode = "list", length = dim(variable_params)[1])
for (i in 1:length(variable_params_list)) {
  variable_params_list[[i]]$beta <- variable_params$beta_sens[i]
  variable_params_list[[i]]$shedding_freq <- variable_params$shedding_sens[i]
  variable_params_list[[i]]$ratio_virus_to_non_virus <- variable_params$ratio_sens[i]
  variable_params_list[[i]]$num_flights <- variable_params$num_flights_sens[i]
  variable_params_list[[i]]$seq_tot <- variable_params$seq_tot_sens[i]
}

## 12 Cores

# Configure the queue 
config <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", template = "32Core" ) #  "fi--dideclusthn" is the small cluster. "fi--didemrchnb" is the big cluster- use when the small one is swamped/busy

# Create the Queue- represents the interface to the actual cluster queue should install all the packages required for the context
run <- didehpc::queue_didehpc(ctx, config = config) 

# Summary of all available clusters and checking various tasks
run$cluster_load(nodes = FALSE)
run$task_list()
run$task_times()

# 32 Core Big Job Running
test_run32 <- run$enqueue(wrapped_parallel(variable_params_list, fixed_params, stoch_seir_dust, NULL))
test_run32$status()
x <- test_run32$result()

saveRDS(x, file = "Q:/environmental_surveillance_modelling/outputs/full_sensitivity_analysis_cluster.rds")


# Test Run
test_run12 <- run$enqueue(wrapped_parallel(variable_params_list, fixed_params, stoch_seir_dust, NULL))
test_run12$status()
x <- test_run12$result()

saveRDS(x, file = "Q:/environmental_surveillance_modelling/outputs/full_sensitivity_analysis_cluster.rds")

## 8 Cores

# Configure the queue 
config8 <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", cores = 8) # template = "32Core" ) #  "fi--dideclusthn" is the small cluster. "fi--didemrchnb" is the big cluster- use when the small one is swamped/busy

# Create the Queue- represents the interface to the actual cluster queue should install all the packages required for the context
run8 <- didehpc::queue_didehpc(ctx, config = config8) 

# Summary of all available clusters and checking various tasks
run8$cluster_load(nodes = FALSE)
run8$task_list()
run8$task_times()

# Test Run
test_run8 <- run8$enqueue(wrapped_parallel(variable_params_list, fixed_params, stoch_seir_dust, NULL))
test_run8$status()
test_run8$result()

## 4 Cores

# Configure the queue 
config4 <- didehpc::didehpc_config(home = "Q:/", cluster = "fi--didemrchnb", cores = 4) # template = "32Core" ) #  "fi--dideclusthn" is the small cluster. "fi--didemrchnb" is the big cluster- use when the small one is swamped/busy

# Create the Queue- represents the interface to the actual cluster queue should install all the packages required for the context
run4 <- didehpc::queue_didehpc(ctx, config = config4) 

# Summary of all available clusters and checking various tasks
run4$cluster_load(nodes = FALSE)
run4$task_list()
run4$task_times()

# Test Run
test_run4 <- run4$enqueue(wrapped_parallel(variable_params_list, fixed_params, stoch_seir_dust, NULL))
test_run4$status()
test_run4$result()


