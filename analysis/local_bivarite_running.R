# Load required libraries
library(tidyverse); library(odin); library(parallel); library(scales); library(tictoc)

# Source helper functions
source("functions/helper_functions.R")
source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R")
options(dplyr.summarise.inform = FALSE)

# Fixed Params for Model Running
stochastic_sim <- 200
fixed_params <- list(
  # Misc Params
  num_reads = 5,
  dt = 0.2, 
  seed = rpois(stochastic_sim, 200) * rpois(stochastic_sim, 10000) * rpois(stochastic_sim, 20),
  stochastic_sim = stochastic_sim, 
  start_infections = 10,
  end = 250,
  
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
R0_sens <- seq(1.5, 3.5, length.out = 25)        
beta_sens <- R0_sens * fixed_params$sigma 
beta_fixed <- beta_sens[R0_sens == 2.5]

# Flight Params - Num Flights (Prop A->B Assumed)
num_flights_sens <- seq(20, 1000, 10)
num_flightsAB_sens <- num_flights_sens * fixed_params$proportion_AB
prop_flying_sens <- num_flights_sens * fixed_params$capacity_per_flight / fixed_params$population_size
prop_flyingAB <- prop_flying_sens * fixed_params$proportion_AB
num_flights_fixed <- num_flights_sens[round(length(num_flights_sens)/2)]

# Shedding Params - Frequency & Ratio Virus:Everything Else
shedding_sens <- seq(0.4, 2.5, 0.1)
shedding_fixed <- 1
ratio_sens <- lseq(1e-07, 1e-04, 50)   
ratio_fixed <- 5e-06

# Sequencing Params - Total Sequencing Done
seq_tot_sens <- round(lseq(10^7, 10^9, 50))
seq_tot_fixed <- 10^8

# Create Overall Set of Parameter Values and Run Bivariate Sensitivity Analyses

## R0 vs SeqTotal Exploration
R0_seqTotal_params <- expand.grid(beta = beta_sens, seq_tot = seq_tot_sens)
R0_seqTotal_params$shedding_freq <- shedding_fixed
R0_seqTotal_params$ratio_virus_to_non_virus <- ratio_fixed
R0_seqTotal_params$num_flights <- num_flights_fixed
R0_seqTotal_params_list <- vector(mode = "list", length = dim(R0_seqTotal_params)[1])
for (i in 1:length(R0_seqTotal_params_list)) {
  R0_seqTotal_params_list[[i]]$beta <- R0_seqTotal_params$beta[i]
  R0_seqTotal_params_list[[i]]$seq_tot <- R0_seqTotal_params$seq_tot[i]
  R0_seqTotal_params_list[[i]]$shedding_freq <- R0_seqTotal_params$shedding_freq[i]
  R0_seqTotal_params_list[[i]]$ratio_virus_to_non_virus <- R0_seqTotal_params$ratio_virus_to_non_virus[i]
  R0_seqTotal_params_list[[i]]$num_flights <- R0_seqTotal_params$num_flights[i]
}
numCores <- 8 
cluster <- makeCluster(numCores)
clusterEvalQ(cluster, {
  library(Rcpp)
  library(odin)
  library(dplyr)
})
tic()
clusterEvalQ(cluster, source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"))
clusterEvalQ(cluster, source("functions/helper_functions.R"))
R0_seqTotal_output <- wrapped_parallel(variable_params_list = R0_seqTotal_params_list,
                                       fixed_params = fixed_params, 
                                       generator = stoch_seir_dust,
                                       filename = "R0_seqTotal_bivariate_sensitivity_analysis.rds",
                                       cluster = cluster)
toc()
stopCluster(cluster)

## Virus Ratio vs SeqTotal Exploration
ratio_seqTotal_params <- expand.grid(ratio_virus_to_non_virus = ratio_sens, seq_tot = seq_tot_sens)
ratio_seqTotal_params$shedding_freq <- shedding_fixed
ratio_seqTotal_params$beta <- beta_fixed
ratio_seqTotal_params$num_flights <- num_flights_fixed
ratio_seqTotal_params_list <- vector(mode = "list", length = dim(ratio_seqTotal_params)[1])
for (i in 1:length(ratio_seqTotal_params_list)) {
  ratio_seqTotal_params_list[[i]]$beta <- ratio_seqTotal_params$beta[i]
  ratio_seqTotal_params_list[[i]]$seq_tot <- ratio_seqTotal_params$seq_tot[i]
  ratio_seqTotal_params_list[[i]]$shedding_freq <- ratio_seqTotal_params$shedding_freq[i]
  ratio_seqTotal_params_list[[i]]$ratio_virus_to_non_virus <- ratio_seqTotal_params$ratio_virus_to_non_virus[i]
  ratio_seqTotal_params_list[[i]]$num_flights <- ratio_seqTotal_params$num_flights[i]
}
numCores <- 8 
cluster <- makeCluster(numCores)
clusterEvalQ(cluster, {
  library(Rcpp)
  library(odin)
  library(dplyr)
})
clusterEvalQ(cluster, source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"))
clusterEvalQ(cluster, source("functions/helper_functions.R"))
ratio_seqTotal_output <- wrapped_parallel(variable_params_list = ratio_seqTotal_params_list,
                                          fixed_params = fixed_params, 
                                          generator = stoch_seir_dust,
                                          filename = "ratio_seqTotal_bivariate_sensitivity_analysis.rds",
                                          cluster = cluster)
stopCluster(cluster)

## Num Flights vs SeqTotal Exploration
flights_seqTotal_params <- expand.grid(num_flights = num_flights_sens, seq_tot = seq_tot_sens)
flights_seqTotal_params$shedding_freq <- shedding_fixed
flights_seqTotal_params$beta <- beta_fixed
flights_seqTotal_params$ratio_virus_to_non_virus <- ratio_fixed
flights_seqTotal_params_list <- vector(mode = "list", length = dim(flights_seqTotal_params)[1])
for (i in 1:length(flights_seqTotal_params_list)) {
  flights_seqTotal_params_list[[i]]$beta <- flights_seqTotal_params$beta[i]
  flights_seqTotal_params_list[[i]]$seq_tot <- flights_seqTotal_params$seq_tot[i]
  flights_seqTotal_params_list[[i]]$shedding_freq <- flights_seqTotal_params$shedding_freq[i]
  flights_seqTotal_params_list[[i]]$ratio_virus_to_non_virus <- flights_seqTotal_params$ratio_virus_to_non_virus[i]
  flights_seqTotal_params_list[[i]]$num_flights <- flights_seqTotal_params$num_flights[i]
}
numCores <- 8 
cluster <- makeCluster(numCores)
clusterEvalQ(cluster, {
  library(Rcpp)
  library(odin)
  library(dplyr)
})
clusterEvalQ(cluster, source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"))
clusterEvalQ(cluster, source("functions/helper_functions.R"))
flights_seqTotal_output <- wrapped_parallel(variable_params_list = flights_seqTotal_params_list,
                                           fixed_params = fixed_params, 
                                           generator = stoch_seir_dust,
                                           filename = "flights_seqTotal_bivariate_sensitivity_analysis.rds",
                                           cluster = cluster)
stopCluster(cluster)
