# Installing Required Libraries (If Required)
install_packages <- FALSE 
if(install_packages) { 
  install.packages("drat")
  drat:::add("mrc-ide")
  install.packages("dde")
  install.packages("odin")
  install.packages("remotes")
  remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE)
  pkgbuild::check_build_tools()
}

# Load required libraries & source helper functions
library(odin); library(dplyr); library(ggplot2) 
source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R") ## output of this is a function "stoch_seir" that can run the specified model given some parameter inputs
options(dplyr.summarise.inform = FALSE)

# Specifying model parameters
params <- list(

  # Epi Params
  beta = 0.3,                             # transmission rate - beta/sigma is R0 here
  gamma = 0.25,                           # 1/gamma = duration of latent period
  sigma = 0.2,                            # 1/sigma = duration of infectious period
  population_size = 10^7,                 # size of population in Location A
  start_infections = 10,                  # starting number of infections
  
  # Flight Params
  num_flights_leaveA = 250,              # number of flights leaving A per day
  num_flights_AtoB = 10,                 # number of flights leaving A per day that fly to B
  num_flights_arriveB = 200,             # total number of flights (including those from A) that arrive at B per day
  capacity_per_flight = 200,              # capacity of a single flight
  
  # Shedding Params
  non_virus_shed = 10^11,                 # total amount of uninteresting NA shed per shedding event 
  ratio_virus_to_non_virus = 10^-6,       # ratio of interesting to uninteresting NA shed per shedding event (i.e. how much pathogen as a proportion of all NA shed)
  shedding_prop = 0.75,                   # proportion of infected individuals who actually shed NA of interest
  shedding_freq = 1,                      # rate of defecation per flight
  
  # Sequencing Params
  sample_volume = 1,                      # volume of wastewater sampled (note just a multiplicative term that doesn't affect relative abundance currently)
  flight_ww_volume = 200,                 # volume of pooled flight wastewater (note just a multiplicative term that doesn't affect relative abundance currently)
  met_bias = 1,                           # metagenomic bias term
  seq_tot = 10^10,                        # sequenced reads per day;
  
  # Misc Params
  dt = 0.2,                               # model timestep (in days)
  end = 500                               # time (in days) to run model for
)

# Sense checking parameters
perc_popA_flying <- 100 * (params$capacity_per_flight * params$num_flights_leaveA)/params$population_size # percentage of population in location A taking a flight each day
perc_flights_fromA_toB <- 100 * params$num_flights_AtoB/params$num_flights_leaveA
perc_flights_fromA_toB_totalB <- 100 * params$num_flights_AtoB/params$num_flights_arriveB
virus_shed_per_event <- params$non_virus_shed * params$ratio_virus_to_non_virus # Amount of virus shed per event

# Create instance of model with your chosen parameters
mod <- stoch_seir$new(
  
  # Epidemiological Parameters
  beta = params$beta, 
  gamma = params$gamma, 
  sigma = params$sigma, 
  population_size = params$population_size, 
  start_infections = params$start_infections,
  
  # Flight Parameters
  capacity_per_flight = params$capacity_per_flight, 
  num_flights_leaveA = params$num_flights_leaveA, 
  num_flights_AtoB = params$num_flights_AtoB, 
  num_flights_arriveB = params$num_flights_arriveB, 
  
  # Sequencing Parameters
  shedding_prop = params$shedding_prop,
  shedding_freq = params$shedding_freq, 
  virus_shed = params$non_virus_shed * params$ratio_virus_to_non_virus, 
  non_virus_shed = params$non_virus_shed, 
  met_bias = params$met_bias, 
  seq_tot = params$seq_tot, 
  samp_frac_aggFlight = params$sample_volume/(params$flight_ww_volume * params$num_flights_arriveB),
  
  # Miscellaenous Parameters
  dt = params$dt)

# Running the model
set.seed(1000)
output <- mod$run(1:(params$end/params$dt))
output <- mod$transform_variables(output)

# Sense Checking Epidemiological Output
plot(output$time, output$S, ylim = c(0, params$population_size), type = "l")
lines(output$time, output$E, ylim = c(0, params$population_size), col = "blue")
lines(output$time, output$I, ylim = c(0, params$population_size), col = "red")
lines(output$time, output$R, ylim = c(0, params$population_size), col = "orange")

# Aggregating Data to Daily Level
raw_df <- data.frame(time = output$time,
                     infectionsAtoB = output$n_inf_flight_AtoB_Out,
                     amount_virus = output$amount_virus_aggFlight_Out,
                     amount_non_virus = output$amount_non_virus_aggFlight_Out)
daily_df <- raw_df %>%
  dplyr::mutate(time2 = midpoints(cut(time, breaks = max(time)))) %>% 
  group_by(time2) %>%
  summarise(daily_infections = sum(infectionsAtoB),
            daily_amount_virus = sum(amount_virus),
            daily_amount_non_virus = sum(amount_non_virus),
            daily_relative_abundance = daily_amount_virus/(daily_amount_virus + daily_amount_non_virus),
            daily_reads_det = daily_relative_abundance * params$seq_tot,
            daily_reads_stoch = rbinom(n = 1, size = params$seq_tot, prob = daily_relative_abundance)) 

# Plotting output for model run
max_time_plotting <- which(daily_df$daily_reads_det == max(daily_df$daily_reads_det))
ggplot(daily_df, aes(x = time2)) +
  geom_bar(aes(y = daily_reads_stoch), stat = "identity", fill = adjustcolor("#F45B69", alpha.f = 0.65)) +
  geom_line(aes(y = daily_reads_det), col = "blue") +
  lims(x = c(1, max_time_plotting)) +
  theme_bw() +
  labs(x = "Time (Days)", y = "Daily Reads")

