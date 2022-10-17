# To run on Mac, if you've never used R:
#   $ brew install R
#   $ echo 'options(repos=structure(c(CRAN="https://cran.r-project.org")))' \
#           >> ~/.Rprofile
#   $ Rscript ./analysis/jeff_example.R
#   $ open Rplots.pdf

# Installing Required Libraries (If Required)
install_packages <- FALSE 
if(install_packages) { 
  install.packages("drat")
  drat:::add("mrc-ide")
  install.packages("pkgbuild")
  install.packages("pkgload")
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("cowplot")
  install.packages("dde")
  install.packages("odin")
  pkgbuild::check_build_tools()
}

# Load required libraries & load the model
library(odin); library(dplyr); library(ggplot2); library(cowplot)
source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R") ## output of this is a function "stoch_seir" that can 
                                                                    ## run the specified model given some parameter inputs
source("functions/helper_functions.R")
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
end_time <- params$end/params$dt
output <- mod$run(1:end_time)
output <- mod$transform_variables(output)

# Aggregating Data to Daily Level
raw_df <- data.frame(time = output$time,
                     new_infections = output$n_EI_Output,
                     currently_infectedA = output$I,
                     infectionsAtoB = output$n_inf_flight_AtoB_Out,
                     amount_virus = output$amount_virus_aggFlight_Out,
                     amount_non_virus = output$amount_non_virus_aggFlight_Out)
daily_df <- raw_df %>%
  dplyr::mutate(time2 = midpoints(cut(time, breaks = max(time)))) %>% 
  group_by(time2) %>%
  summarise(daily_infections = sum(new_infections),
            daily_infections_AtoB = sum(infectionsAtoB),
            daily_prevalence_infection = 100 * mean(currently_infectedA)/params$population_size,
            daily_amount_virus = sum(amount_virus),
            daily_amount_non_virus = sum(amount_non_virus),
            daily_relative_abundance = daily_amount_virus/(daily_amount_virus + daily_amount_non_virus),
            daily_reads_det = daily_relative_abundance * params$seq_tot,
            daily_reads_stoch = rbinom(n = 1, size = params$seq_tot, prob = daily_relative_abundance)) %>%
  mutate(cumulative_incidence = 100 * cumsum(daily_infections)/params$population_size)

# Plotting output for model run
max_time_plotting <- which(daily_df$daily_reads_det == max(daily_df$daily_reads_det))

reads_plot <- ggplot(daily_df, aes(x = time2)) +
  geom_bar(aes(y = daily_reads_stoch), stat = "identity", fill = adjustcolor("#F45B69", alpha.f = 0.65)) +
  geom_line(aes(y = daily_reads_det), col = "blue") +
  lims(x = c(1, max_time_plotting)) +
  theme_bw() +
  labs(x = "Time (Days)", y = "Daily Reads")

ever_inf_plot <- ggplot(daily_df, aes(x = time2)) +
  geom_line(aes(y = cumulative_incidence), col = "blue") +
  lims(x = c(1, max_time_plotting), y = c(0, daily_df$cumulative_incidence[max_time_plotting])) +
  theme_bw() +
  labs(x = "Time (Days)", y = "Cumulative Incidence (% Ever Infected)")

curr_inf_plot <- ggplot(daily_df, aes(x = time2)) +
  geom_line(aes(y = daily_prevalence_infection), col = "blue") +
  lims(x = c(1, max_time_plotting), y = c(0, daily_df$daily_prevalence_infection[max_time_plotting])) +
  theme_bw() +
  labs(x = "Time (Days)", y = "Daily Prevalence Of Infection (%)")

plot_grid(reads_plot, ever_inf_plot, curr_inf_plot, nrow = 1, rel_widths = c(2, 1, 1))

# Running the model multiple times
n_iter <- 10
seeds <- rpois(n = n_iter, lambda = 10^7)
output_storage <- vector(mode = "list", length = n_iter)
for (i in 1:n_iter) {
  set.seed(seeds[i])
  
  output <- mod$run(1:end_time)
  output <- mod$transform_variables(output)
  
  raw_df <- data.frame(time = output$time,
                       new_infections = output$n_EI_Output,
                       currently_infectedA = output$I,
                       infectionsAtoB = output$n_inf_flight_AtoB_Out,
                       amount_virus = output$amount_virus_aggFlight_Out,
                       amount_non_virus = output$amount_non_virus_aggFlight_Out)
  daily_df <- raw_df %>%
    dplyr::mutate(time2 = midpoints(cut(time, breaks = max(time)))) %>% 
    group_by(time2) %>%
    summarise(daily_infections = sum(new_infections),
              daily_infections_AtoB = sum(infectionsAtoB),
              daily_prevalence_infection = 100 * mean(currently_infectedA)/params$population_size,
              daily_amount_virus = sum(amount_virus),
              daily_amount_non_virus = sum(amount_non_virus),
              daily_relative_abundance = daily_amount_virus/(daily_amount_virus + daily_amount_non_virus),
              daily_reads_det = daily_relative_abundance * params$seq_tot,
              daily_reads_stoch = rbinom(n = 1, size = params$seq_tot, prob = daily_relative_abundance)) %>%
    mutate(cumulative_incidence = 100 * cumsum(daily_infections)/params$population_size)
  daily_df$iteration <- i
  
  output_storage[[i]] <- daily_df
  
}

multi_run_output <- bind_rows(output_storage) %>%
  relocate(iteration) %>%
  mutate(iteration = paste0("Iteration ", iteration))
multi_run_lines <- ggplot(multi_run_output, aes(x = time2, col = factor(iteration))) +
  geom_line(aes(y = daily_reads_det), alpha = 0.5) +
  labs(x = "Time (Days)", y = "Daily Reads (Continuous)") +
  theme(legend.position = "none")
multi_run_bars <- ggplot(multi_run_output, aes(x = time2, fill = factor(iteration))) +
  geom_bar(aes(y = daily_reads_stoch), stat = "identity") +
  facet_wrap(~factor(iteration), scales = "free_y") +
  lims(x = c(50, 225)) +
  labs(x = "Time (Days)", y = "Daily Reads (Counts)") +
  theme(legend.position = "none")

plot_grid(multi_run_lines, multi_run_bars, rel_widths = c(1, 2))
