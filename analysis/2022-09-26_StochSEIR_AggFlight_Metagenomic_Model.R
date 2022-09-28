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
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel); library(patchwork)
source("functions/helper_functions.R")
options(dplyr.summarise.inform = FALSE)

# Stochastic SEIR
stoch_seir_dust <- odin::odin({
  
  ################## Epidemiological Model ##################
  
  ## Epidemiological Parameters
  beta <- user()              # probability of a contact successfully transmitting the disease
  gamma <- user()             # rate of transition from Exposed -> Infectious (incubation period)
  sigma <- user()             # rate of transition from Infectious -> Recovered (rate of recovery)
  population_size <- user()   # overall size of population
  start_infections <- user()  # starting number of infections (in Exposed compartment)
  
  ## Converting Epidemiological Rates to Probabilities of Leaving Each Compartment
  lambda <- ((beta * I) / N) 
  p_SE <- 1 - exp(-lambda* dt) # S to E 
  p_EI <- 1 - exp(-gamma * dt) # E to I 
  p_IR <- 1 - exp(-sigma * dt) # I to R
  
  ### Binomial Draws for Number of People Moving Between Compartments
  n_SE <- rbinom(S, p_SE)    # number of individuals infected at each timestep (S->E)
  n_EI <- rbinom(E, p_EI)    # number of individuals becoming infectious at each timestep (E->I)
  n_IR <- rbinom(I, p_IR)    # number of individuals recovering at each timestep (I->R)
  
  ### Stochastic Model Updates for Epidemiological States (S, E, I & R) & Quantities of Interest (new infectious)
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR - n_inf_flight # change this
  update(R) <- R + n_IR         # NOTE: need to eventually return n_inf_flight to here - still TODO
  update(N) <- S + E + I + R
  update(n_SE_Output) <- n_SE   # odin.dust doesn't let you output calculated quantities like n_EI & n_IR without making 
  update(n_EI_Output) <- n_EI   # odin.dust doesn't let you output calculated quantities like n_EI & n_IR without making 
  update(n_IR_Output) <- n_IR   # them tracked variables, which requires having initial() and update() calls for them
  
  ## Initial Values for Epidemiological States & Quantities of Interest
  initial(S) <- population_size - start_infections  
  initial(E) <- start_infections                    
  initial(I) <- 0
  initial(R) <- 0
  initial(N) <- S + E + I + R
  initial(n_SE_Output) <- 0
  initial(n_EI_Output) <- 0
  initial(n_IR_Output) <- 0
  
  ################## Airport/Airplane Surveillance Model ##################
  
  ## Airport/Airplane Parameters
  num_flights <- user()            # Number of flights per day
  num_flightsAB <- user()          # Number of flights per day from Location A to NAO Location (Location B)
  capacity_per_flight <- user()
  
  ## Calculating the number of infected people taking a flight based on number of infections and airport/airplane parameters
  n_inf_flight <- rhyper(I, population_size - I, dt * num_flights * capacity_per_flight) 
  n_inf_flightAB <- rhyper(n_inf_flight, num_flights * capacity_per_flight * dt - n_inf_flight, num_flightsAB * capacity_per_flight * dt) 

  ### Stochastic Model Updates for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  update(n_inf_flightOut) <- n_inf_flight
  update(n_inf_flightABOut) <- n_inf_flightAB

  ### Initial Values for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  initial(n_inf_flightOut) <- 0
  initial(n_inf_flightABOut) <- 0

  ################## Metagenomic Sequencing Model ##################
  
  ## Metagenomic and Sequencing Parameters
  shedding_freq <- user()              # Average number of defecation events per person per flight
  virus_shed <- user()                 # Average amount of material shed per event for our virus of interest (defecation)
  non_virus_shed <- user()             # Average amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  met_bias <- user()                   # Bias term for the metagenomic model
  seq_tot <- user()                    # Total amount of sequencing that is done
  samp_frac_aggFlight <- user()      # Fraction of all flight's wastewater that would be sampled if individual flight wastewater pooled and then sampled (via triturator) 
  
  ## AGGREGATED FLIGHT CALCULATIONS (E.G. SAMPLING FROM TRITURATOR)

  # Calculating the number of shedding events from infected and uninfected individuals on each airplane
  infected_indiv_shedding_events <- rpois(n_inf_flightAB * shedding_freq)
  uninfected_indiv_shedding_events <- rpois(((capacity_per_flight * num_flightsAB) - n_inf_flightAB) * shedding_freq)
  
  ### Calculating amount and concentration of nucleic acid shed into aggregated flight wastewater
  amount_virus_aggFlight <- infected_indiv_shedding_events * virus_shed 
  amount_non_virus_aggFlight <- (uninfected_indiv_shedding_events + infected_indiv_shedding_events) * non_virus_shed

  ### DETERMINISTIC - Converting this nucleic acid abundance into sequencing reads for aggregated flight wastewater
  sample_amount_virus_aggFlight_det <- amount_virus_aggFlight * samp_frac_aggFlight
  sample_amount_non_virus_aggFlight_det <- amount_non_virus_aggFlight * samp_frac_aggFlight
  seq_reads_virus_aggFlight_det <- if(sample_amount_virus_aggFlight_det == 0 && sample_amount_non_virus_aggFlight_det == 0) 0 else seq_tot * (sample_amount_virus_aggFlight_det * met_bias)/((sample_amount_virus_aggFlight_det * met_bias) + sample_amount_non_virus_aggFlight_det)
  seq_reads_non_virus_aggFlight_det <- seq_tot - seq_reads_virus_aggFlight_det
  
  ### STOCHASTIC - Converting this nucleic acid abundance into sequencing reads for aggregated flight wastewater
  sample_amount_virus_aggFlight_stoch <- rbinom(amount_virus_aggFlight, samp_frac_aggFlight)
  sample_amount_non_virus_aggFlight_stoch <- rbinom(amount_non_virus_aggFlight, samp_frac_aggFlight)
  seq_reads_virus_aggFlight_stoch <- if(sample_amount_virus_aggFlight_stoch == 0 && sample_amount_non_virus_aggFlight_stoch == 0) 0 else seq_tot * (sample_amount_virus_aggFlight_stoch * met_bias)/((sample_amount_virus_aggFlight_stoch * met_bias) + sample_amount_non_virus_aggFlight_stoch)
  seq_reads_non_virus_aggFlight_stoch <- seq_tot - seq_reads_virus_aggFlight_stoch

  ### Stochastic Model Updates for Outputs Relevant to Metagenomic Sequencing
  update(amount_virus_aggFlight_Out) <- amount_virus_aggFlight
  update(amount_non_virus_aggFlight_Out) <- amount_non_virus_aggFlight
  update(sample_amount_virus_aggFlight_det_Out) <- sample_amount_virus_aggFlight_det
  update(sample_amount_non_virus_aggFlight_det_Out) <- sample_amount_non_virus_aggFlight_det
  update(seq_reads_virus_aggFlight_det_Out) <- seq_reads_virus_aggFlight_det
  update(seq_reads_non_virus_aggFlight_det_Out) <- seq_reads_non_virus_aggFlight_det
  update(sample_amount_virus_aggFlight_stoch_Out) <- sample_amount_virus_aggFlight_stoch
  update(sample_amount_non_virus_aggFlight_stoch_Out) <- sample_amount_non_virus_aggFlight_stoch
  update(seq_reads_virus_aggFlight_stoch_Out) <- seq_reads_virus_aggFlight_stoch
  update(seq_reads_non_virus_aggFlight_stoch_Out) <- seq_reads_non_virus_aggFlight_stoch
  
  ### Initial Values for Outputs Relevant to Metagenomic Sequencing
  initial(amount_virus_aggFlight_Out) <- 0
  initial(amount_non_virus_aggFlight_Out) <- 0
  initial(sample_amount_virus_aggFlight_det_Out) <- 0
  initial(sample_amount_non_virus_aggFlight_det_Out) <- 0
  initial(seq_reads_virus_aggFlight_det_Out) <- 0
  initial(seq_reads_non_virus_aggFlight_det_Out) <- 0
  initial(sample_amount_virus_aggFlight_stoch_Out) <- 0
  initial(sample_amount_non_virus_aggFlight_stoch_Out) <- 0
  initial(seq_reads_virus_aggFlight_stoch_Out) <- 0
  initial(seq_reads_non_virus_aggFlight_stoch_Out) <- 0

  ################# Miscellaneous Model Stuff ##################
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt 
  
})

# Fixed Parameters 
stochastic_sim <- 150
fixed_params <- list(# Misc Params
                     num_reads = 5,
                     dt = 0.2, 
                     end = 250,
                     seed = rpois(stochastic_sim, 200) * rpois(stochastic_sim, 200) * rpois(stochastic_sim, 20),
                     stochastic_sim = stochastic_sim, 
                     start_infections = 10, 
                     
                     # Epi Params
                     beta = 0.6,
                     gamma = 0.25,
                     sigma = 0.2,
                     population_size = 10^6, 
                     
                     # Flight Params
                     num_flights = 40,
                     num_flightsAB = 20,
                     vol_flight_ww = 800, 
                     sample_flight_ww = 1, 
                     capacity_per_flight = 200, 
                     
                     # Shedding Params
                     non_virus_shed = 2000000, 
                     shedding_freq = 1,
                     ratio_virus_to_non_virus = 1/10^5,
                     
                     # Sequencing Params
                     met_bias = 1,
                     seq_tot = 10^9)

# Sense Check Parameters
R0_sc <- fixed_params$beta/fixed_params$sigma  # R0
flight_sc <- 100 * fixed_params$num_flightsAB * fixed_params$capacity_per_flight / fixed_params$population_size  # % of Pop Travelling A->B every day
virus_shed_sc <- fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus # Amount of virus shed per event

new_run <- TRUE
colours <- c("#88D18A", "#788BFF", "#5B5750", "#F45B69", "#F18F01")

#########################################################################
#####                    Beta Sensitivity Analysis                  #####
#########################################################################

# Generating values to vary beta over and dataframe to store outputs from model running
R0 <- seq(1.5, 3.5, 0.05)         
beta_sens <- R0 * fixed_params$sigma 
R0_df <- data.frame(beta = beta_sens, R0 = beta_sens/fixed_params$sigma)
beta_output <- data.frame(beta = rep(beta_sens, each = fixed_params$stochastic_sim), 
                          stochastic_realisation = 1:fixed_params$stochastic_sim, 
                          num_reads = fixed_params$num_reads, 
                          time_to_detection = NA,
                          flight_infections = NA,
                          flight_prevalence = NA,
                          community_prevalence = NA,
                          cumulative_incidence = NA)

# Running model and generating outputs for range of beta values
if (new_run) {
  counter <- 1
  for (i in 1:length(beta_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir_dust$new(
        
        # Epidemiological Parameters
        beta = beta_sens[i], gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB, 
        
        # Sequencing Parameters
        shedding_freq = fixed_params$shedding_freq, virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      ttd_metrics <- ttd_fun(mod, output2, fixed_params$num_reads)
      
      beta_output$time_to_detection[counter] <- ttd_metrics$time
      beta_output$flight_infections[counter] <- ttd_metrics$daily_flightAB_infections
      beta_output$flight_prevalence[counter] <- ttd_metrics$daily_flight_AB_prevalence
      beta_output$community_prevalence[counter] <- ttd_metrics$daily_community_prevalence
      beta_output$cumulative_incidence[counter] <- if(!is.na(ttd_metrics$time)) sum(output2$n_SE_Output[1:(ttd_metrics$time/fixed_params$dt)])/fixed_params$population_size else NA
      
      counter <- counter + 1   
    }
  }
  saveRDS(beta_output, file = "outputs/agg_beta_sensitivity_analysis.rds")
}  else {
  beta_output <- readRDS("outputs/agg_beta_sensitivity_analysis.rds")
}

# Summarising and plotting the output from beta sensitivity analysis
beta_df_summary <- beta_output %>%
  left_join(R0_df, by = "beta") %>%
  group_by(beta, R0) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections, na.rm = TRUE),
            lower_inf = quantile(flight_infections, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.95),
            lower_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection)))/fixed_params$stochastic_sim) %>%
  ungroup() %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         lower_inf = ifelse(is.infinite(lower_inf), NaN, lower_inf),
         upper_inf = ifelse(is.infinite(upper_inf), NaN, upper_inf),
         lower_flight_prev = ifelse(is.infinite(lower_flight_prev), NaN, lower_flight_prev),
         upper_flight_prev = ifelse(is.infinite(upper_flight_prev), NaN, upper_flight_prev),
         lower_comm_prev = ifelse(is.infinite(lower_comm_prev), NaN, lower_comm_prev),
         upper_comm_prev = ifelse(is.infinite(upper_comm_prev), NaN, upper_comm_prev))

beta_cuminf <- ggplot(data = beta_df_summary) +
  geom_line(aes(x = R0, y = avg_cuminf), col = colours[1]) +
  geom_ribbon(aes(x = R0, y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - R0", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

beta_ttd <- ggplot(data = beta_df_summary) +
  geom_line(aes(x = R0, y = avg_ttd), col = colours[1]) +
  geom_ribbon(aes(x = R0, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - R0", y = "Time to Detection (Days)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

beta_infs <- ggplot(data = beta_df_summary) +
  geom_line(aes(x = R0, y = avg_inf), col = colours[1]) +
  geom_ribbon(aes(x = R0, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - R0", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

beta_flight_prev <- ggplot(data = beta_df_summary) +
  geom_line(aes(x = R0, y = avg_flight_prev), col = colours[1]) +
  geom_ribbon(aes(x = R0, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - R0", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

beta_perc_reached <- ggplot(data = beta_df_summary) +
  geom_bar(aes(x = R0, y = perc_reached), col = "#4A4844",
           fill = adjustcolor(colours[1], alpha.f = 0.5), stat = "identity") +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  lims(y = c(0, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

beta_cuminf_plot <- beta_perc_reached + beta_cuminf +
  plot_layout(nrow = 2, heights = c(1, 6))
beta_ttd_plot <- beta_perc_reached + beta_ttd +
  plot_layout(nrow = 2, heights = c(1, 6))
beta_inf_plot <- beta_perc_reached + beta_infs +
  plot_layout(nrow = 2, heights = c(1, 6))
beta_prev_plot <- beta_perc_reached + beta_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 6))
beta_plot <- cowplot::plot_grid(beta_ttd_plot, beta_cuminf_plot, beta_inf_plot, beta_prev_plot, nrow = 1)

#########################################################################
#####               Seq Total Sensitivity Analysis                  #####
#########################################################################

# Generating values to vary seq total over and dataframe to store outputs from model running
seq_sens <- round(lseq(10^5, 10^9, 40))
seq_output <- data.frame(seq_total = rep(seq_sens, each = fixed_params$stochastic_sim), 
                         stochastic_realisation = 1:fixed_params$stochastic_sim, 
                         num_reads = fixed_params$num_reads, 
                         time_to_detection = NA,
                         cumulative_incidence = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(seq_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir_dust$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB, 
        
        # Sequencing Parameters
        shedding_freq = fixed_params$shedding_freq, virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
        seq_tot = seq_sens[i], samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      time_to_detection <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      seq_output$time_to_detection[counter] <- time_to_detection
      seq_output$cumulative_incidence[counter] <- if(!is.na(time_to_detection)) sum(output2$n_SE_Output[1:(time_to_detection/fixed_params$dt)])/fixed_params$population_size else NA
      
      counter <- counter + 1   
    }
  }
  saveRDS(seq_output, file = "outputs/agg_seqTotal_sensitivity_analysis.rds") 
} else {
  seq_output <- readRDS("outputs/agg_seqTotal_sensitivity_analysis.rds")
}

# Summarising and plotting the output from seq total sensitivity analysis
seq_df_summary <- seq_output %>%
  group_by(seq_total) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections, na.rm = TRUE),
            lower_inf = quantile(flight_infections, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.95),
            lower_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection)))/fixed_params$stochastic_sim) %>%
  ungroup() %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         lower_inf = ifelse(is.infinite(lower_inf), NaN, lower_inf),
         upper_inf = ifelse(is.infinite(upper_inf), NaN, upper_inf),
         lower_flight_prev = ifelse(is.infinite(lower_flight_prev), NaN, lower_flight_prev),
         upper_flight_prev = ifelse(is.infinite(upper_flight_prev), NaN, upper_flight_prev),
         lower_comm_prev = ifelse(is.infinite(lower_comm_prev), NaN, lower_comm_prev),
         upper_comm_prev = ifelse(is.infinite(upper_comm_prev), NaN, upper_comm_prev))

seq_cuminf <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_cuminf), col = colours[2]) +
  geom_ribbon(aes(x = seq_total, y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[2]) +
  labs(x = "Parameter Value - Seq Total", y = "Cumulative Incidence (%)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_ttd <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_ttd), col = colours[2]) +
  geom_ribbon(aes(x = seq_total, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[2]) +
  labs(x = "Parameter Value - Seq Total", y = "Time to Detection (Days)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_infs <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_inf), col = colours[1]) +
  geom_ribbon(aes(x = seq_total, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Seq Total", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_flight_prev <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_flight_prev), col = colours[1]) +
  geom_ribbon(aes(x = seq_total, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Seq Total", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_perc_reached <- ggplot(data = seq_df_summary) +
  geom_bar(aes(x = seq_total, y = perc_reached), col = "#4A4844",
           fill = adjustcolor(colours[2], alpha.f = 0.5), stat = "identity") +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  lims(y = c(0, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

seq_cuminf_plot <- seq_perc_reached + seq_cuminf +
  plot_layout(nrow = 2, heights = c(1, 6))
seq_ttd_plot <- seq_perc_reached + seq_ttd +
  plot_layout(nrow = 2, heights = c(1, 6))
seq_inf_plot <- seq_perc_reached + seq_infs +
  plot_layout(nrow = 2, heights = c(1, 6))
seq_prev_plot <- seq_perc_reached + seq_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 6))
seq_plot <- cowplot::plot_grid(seq_ttd_plot, seq_cuminf_plot, seq_inf_plot, seq_prev_plot, nrow = 1)

#########################################################################
#####           Shedding Frequency Sensitivity Analysis             #####
#########################################################################

# Generating values to vary shedding frequency over and dataframe to store outputs from model running
shed_sens <- seq(0.5, 2.5, 0.05)
shed_output <- data.frame(shed_total = rep(shed_sens, each = fixed_params$stochastic_sim), 
                          stochastic_realisation = 1:fixed_params$stochastic_sim, 
                          num_reads = fixed_params$num_reads, 
                          time_to_detection = NA,
                          cumulative_incidence = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(shed_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir_dust$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB, 
        
        # Sequencing Parameters
        shedding_freq = shed_sens[i], virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      time_to_detection <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      shed_output$time_to_detection[counter] <- time_to_detection
      shed_output$cumulative_incidence[counter] <- if(!is.na(time_to_detection)) sum(output2$n_SE_Output[1:(time_to_detection/fixed_params$dt)])/fixed_params$population_size else NA
      
      counter <- counter + 1   
    }
  }
  saveRDS(shed_output, file = "outputs/agg_shedFreq_sensitivity_analysis.rds")
} else {
  shed_output <- readRDS("outputs/agg_shedFreq_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
shed_df_summary <- shed_output %>%
  group_by(shed_total) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections, na.rm = TRUE),
            lower_inf = quantile(flight_infections, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.95),
            lower_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection)))/fixed_params$stochastic_sim) %>%
  ungroup() %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         lower_inf = ifelse(is.infinite(lower_inf), NaN, lower_inf),
         upper_inf = ifelse(is.infinite(upper_inf), NaN, upper_inf),
         lower_flight_prev = ifelse(is.infinite(lower_flight_prev), NaN, lower_flight_prev),
         upper_flight_prev = ifelse(is.infinite(upper_flight_prev), NaN, upper_flight_prev),
         lower_comm_prev = ifelse(is.infinite(lower_comm_prev), NaN, lower_comm_prev),
         upper_comm_prev = ifelse(is.infinite(upper_comm_prev), NaN, upper_comm_prev))

shed_cuminf <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_cuminf), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Frequency", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_ttd <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_ttd), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Frequency", y = "Time to Detection (Days)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_infs <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_inf), col = colours[1]) +
  geom_ribbon(aes(x = shed_total, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Shedding Frequency", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_flight_prev <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_flight_prev), col = colours[1]) +
  geom_ribbon(aes(x = shed_total, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Shedding Frequency", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_perc_reached <- ggplot(data = shed_df_summary) +
  geom_bar(aes(x = shed_total, y = perc_reached), col = "#4A4844",
           fill = adjustcolor(colours[3], alpha.f = 0.5), stat = "identity") +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  lims(y = c(0, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

shed_cuminf_plot <- shed_perc_reached + shed_cuminf +
  plot_layout(nrow = 2, heights = c(1, 6))
shed_ttd_plot <- shed_perc_reached + shed_ttd +
  plot_layout(nrow = 2, heights = c(1, 6))
shed_inf_plot <- shed_perc_reached + shed_infs +
  plot_layout(nrow = 2, heights = c(1, 6))
shed_prev_plot <- shed_perc_reached + shed_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 6))
shed_plot <- cowplot::plot_grid(shed_ttd_plot, shed_cuminf_plot, shed_inf_plot, shed_prev_plot, nrow = 1)

#########################################################################
#####             Flight Related Sensitivity Analysis               #####
#########################################################################

# Generating values to vary shedding frequency over and dataframe to store outputs from model running
proportion_AB <- 0.1
num_flights_sens <- seq(10, 300, 1/proportion_AB)
prop_flying <- num_flights_sens * fixed_params$capacity_per_flight / fixed_params$population_size
prop_flyingAB <- prop_flying * proportion_AB
flight_vars_df <- data.frame(num_flights = num_flights_sens, 
                             proportion_AB = proportion_AB, 
                             prop_flying = 100 * prop_flying, 
                             prop_flyingAB = prop_flyingAB * 100) 
flight_output <- data.frame(num_flights = rep(num_flights_sens, each = fixed_params$stochastic_sim), 
                            stochastic_realisation = 1:fixed_params$stochastic_sim, 
                            num_reads = fixed_params$num_reads, 
                            time_to_detection = NA,
                            cumulative_incidence = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(num_flights_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir_dust$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, num_flights = num_flights_sens[i], num_flightsAB = round(num_flights_sens[i] * proportion_AB), 
        
        # Sequencing Parameters
        shedding_freq = fixed_params$shedding_freq, virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      time_to_detection <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      flight_output$time_to_detection[counter] <- time_to_detection
      flight_output$cumulative_incidence[counter] <- if(!is.na(time_to_detection)) sum(output2$n_SE_Output[1:(time_to_detection/fixed_params$dt)])/fixed_params$population_size else NA
      
      counter <- counter + 1   
    }
  }
  saveRDS(flight_output, file = "outputs/agg_flights_sensitivity_analysis.rds")
} else {
  flight_output <- readRDS("outputs/agg_flights_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
flight_df_summary <- flight_output %>%
  left_join(flight_vars_df, by = "num_flights") %>%
  group_by(num_flights, proportion_AB, prop_flying, prop_flyingAB) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections, na.rm = TRUE),
            lower_inf = quantile(flight_infections, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.95),
            lower_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection)))/fixed_params$stochastic_sim) %>%
  ungroup() %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         lower_inf = ifelse(is.infinite(lower_inf), NaN, lower_inf),
         upper_inf = ifelse(is.infinite(upper_inf), NaN, upper_inf),
         lower_flight_prev = ifelse(is.infinite(lower_flight_prev), NaN, lower_flight_prev),
         upper_flight_prev = ifelse(is.infinite(upper_flight_prev), NaN, upper_flight_prev),
         lower_comm_prev = ifelse(is.infinite(lower_comm_prev), NaN, lower_comm_prev),
         upper_comm_prev = ifelse(is.infinite(upper_comm_prev), NaN, upper_comm_prev))

flight_cuminf <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = prop_flyingAB , y = avg_cuminf), col = colours[4]) +
  geom_ribbon(aes(x = prop_flyingAB , y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flying A->B", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_ttd <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = prop_flyingAB , y = avg_ttd), col = colours[4]) +
  geom_ribbon(aes(x = prop_flyingAB, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flying A->B", y = "Time to Detection (Days)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_infs <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = prop_flyingAB, y = avg_inf), col = colours[1]) +
  geom_ribbon(aes(x = prop_flyingAB, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - % Flying A->B", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_flight_prev <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = prop_flyingAB, y = avg_flight_prev), col = colours[1]) +
  geom_ribbon(aes(x = prop_flyingAB, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - % Flying A->B", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_perc_reached <- ggplot(data = flight_df_summary) +
  geom_bar(aes(x = prop_flyingAB , y = perc_reached), col = "#4A4844",
           fill = adjustcolor(colours[4], alpha.f = 0.5), stat = "identity") +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  lims(y = c(0, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

flight_cuminf_plot <- flight_perc_reached + flight_cuminf +
  plot_layout(nrow = 2, heights = c(1, 6))
flight_ttd_plot <- flight_perc_reached + flight_ttd +
  plot_layout(nrow = 2, heights = c(1, 6))
flight_inf_plot <- flight_perc_reached + flight_infs +
  plot_layout(nrow = 2, heights = c(1, 6))
flight_prev_plot <- flight_perc_reached + flight_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 6))
flight_plot <- cowplot::plot_grid(flight_ttd_plot, flight_cuminf_plot, flight_inf_plot, flight_prev_plot, nrow = 1)

#########################################################################
#####        Virus Shed Amount Related Sensitivity Analysis         #####
#########################################################################

# Generating values to vary shedding amount over and dataframe to store outputs from model running
ratio_sens <- lseq(1e-07, 1e-04, 40)
ratio_output <- data.frame(viral_ratio = rep(ratio_sens, each = fixed_params$stochastic_sim), 
                           stochastic_realisation = 1:fixed_params$stochastic_sim, 
                           num_reads = fixed_params$num_reads, 
                           time_to_detection = NA,
                           cumulative_incidence = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(ratio_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir_dust$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB, 
        
        # Sequencing Parameters
        shedding_freq = fixed_params$shedding_freq, virus_shed = fixed_params$non_virus_shed * ratio_sens[i], 
        non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      time_to_detection <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      ratio_output$time_to_detection[counter] <- time_to_detection
      ratio_output$cumulative_incidence[counter] <- if(!is.na(time_to_detection)) sum(output2$n_SE_Output[1:(time_to_detection/fixed_params$dt)])/fixed_params$population_size else NA
      
      counter <- counter + 1   
    }
  }
  saveRDS(ratio_output, file = "outputs/agg_viralRatio_sensitivity_analysis.rds")
} else{
  ratio_output <- readRDS("outputs/agg_viralRatio_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
ratio_df_summary <- ratio_output %>%
  group_by(viral_ratio) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections, na.rm = TRUE),
            lower_inf = quantile(flight_infections, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence, na.rm = TRUE, 0.95),
            lower_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection)))/fixed_params$stochastic_sim) %>%
  ungroup() %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         lower_inf = ifelse(is.infinite(lower_inf), NaN, lower_inf),
         upper_inf = ifelse(is.infinite(upper_inf), NaN, upper_inf),
         lower_flight_prev = ifelse(is.infinite(lower_flight_prev), NaN, lower_flight_prev),
         upper_flight_prev = ifelse(is.infinite(upper_flight_prev), NaN, upper_flight_prev),
         lower_comm_prev = ifelse(is.infinite(lower_comm_prev), NaN, lower_comm_prev),
         upper_comm_prev = ifelse(is.infinite(upper_comm_prev), NaN, upper_comm_prev))

ratio_cuminf <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_cuminf), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Viral Ratio", y = "Cumulative Incidence (%)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_ttd <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_ttd), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Viral Ratio", y = "Time to Detection (Days)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_infs <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_inf), col = colours[1]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Viral Ratio", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_ratio_prev <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_ratio_prev), col = colours[1]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_ratio_prev, ymin = lower_ratio_prev, ymax = upper_ratio_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - Viral Ratio", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_perc_reached <- ggplot(data = ratio_df_summary) +
  geom_bar(aes(x = viral_ratio, y = perc_reached), col = "#4A4844",
           fill = adjustcolor(colours[5], alpha.f = 0.5), stat = "identity") +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  lims(y = c(0, 100)) +
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

ratio_cuminf_plot <- ratio_perc_reached + ratio_cuminf +
  plot_layout(nrow = 2, heights = c(1, 6))
ratio_ttd_plot <- ratio_perc_reached + ratio_ttd +
  plot_layout(nrow = 2, heights = c(1, 6))
ratio_inf_plot <- ratio_perc_reached + ratio_infs +
  plot_layout(nrow = 2, heights = c(1, 6))
ratio_prev_plot <- ratio_perc_reached + ratio_ratio_prev +
  plot_layout(nrow = 2, heights = c(1, 6))
ratio_plot <- cowplot::plot_grid(ratio_ttd_plot, ratio_cuminf_plot, ratio_inf_plot, ratio_prev_plot, nrow = 1)

# Summarising and Plotting All the Plots Together
beta_plot + seq_plot + shed_plot + ratio_plot + flight_plot +
  plot_layout(nrow = 2)


#####

#### WHY ARE WE GETTING A BUNCH OF ZEROS HERE - WHAT'S UP WITH THAT???
#### IS THE BINOMIAL INTRODUCING TOO MUCH STOCHASTICITY???
#### only because the amount of viral shed was horrifically low i.e. 20
#### so binomial is returning a tonne of zeroes - when it occasionally has something, it's
#### a 1, which is translated into 1000s of reads with the non-sensical param combo we have currently
#### deterministic version is far better behaved
# plot(daily_df$daily_reads_stoch, type = "l")
# lines(daily_df$daily_reads_det, col = "red")
# plot(daily_df$daily_flightAB_infections, type = "l")
# plot(daily_df$daily_flightAB_prevalence, type = "l")
# plot(daily_df$daily_community_prevalence, type = "l")
# 
# plot(mod_output$n_inf_flightABOut, type = "l")
# plot(mod_output$amount_virus_aggFlight, type = "l")
# plot(mod_output$sample_amount_virus_aggFlight_stoch, type = "l")
# lines(mod_output$sample_amount_virus_aggFlight_det, col = "red")
# 
# sum(mod_output$sample_amount_virus_aggFlight_stoch)
# sum(mod_output$sample_amount_virus_aggFlight_det)
# 
# plot(mod_output$seq_reads_virus_aggFlight_stoch, type = "l")
# lines(mod_output$seq_reads_virus_aggFlight_det, col = "red")