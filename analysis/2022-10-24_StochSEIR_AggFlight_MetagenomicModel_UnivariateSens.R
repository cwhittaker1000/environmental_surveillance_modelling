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
source("models/2022-10-19_StochSEIR_AggFlight_Metagenomic_Model.R")
options(dplyr.summarise.inform = FALSE)

# Specifying model parameters
stochastic_sim <- 200                     # number of stochastic simulations to run
fixed_params <- list(
  
  # Epi Params
  beta = 0.5,                             # transmission rate - beta/sigma is R0 here
  gamma = 0.25,                           # 1/gamma = duration of latent period
  sigma = 0.2,                            # 1/sigma = duration of infectious period
  population_size = 10^7,                 # size of population in Location A
  start_infections = 20,                  # starting number of infections
  
  # Flight Params
  num_flights_leaveA = 200,               # number of flights leaving A per day
  num_flights_AtoB = 10,                  # number of flights leaving A per day that fly to B
  num_flights_arriveB = 200,              # total number of flights (including those from A) that arrive at B per day
  capacity_per_flight = 200,              # capacity of a single flight
  
  # Shedding Params
  non_virus_shed = 2 * 10^11,             # total amount of uninteresting NA shed per shedding event 
  ratio_virus_to_non_virus = 2 * 10^-5,   # ratio of interesting to uninteresting NA shed per shedding event (i.e. how much pathogen as a proportion of all NA shed)
  shedding_prop = 0.5,                    # proportion of infected individuals who actually shed NA of interest
  shedding_freq = 1,                      # rate of defecation per flight
  
  # Sequencing Params
  sample_volume = 1,                      # volume of wastewater sampled (note just a multiplicative term that doesn't affect relative abundance currently)
  flight_ww_volume = 800,                 # volume of pooled flight wastewater (note just a multiplicative term that doesn't affect relative abundance currently)
  met_bias = 1,                           # metagenomic bias term
  seq_tot = 10^9,                         # sequenced reads per day;
  num_reads = 5,                          # number of reads for "successful" detection
  
  # Misc Params
  dt = 0.2,                               # model timestep (in days)
  end = 500,                              # time (in days) to run model for
  seed = rpois(stochastic_sim, 10^9),     # stochastic seed for each simulation
  stochastic_sim = stochastic_sim         # number of stochastic sims
)

# Infilling Flight Parameters
prop_pop_flying <- 0.004                      # proportion of population in Location A taking flight every day 
prop_pop_flying_to_AB <- 0.05                 # proportion of flight taking population travelling Location A -> Location B (location with NAO)
num_flights_leaveA <- floor((prop_pop_flying * fixed_params$population_size) / fixed_params$capacity_per_flight)
num_flights_AtoB <- round(num_flights_leaveA * prop_pop_flying_to_AB)
fixed_params$num_flights_leaveA <- num_flights_leaveA
fixed_params$num_flights_AtoB <- num_flights_AtoB

# Sense Checking Parameters
1 / (prop_pop_flying * prop_pop_flying_to_AB)                       # 1 in every X of Location A population taking flight to Location B every day
num_flights_AtoB/fixed_params$num_flights_arriveB                   # what proportion of flights arriving at Location B come from Location A
fixed_params$beta/fixed_params$sigma                                # R0
fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus # Amount of virus shed per event

## for calculating the implied daily growth rate from the parameters
a <- 1
b <- fixed_params$sigma + fixed_params$gamma
c <- -(fixed_params$sigma * fixed_params$gamma) * ((fixed_params$beta/fixed_params$sigma) - 1)
implied_growth_rate <- (-b + sqrt(b^2 - 4 * a * c))/(2*a)
implied_growth_rate_per_day <- 100 * (exp(implied_growth_rate) - 1) # percentage increase day on day

new_run <- TRUE
colours <- c("#88D18A", "#788BFF", "#5B5750", "#F45B69", "#F18F01")

#########################################################################
#####                    Beta Sensitivity Analysis                  #####
#########################################################################

# Generating values to vary beta over and dataframe to store outputs from model running
R0 <- seq(1.5, 3.5, 0.1)         
beta_sens <- R0 * fixed_params$sigma 
R0_df <- data.frame(beta = beta_sens, R0 = beta_sens/fixed_params$sigma)
beta_output <- data.frame(beta = rep(beta_sens, each = fixed_params$stochastic_sim), 
                          stochastic_realisation = 1:fixed_params$stochastic_sim, 
                          num_reads = fixed_params$num_reads, 
                          time_to_detection_det = NA,
                          flight_infections_det = NA,
                          flight_prevalence_det = NA,
                          community_prevalence_det = NA,
                          cumulative_incidence_det = NA,
                          time_to_detection_stoch = NA,
                          flight_infections_stoch = NA,
                          flight_prevalence_stoch = NA,
                          community_prevalence_stoch = NA,
                          cumulative_incidence_stoch = NA)

# Running model and generating outputs for range of beta values
if (new_run) {
  counter <- 1
  for (i in 1:length(beta_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      
      mod <- stoch_seir$new(
        
        # Epidemiological Parameters
        beta = beta_sens[i], 
        gamma = fixed_params$gamma, 
        sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, 
        start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, 
        num_flights_leaveA = fixed_params$num_flights_leaveA, 
        num_flights_AtoB = fixed_params$num_flights_AtoB, 
        num_flights_arriveB = fixed_params$num_flights_arriveB, 
        
        # Shedding Parameters
        non_virus_shed = fixed_params$non_virus_shed, 
        virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        shedding_prop = fixed_params$shedding_prop,
        shedding_freq = fixed_params$shedding_freq, 
        
        # Sequencing Parameters
        met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, 
        samp_frac_aggFlight = fixed_params$sample_volume/(fixed_params$flight_ww_volume * fixed_params$num_flights_arriveB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      daily_output <- agg_to_daily(output2, fixed_params, fixed_params$seq_tot)
      time_to_detection <- daily_output  %>%
        dplyr::summarise(time_to_detection_det = time2[min(which(daily_reads_det >= fixed_params$num_reads))], 
                         time_to_detection_stoch = time2[min(which(daily_reads_stoch >= fixed_params$num_reads))])
      other_metrics_det <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_det)
      other_metrics_stoch <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_stoch)
        
      # Summarising Results
      ttd_det_success <- ifelse(!is.na(time_to_detection$time_to_detection_det), TRUE, FALSE)
      beta_output$time_to_detection_det[counter] <- time_to_detection$time_to_detection_det
      beta_output$flight_infections_det[counter] <- if (ttd_det_success) other_metrics_det$daily_infections_AtoB else NA
      beta_output$flight_prevalence_det[counter] <- if (ttd_det_success) 100 * other_metrics_det$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      beta_output$community_prevalence_det[counter] <- if (ttd_det_success) other_metrics_det$daily_prevalence_infection else NA
      beta_output$cumulative_incidence_det[counter] <- if (ttd_det_success) other_metrics_det$cumulative_incidence else NA
      
      ttd_stoch_success <- ifelse(!is.na(time_to_detection$time_to_detection_stoch), TRUE, FALSE)
      beta_output$time_to_detection_stoch[counter] <- time_to_detection$time_to_detection_stoch
      beta_output$flight_infections_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_infections_AtoB else NA
      beta_output$flight_prevalence_stoch[counter] <- if (ttd_stoch_success) 100 * other_metrics_stoch$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      beta_output$community_prevalence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_prevalence_infection else NA
      beta_output$cumulative_incidence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$cumulative_incidence else NA

      counter <- counter + 1   
    }
  }
  beta_output <- list(beta_output = beta_output, fixed_params = fixed_params)
  saveRDS(beta_output, file = "outputs/agg_beta_sensitivity_analysis.rds")
}  else {
  beta_output <- readRDS("outputs/agg_beta_sensitivity_analysis.rds")
}

# Summarising and plotting the output from beta sensitivity analysis
beta_df_summary <- beta_output$beta_output %>%
  left_join(R0_df, by = "beta") %>%
  group_by(beta, R0) %>%
  summarise(avg_ttd = mean(time_to_detection_det, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence_det, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections_det, na.rm = TRUE),
            lower_inf = quantile(flight_infections_det, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections_det, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence_det, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.95),
            avg_comm_prev = mean(community_prevalence_det, na.rm = TRUE),
            lower_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection_det))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection_det)))/fixed_params$stochastic_sim) %>%
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

beta_comm_prev <- ggplot(data = beta_df_summary) +
  geom_line(aes(x = R0, y = avg_comm_prev), col = colours[1]) +
  geom_ribbon(aes(x = R0, y = avg_comm_prev, ymin = lower_comm_prev, ymax = upper_comm_prev), 
              alpha = 0.2, fill = colours[1]) +
  labs(x = "Parameter Value - R0", y = "Location A Prev (%) At ToD") +
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
  plot_layout(nrow = 2, heights = c(1, 4.5))
beta_ttd_plot <- beta_perc_reached + beta_ttd +
  plot_layout(nrow = 2, heights = c(1, 4.5))
beta_inf_plot <- beta_perc_reached + beta_infs +
  plot_layout(nrow = 2, heights = c(1, 4.5))
beta_flight_prev_plot <- beta_perc_reached + beta_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
beta_comm_prev_plot <- beta_perc_reached + beta_comm_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
# note: reason comm prev is lower than flight prev is because the random day that achieves success
#       first is the day that has a higher than representative prevalence on the flight
beta_plot <- cowplot::plot_grid(beta_ttd_plot, beta_cuminf_plot, beta_inf_plot, beta_flight_prev_plot, beta_comm_prev_plot, nrow = 2)
ggsave(filename = "figures/univariate/beta_univariate_sensitivity.pdf",
       plot = beta_plot,
       width = 10.5,
       height = 7)

#########################################################################
#####               Seq Total Sensitivity Analysis                  #####
#########################################################################

# Generating values to vary seq total over and dataframe to store outputs from model running
seq_sens <- round(lseq(10^8, 10^10, 40))
seq_output <- data.frame(seq_total = rep(seq_sens, each = fixed_params$stochastic_sim), 
                         stochastic_realisation = 1:fixed_params$stochastic_sim, 
                         num_reads = fixed_params$num_reads, 
                         time_to_detection_det = NA,
                         flight_infections_det = NA,
                         flight_prevalence_det = NA,
                         community_prevalence_det = NA,
                         cumulative_incidence_det = NA,
                         time_to_detection_stoch = NA,
                         flight_infections_stoch = NA,
                         flight_prevalence_stoch = NA,
                         community_prevalence_stoch = NA,
                         cumulative_incidence_stoch = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(seq_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, 
        gamma = fixed_params$gamma, 
        sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, 
        start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, 
        num_flights_leaveA = fixed_params$num_flights_leaveA, 
        num_flights_AtoB = fixed_params$num_flights_AtoB, 
        num_flights_arriveB = fixed_params$num_flights_arriveB, 
        
        # Shedding Parameters
        non_virus_shed = fixed_params$non_virus_shed, 
        virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        shedding_prop = fixed_params$shedding_prop,
        shedding_freq = fixed_params$shedding_freq, 
        
        # Sequencing Parameters
        met_bias = fixed_params$met_bias, 
        seq_tot = seq_sens[i], 
        samp_frac_aggFlight = fixed_params$sample_volume/(fixed_params$flight_ww_volume * fixed_params$num_flights_arriveB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      daily_output <- agg_to_daily(output2, fixed_params, seq_sens[i])
      time_to_detection <- daily_output  %>%
        dplyr::summarise(time_to_detection_det = time2[min(which(daily_reads_det >= fixed_params$num_reads))], 
                         time_to_detection_stoch = time2[min(which(daily_reads_stoch >= fixed_params$num_reads))])
      other_metrics_det <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_det)
      other_metrics_stoch <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_stoch)
      
      # Summarising Results
      ttd_det_success <- ifelse(!is.na(time_to_detection$time_to_detection_det), TRUE, FALSE)
      seq_output$time_to_detection_det[counter] <- time_to_detection$time_to_detection_det
      seq_output$flight_infections_det[counter] <- if (ttd_det_success) other_metrics_det$daily_infections_AtoB else NA
      seq_output$flight_prevalence_det[counter] <- if (ttd_det_success) 100 * other_metrics_det$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      seq_output$community_prevalence_det[counter] <- if (ttd_det_success) other_metrics_det$daily_prevalence_infection else NA
      seq_output$cumulative_incidence_det[counter] <- if (ttd_det_success) other_metrics_det$cumulative_incidence else NA
      
      ttd_stoch_success <- ifelse(!is.na(time_to_detection$time_to_detection_stoch), TRUE, FALSE)
      seq_output$time_to_detection_stoch[counter] <- time_to_detection$time_to_detection_stoch
      seq_output$flight_infections_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_infections_AtoB else NA
      seq_output$flight_prevalence_stoch[counter] <- if (ttd_stoch_success) 100 * other_metrics_stoch$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      seq_output$community_prevalence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_prevalence_infection else NA
      seq_output$cumulative_incidence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$cumulative_incidence else NA
      
      counter <- counter + 1   
    }
  }
  seq_output <- list(seq_output = seq_output, fixed_params = fixed_params)
  saveRDS(seq_output, file = "outputs/agg_seqTotal_sensitivity_analysis.rds")
} else {
  seq_output <- readRDS("outputs/agg_seqTotal_sensitivity_analysis.rds")
}

# Summarising and plotting the output from seq total sensitivity analysis
seq_df_summary <- seq_output$seq_output %>%
  group_by(seq_total) %>%
  summarise(avg_ttd = mean(time_to_detection_det, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence_det, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections_det, na.rm = TRUE),
            lower_inf = quantile(flight_infections_det, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections_det, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence_det, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.95),
            avg_comm_prev = mean(community_prevalence_det, na.rm = TRUE),
            lower_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection_det))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection_det)))/fixed_params$stochastic_sim) %>%
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
  geom_line(aes(x = seq_total, y = avg_inf), col = colours[2]) +
  geom_ribbon(aes(x = seq_total, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[2]) +
  labs(x = "Parameter Value - Seq Total", y = "# Infections At ToD") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_flight_prev <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_flight_prev), col = colours[2]) +
  geom_ribbon(aes(x = seq_total, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[2]) +
  labs(x = "Parameter Value - Seq Total", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

seq_comm_prev <- ggplot(data = seq_df_summary) +
  geom_line(aes(x = seq_total, y = avg_comm_prev), col = colours[2]) +
  geom_ribbon(aes(x = seq_total, y = avg_comm_prev, ymin = lower_comm_prev, ymax = upper_comm_prev), 
              alpha = 0.2, fill = colours[2]) +
  labs(x = "Parameter Value - Seq Total", y = "Location A Prev (%) At ToD") +
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
  plot_layout(nrow = 2, heights = c(1, 4.5))
seq_ttd_plot <- seq_perc_reached + seq_ttd +
  plot_layout(nrow = 2, heights = c(1, 4.5))
seq_inf_plot <- seq_perc_reached + seq_infs +
  plot_layout(nrow = 2, heights = c(1, 4.5))
seq_flight_prev_plot <- seq_perc_reached + seq_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
seq_comm_prev_plot <- seq_perc_reached + seq_comm_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
# note: reason comm prev is lower than flight prev is because the random day that achieves success
#       first is the day that has a higher than representative prevalence on the flight
seq_plot <- cowplot::plot_grid(seq_ttd_plot, seq_cuminf_plot, seq_inf_plot, seq_flight_prev_plot, seq_comm_prev_plot, nrow = 2)
ggsave(filename = "figures/univariate/seqTotal_univariate_sensitivty.pdf",
       plot = beta_plot,
       width = 10.5,
       height = 7)


#########################################################################
#####           Shedding Frequency Sensitivity Analysis             #####
#########################################################################

# Generating values to vary shedding frequency over and dataframe to store outputs from model running
shed_prop_sens <- seq(0.05, 1, 0.05)
shed_output <- data.frame(shed_total = rep(shed_prop_sens, each = fixed_params$stochastic_sim), 
                          stochastic_realisation = 1:fixed_params$stochastic_sim, 
                          num_reads = fixed_params$num_reads, 
                          time_to_detection_det = NA,
                          flight_infections_det = NA,
                          flight_prevalence_det = NA,
                          community_prevalence_det = NA,
                          cumulative_incidence_det = NA,
                          time_to_detection_stoch = NA,
                          flight_infections_stoch = NA,
                          flight_prevalence_stoch = NA,
                          community_prevalence_stoch = NA,
                          cumulative_incidence_stoch = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(shed_prop_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, 
        gamma = fixed_params$gamma, 
        sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, 
        start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, 
        num_flights_leaveA = fixed_params$num_flights_leaveA, 
        num_flights_AtoB = fixed_params$num_flights_AtoB, 
        num_flights_arriveB = fixed_params$num_flights_arriveB, 
        
        # Shedding Parameters
        non_virus_shed = fixed_params$non_virus_shed, 
        virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        shedding_prop = shed_prop_sens[i],
        shedding_freq = fixed_params$shedding_freq, 
        
        # Sequencing Parameters
        met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, 
        samp_frac_aggFlight = fixed_params$sample_volume/(fixed_params$flight_ww_volume * fixed_params$num_flights_arriveB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      daily_output <- agg_to_daily(output2, fixed_params, fixed_params$seq_tot)
      time_to_detection <- daily_output  %>%
        dplyr::summarise(time_to_detection_det = time2[min(which(daily_reads_det >= fixed_params$num_reads))], 
                         time_to_detection_stoch = time2[min(which(daily_reads_stoch >= fixed_params$num_reads))])
      other_metrics_det <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_det)
      other_metrics_stoch <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_stoch)
      
      # Summarising Results
      ttd_det_success <- ifelse(!is.na(time_to_detection$time_to_detection_det), TRUE, FALSE)
      shed_output$time_to_detection_det[counter] <- time_to_detection$time_to_detection_det
      shed_output$flight_infections_det[counter] <- if (ttd_det_success) other_metrics_det$daily_infections_AtoB else NA
      shed_output$flight_prevalence_det[counter] <- if (ttd_det_success) 100 * other_metrics_det$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      shed_output$community_prevalence_det[counter] <- if (ttd_det_success) other_metrics_det$daily_prevalence_infection else NA
      shed_output$cumulative_incidence_det[counter] <- if (ttd_det_success) other_metrics_det$cumulative_incidence else NA
      
      ttd_stoch_success <- ifelse(!is.na(time_to_detection$time_to_detection_stoch), TRUE, FALSE)
      shed_output$time_to_detection_stoch[counter] <- time_to_detection$time_to_detection_stoch
      shed_output$flight_infections_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_infections_AtoB else NA
      shed_output$flight_prevalence_stoch[counter] <- if (ttd_stoch_success) 100 * other_metrics_stoch$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      shed_output$community_prevalence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_prevalence_infection else NA
      shed_output$cumulative_incidence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$cumulative_incidence else NA
      
      counter <- counter + 1   
    }
  }
  shed_output <- list(shed_output = shed_output, fixed_params = fixed_params)
  saveRDS(shed_output, file = "outputs/agg_shedProp_sensitivity_analysis.rds")
} else {
  shed_output <- readRDS("outputs/agg_shedProp_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
shed_df_summary <- shed_output$shed_output %>%
  group_by(shed_total) %>%
  summarise(avg_ttd = mean(time_to_detection_det, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence_det, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections_det, na.rm = TRUE),
            lower_inf = quantile(flight_infections_det, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections_det, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence_det, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.95),
            avg_comm_prev = mean(community_prevalence_det, na.rm = TRUE),
            lower_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection_det))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection_det)))/fixed_params$stochastic_sim) %>%
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
  labs(x = "Parameter Value - Shedding Proportion", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_ttd <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_ttd), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Proportion", y = "Time to Detection (Days)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_infs <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_inf), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Proportion", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_flight_prev <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_flight_prev), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Proportion", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

shed_comm_prev <- ggplot(data = shed_df_summary) +
  geom_line(aes(x = shed_total, y = avg_comm_prev), col = colours[3]) +
  geom_ribbon(aes(x = shed_total, y = avg_comm_prev, ymin = lower_comm_prev, ymax = upper_comm_prev), 
              alpha = 0.2, fill = colours[3]) +
  labs(x = "Parameter Value - Shedding Proportion", y = "Location A Prev (%) At ToD") +
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
  plot_layout(nrow = 2, heights = c(1, 4.5))
shed_ttd_plot <- shed_perc_reached + shed_ttd +
  plot_layout(nrow = 2, heights = c(1, 4.5))
shed_inf_plot <- shed_perc_reached + shed_infs +
  plot_layout(nrow = 2, heights = c(1, 4.5))
shed_flight_prev_plot <- shed_perc_reached + shed_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
shed_comm_prev_plot <- shed_perc_reached + shed_comm_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
shed_plot <- cowplot::plot_grid(shed_ttd_plot, shed_cuminf_plot, shed_inf_plot, shed_flight_prev_plot, shed_comm_prev_plot, nrow = 2)
shed_plot
ggsave(filename = "figures/univariate/shedProp_univariate_sensitivty.pdf",
       plot = shed_plot,
       width = 10.5,
       height = 7)


#########################################################################
#####             Flight Related Sensitivity Analysis               #####
#########################################################################
prop_AB_overall_sens <- seq(0.01, 0.20, 0.01)
num_flights_AB_sens <- fixed_params$num_flights_leaveA * prop_AB_overall_sens
flight_output <- data.frame(num_flightsAB = rep(num_flights_AB_sens, each = fixed_params$stochastic_sim), 
                            stochastic_realisation = 1:fixed_params$stochastic_sim, 
                            num_reads = fixed_params$num_reads, 
                            time_to_detection_det = NA,
                            flight_infections_det = NA,
                            flight_prevalence_det = NA,
                            community_prevalence_det = NA,
                            cumulative_incidence_det = NA,
                            time_to_detection_stoch = NA,
                            flight_infections_stoch = NA,
                            flight_prevalence_stoch = NA,
                            community_prevalence_stoch = NA,
                            cumulative_incidence_stoch = NA)
flight_vars_df <- data.frame(num_flightsAB = num_flights_AB_sens,
                             perc_flightsAB = 100 * prop_AB_overall_sens,
                             perc_indivs_flyAB = 100 * prop_AB_overall_sens * prop_pop_flying,
                             perc_flights_BfromA = 100 * num_flights_AB_sens/fixed_params$num_flights_arriveB) 


# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(num_flights_AB_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, 
        gamma = fixed_params$gamma, 
        sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, 
        start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, 
        num_flights_leaveA = fixed_params$num_flights_leaveA, 
        num_flights_AtoB = num_flights_AB_sens[i],
        num_flights_arriveB = fixed_params$num_flights_arriveB, 
        
        # Shedding Parameters
        non_virus_shed = fixed_params$non_virus_shed, 
        virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
        shedding_prop = fixed_params$shedding_prop,
        shedding_freq = fixed_params$shedding_freq, 
        
        # Sequencing Parameters
        met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, 
        samp_frac_aggFlight = fixed_params$sample_volume/(fixed_params$flight_ww_volume * fixed_params$num_flights_arriveB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      daily_output <- agg_to_daily(output2, fixed_params, fixed_params$seq_tot)
      time_to_detection <- daily_output  %>%
        dplyr::summarise(time_to_detection_det = time2[min(which(daily_reads_det >= fixed_params$num_reads))], 
                         time_to_detection_stoch = time2[min(which(daily_reads_stoch >= fixed_params$num_reads))])
      other_metrics_det <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_det)
      other_metrics_stoch <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_stoch)
      
      # Summarising Results
      ttd_det_success <- ifelse(!is.na(time_to_detection$time_to_detection_det), TRUE, FALSE)
      flight_output$time_to_detection_det[counter] <- time_to_detection$time_to_detection_det
      flight_output$flight_infections_det[counter] <- if (ttd_det_success) other_metrics_det$daily_infections_AtoB else NA
      flight_output$flight_prevalence_det[counter] <- if (ttd_det_success) 100 * other_metrics_det$daily_infections_AtoB/(num_flights_AB_sens[i] * fixed_params$capacity_per_flight) else NA
      flight_output$community_prevalence_det[counter] <- if (ttd_det_success) other_metrics_det$daily_prevalence_infection else NA
      flight_output$cumulative_incidence_det[counter] <- if (ttd_det_success) other_metrics_det$cumulative_incidence else NA
      
      ttd_stoch_success <- ifelse(!is.na(time_to_detection$time_to_detection_stoch), TRUE, FALSE)
      flight_output$time_to_detection_stoch[counter] <- time_to_detection$time_to_detection_stoch
      flight_output$flight_infections_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_infections_AtoB else NA
      flight_output$flight_prevalence_stoch[counter] <- if (ttd_stoch_success) 100 * other_metrics_stoch$daily_infections_AtoB/(num_flights_AB_sens[i] * fixed_params$capacity_per_flight) else NA
      flight_output$community_prevalence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_prevalence_infection else NA
      flight_output$cumulative_incidence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$cumulative_incidence else NA
      
      counter <- counter + 1     
    }
  }
  flight_output <- list(flight_output = flight_output, fixed_params = fixed_params)
  saveRDS(flight_output, file = "outputs/agg_flightsAB_sensitivity_analysis.rds")
} else {
  flight_output <- readRDS("outputs/agg_flightsAB_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
flight_df_summary <- flight_output$flight_output %>%
  left_join(flight_vars_df, by = "num_flightsAB") %>%
  group_by(num_flightsAB, perc_flightsAB, perc_indivs_flyAB, perc_flights_BfromA) %>%
  summarise(avg_ttd = mean(time_to_detection_det, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence_det, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections_det, na.rm = TRUE),
            lower_inf = quantile(flight_infections_det, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections_det, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence_det, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.95),
            avg_comm_prev = mean(community_prevalence_det, na.rm = TRUE),
            lower_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection_det))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection_det)))/fixed_params$stochastic_sim) %>%
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
  geom_line(aes(x = perc_flights_BfromA , y =  avg_cuminf), col = colours[4]) +
  geom_ribbon(aes(x = perc_flights_BfromA , y = avg_cuminf , ymin = lower_cuminf, ymax = upper_cuminf), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flights Arriving at B From A", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_ttd <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = perc_flights_BfromA , y = avg_ttd), col = colours[4]) +
  geom_ribbon(aes(x = perc_flights_BfromA, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flights Arriving at B From A", y = "Time to Detection (Days)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_infs <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = perc_flights_BfromA, y = avg_inf), col = colours[4]) +
  geom_ribbon(aes(x = perc_flights_BfromA, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flights Arriving at B From A", y = "# Infections At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_flight_prev <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = perc_flights_BfromA, y = avg_flight_prev), col = colours[4]) +
  geom_ribbon(aes(x = perc_flights_BfromA, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - % Flights Arriving at B From A", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_comm_prev <- ggplot(data = flight_df_summary) +
  geom_line(aes(x = perc_flights_BfromA, y = avg_comm_prev), col = colours[4]) +
  geom_ribbon(aes(x = perc_flights_BfromA, y = avg_comm_prev, ymin = lower_comm_prev, ymax = upper_comm_prev), 
              alpha = 0.2, fill = colours[4]) +
  labs(x = "Parameter Value - Flights Arriving at B From A", y = "Location A Prev (%) At ToD") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

flight_perc_reached <- ggplot(data = flight_df_summary) +
  geom_bar(aes(x = perc_flights_BfromA , y = perc_reached), col = "#4A4844",
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
  plot_layout(nrow = 2, heights = c(1, 4.5))
flight_ttd_plot <- flight_perc_reached + flight_ttd +
  plot_layout(nrow = 2, heights = c(1, 4.5))
flight_inf_plot <- flight_perc_reached + flight_infs +
  plot_layout(nrow = 2, heights = c(1, 4.5))
flight_flight_prev_plot <- flight_perc_reached + flight_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
flight_comm_prev_plot <- flight_perc_reached + flight_comm_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
flight_plot <- cowplot::plot_grid(flight_ttd_plot, flight_cuminf_plot, flight_inf_plot, flight_flight_prev_plot, flight_comm_prev_plot, nrow = 2)
flight_plot
ggsave(filename = "figures/univariate/flightsAB_univariate_sensitivty.pdf",
       plot = flight_plot,
       width = 10.5,
       height = 7)

#########################################################################
#####        Virus Shed Amount Related Sensitivity Analysis         #####
#########################################################################

# Generating values to vary shedding amount over and dataframe to store outputs from model running
ratio_sens <- lseq(10^-7, 10^-4, 40)
ratio_output <- data.frame(viral_ratio = rep(ratio_sens, each = fixed_params$stochastic_sim), 
                           stochastic_realisation = 1:fixed_params$stochastic_sim, 
                           num_reads = fixed_params$num_reads, 
                           time_to_detection_det = NA,
                           flight_infections_det = NA,
                           flight_prevalence_det = NA,
                           community_prevalence_det = NA,
                           cumulative_incidence_det = NA,
                           time_to_detection_stoch = NA,
                           flight_infections_stoch = NA,
                           flight_prevalence_stoch = NA,
                           community_prevalence_stoch = NA,
                           cumulative_incidence_stoch = NA)

# Running model and generating outputs for range of seq total values
if (new_run) {
  counter <- 1
  for (i in 1:length(ratio_sens)) {
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Running the Model
      set.seed(fixed_params$seed[j])
      mod <- stoch_seir$new(
        
        # Epidemiological Parameters
        beta = fixed_params$beta, 
        gamma = fixed_params$gamma, 
        sigma = fixed_params$sigma, 
        population_size = fixed_params$population_size, 
        start_infections = fixed_params$start_infections,
        
        # Flight Parameters
        capacity_per_flight = fixed_params$capacity_per_flight, 
        num_flights_leaveA = fixed_params$num_flights_leaveA, 
        num_flights_AtoB = fixed_params$num_flights_AtoB,
        num_flights_arriveB = fixed_params$num_flights_arriveB, 
        
        # Shedding Parameters
        non_virus_shed = fixed_params$non_virus_shed, 
        virus_shed = fixed_params$non_virus_shed * ratio_sens[i], 
        shedding_prop = fixed_params$shedding_prop,
        shedding_freq = fixed_params$shedding_freq, 
        
        # Sequencing Parameters
        met_bias = fixed_params$met_bias, 
        seq_tot = fixed_params$seq_tot, 
        samp_frac_aggFlight = fixed_params$sample_volume/(fixed_params$flight_ww_volume * fixed_params$num_flights_arriveB),
        
        # Miscellaenous Parameters
        dt = fixed_params$dt)
      
      # Extracting Output
      output <- mod$run(1:(fixed_params$end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Calculating Time to Detection
      daily_output <- agg_to_daily(output2, fixed_params, fixed_params$seq_tot)
      time_to_detection <- daily_output  %>%
        dplyr::summarise(time_to_detection_det = time2[min(which(daily_reads_det >= fixed_params$num_reads))], 
                         time_to_detection_stoch = time2[min(which(daily_reads_stoch >= fixed_params$num_reads))])
      other_metrics_det <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_det)
      other_metrics_stoch <- daily_output %>%
        filter(time2 == time_to_detection$time_to_detection_stoch)
      
      # Summarising Results
      ttd_det_success <- ifelse(!is.na(time_to_detection$time_to_detection_det), TRUE, FALSE)
      ratio_output$time_to_detection_det[counter] <- time_to_detection$time_to_detection_det
      ratio_output$flight_infections_det[counter] <- if (ttd_det_success) other_metrics_det$daily_infections_AtoB else NA
      ratio_output$flight_prevalence_det[counter] <- if (ttd_det_success) 100 * other_metrics_det$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      ratio_output$community_prevalence_det[counter] <- if (ttd_det_success) other_metrics_det$daily_prevalence_infection else NA
      ratio_output$cumulative_incidence_det[counter] <- if (ttd_det_success) other_metrics_det$cumulative_incidence else NA
      
      ttd_stoch_success <- ifelse(!is.na(time_to_detection$time_to_detection_stoch), TRUE, FALSE)
      ratio_output$time_to_detection_stoch[counter] <- time_to_detection$time_to_detection_stoch
      ratio_output$flight_infections_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_infections_AtoB else NA
      ratio_output$flight_prevalence_stoch[counter] <- if (ttd_stoch_success) 100 * other_metrics_stoch$daily_infections_AtoB/(fixed_params$num_flights_AtoB * fixed_params$capacity_per_flight) else NA
      ratio_output$community_prevalence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$daily_prevalence_infection else NA
      ratio_output$cumulative_incidence_stoch[counter] <- if (ttd_stoch_success) other_metrics_stoch$cumulative_incidence else NA
      
      counter <- counter + 1   
    }
  }
  ratio_output <- list(ratio_output = ratio_output, fixed_params = fixed_params)
  saveRDS(ratio_output, file = "outputs/agg_viralRatio_sensitivity_analysis.rds")
} else{
  ratio_output <- readRDS("outputs/agg_viralRatio_sensitivity_analysis.rds")
}

# Summarising and plotting the output from shedding sensitivity analysis
ratio_df_summary <- ratio_output$ratio_output %>%
  group_by(viral_ratio) %>%
  summarise(avg_ttd = mean(time_to_detection_det, na.rm = TRUE),
            lower_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.05),
            upper_ttd = quantile(time_to_detection_det, na.rm = TRUE, 0.95),
            avg_cuminf = mean(cumulative_incidence_det, na.rm = TRUE),
            lower_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.05),
            upper_cuminf = quantile(cumulative_incidence_det, na.rm = TRUE, 0.95),
            avg_inf = mean(flight_infections_det, na.rm = TRUE),
            lower_inf = quantile(flight_infections_det, na.rm = TRUE, 0.05),
            upper_inf = quantile(flight_infections_det, na.rm = TRUE, 0.95),
            avg_flight_prev = mean(flight_prevalence_det, na.rm = TRUE),
            lower_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.05),
            upper_flight_prev = quantile(flight_prevalence_det, na.rm = TRUE, 0.95),
            avg_comm_prev = mean(community_prevalence_det, na.rm = TRUE),
            lower_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.05),
            upper_comm_prev = quantile(community_prevalence_det, na.rm = TRUE, 0.95),
            num_reached = paste0(fixed_params$stochastic_sim - sum(is.na(time_to_detection_det))),
            perc_reached = 100 * (fixed_params$stochastic_sim - sum(is.na(time_to_detection_det)))/fixed_params$stochastic_sim) %>%
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
  labs(x = "Parameter Value - Ratio of Virus:Non-Virus NA", y = "Cumulative Incidence (%)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_ttd <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_ttd), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Ratio of Virus:Non-Virus NA", y = "Time to Detection (Days)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_infs <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_inf), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_inf, ymin = lower_inf, ymax = upper_inf), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Ratio of Virus:Non-Virus NA", y = "# Infections At ToD") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_flight_prev <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_flight_prev), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_flight_prev, ymin = lower_flight_prev, ymax = upper_flight_prev), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Ratio of Virus:Non-Virus NA", y = "Flight Prev (%) At ToD") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

ratio_comm_prev <- ggplot(data = ratio_df_summary) +
  geom_line(aes(x = viral_ratio, y = avg_comm_prev), col = colours[5]) +
  geom_ribbon(aes(x = viral_ratio, y = avg_comm_prev, ymin = lower_comm_prev, ymax = upper_comm_prev), 
              alpha = 0.2, fill = colours[5]) +
  labs(x = "Parameter Value - Ratio of Virus:Non-Virus NA", y = "Location A Prev (%) At ToD") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
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
  plot_layout(nrow = 2, heights = c(1, 4.5))
ratio_ttd_plot <- ratio_perc_reached + ratio_ttd +
  plot_layout(nrow = 2, heights = c(1, 4.5))
ratio_inf_plot <- ratio_perc_reached + ratio_infs +
  plot_layout(nrow = 2, heights = c(1, 4.5))
ratio_flight_prev_plot <- ratio_perc_reached + ratio_flight_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
ratio_comm_prev_plot <- ratio_perc_reached + ratio_comm_prev +
  plot_layout(nrow = 2, heights = c(1, 4.5))
ratio_plot <- cowplot::plot_grid(ratio_ttd_plot, ratio_cuminf_plot, ratio_inf_plot, ratio_flight_prev_plot, ratio_comm_prev_plot, nrow = 2)
ratio_plot
ggsave(filename = "figures/univariate/virusRatio_univariate_sensitivty.pdf",
       plot = ratio_plot,
       width = 10.5,
       height = 7)

# Summarising and Plotting All the Plots Together
beta_plot <- cowplot::plot_grid(beta_ttd_plot, beta_cuminf_plot, beta_inf_plot, beta_flight_prev_plot, beta_comm_prev_plot, nrow = 1)
seq_plot <- cowplot::plot_grid(seq_ttd_plot, seq_cuminf_plot, seq_inf_plot, seq_flight_prev_plot, seq_comm_prev_plot, nrow = 1)
shed_plot <- cowplot::plot_grid(shed_ttd_plot, shed_cuminf_plot, shed_inf_plot, shed_flight_prev_plot, shed_comm_prev_plot, nrow = 1)
flight_plot <- cowplot::plot_grid(flight_ttd_plot, flight_cuminf_plot, flight_inf_plot, flight_flight_prev_plot, flight_comm_prev_plot, nrow = 1)
ratio_plot <- cowplot::plot_grid(ratio_ttd_plot, ratio_cuminf_plot, ratio_inf_plot, ratio_flight_prev_plot, ratio_comm_prev_plot, nrow = 1)

all_plots <- beta_plot + seq_plot + shed_plot + ratio_plot + flight_plot +
  plot_layout(nrow = 5)
ggsave(filename = "figures/univariate/allPlots_univariate_sensitivty.pdf",
       plot = all_plots,
       width = 17.5,
       height = 17.5)

## Code below illustrates some weirdness that as prop_flightsAB increases, so does time to infection
## crucial insight here is that as we hold num_flights constant and increase num_flights AB,
# we increase overall number of people who are infected flying to location B, but the proportion of all
# travellers they represent stays the same (e.g. double number of flights A->B = double number of infections, but also
# double number of uninfected passengers). So within the context of the current model, we infer an increasing relationship.
# In reality, having more people contributing to the system is probably useful, but within the current framework,
# that's not being shown. 

# set.seed(fixed_params$seed[j])
# mod1 <- stoch_seir_dust$new(
#   
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
#   
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, 
#   num_flights = fixed_params$num_flights, num_flightsAB = 50, 
#   
#   # Sequencing Parameters
#   shedding_prop = fixed_params$shedding_prop, shedding_freq = fixed_params$shedding_freq, 
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
#   seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output1 <- mod1$run(1:(fixed_params$end/fixed_params$dt))
# output1 <- mod1$transform_variables(output1)
# ttd_fun(mod1, output1, fixed_params$num_reads)
# 
# set.seed(fixed_params$seed[j])
# mod2 <- stoch_seir_dust$new(
#   
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
#   
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, 
#   num_flights = fixed_params$num_flights, num_flightsAB = 100, 
#   
#   # Sequencing Parameters
#   shedding_prop = fixed_params$shedding_prop, shedding_freq = fixed_params$shedding_freq, 
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
#   seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output2 <- mod2$run(1:(fixed_params$end/fixed_params$dt))
# output2 <- mod2$transform_variables(output2)
# ttd_fun(mod2, output2, fixed_params$num_reads)
# 
# plot(output2$n_inf_flightABOut, type = "l", col = "red")
# lines(output1$n_inf_flightABOut, type = "l")
# 
# plot(output2$infected_indiv_shedding_events_Out, type = "l", col = "red")
# lines(output1$infected_indiv_shedding_events_Out, type = "l")
# 
# plot(output2$uninfected_indiv_shedding_events_Out, type = "l", col = "red")
# lines(output1$uninfected_indiv_shedding_events_Out, type = "l")

#### Important

## Code below illustrates some weirdness when shedding freq goes really low, and actually
## we end up predicting earlier detection than when shedding freq is really high.
## Explanation is that at really low levels of shedding (which gets applied to both non viral and viral material)
## you get a lot more variation in the ratio of infected to uninfected shedding events

## Below key insight comes from:
# plot(output2$infected_indiv_shedding_events_Out/output2$uninfected_indiv_shedding_events_Out, col = "red")
# lines(output1$infected_indiv_shedding_events_Out/output1$uninfected_indiv_shedding_events_Out)

# output2 which has higher shedding frequency has a far less variable proportion over time (smoother curve);
# output1 is a lot more variable - often it's 0, but sometimes, the proportion of infected:uninfected is a lot higher than 
# the mean (because lambda of poisson is small, so lots of variation relative to the mean); so you get higher ratio of infected to
# uninfected shedding events, a lot more viral NA (in relative terms) in the wastewater, and as a result, lots more 
# sequencing reads (for that specific)

# set.seed(10)
# mod1 <- stoch_seir_dust$new(
# 
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma,
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
# 
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB,
# 
#   # Sequencing Parameters
#   shedding_freq = 0.05, shedding_prop = 1,
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus,
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias,
#   seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
# 
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output1 <- mod1$run(1:(fixed_params$end/fixed_params$dt))
# output1 <- mod1$transform_variables(output1)
# ttd_fun(mod1, output1, fixed_params$num_reads)
# 
# df1 <- data.frame(time = output1$time, 
#                  reads_det = output1$seq_reads_virus_aggFlight_det_Out,
#                  reads_stoch = output1$seq_reads_virus_aggFlight_stoch_Out,
#                  flightAB_infections = output1$n_inf_flightABOut,
#                  community_prevalence = output1$I)
# daily_df1 <- df1 %>%
#   dplyr::mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   dplyr::mutate(time3 = midpoints(time2)) %>%
#   group_by(time2, time3) %>%
#   summarise(daily_reads_det = sum(reads_det),
#             daily_reads_stoch = sum(reads_stoch),
#             daily_flightAB_infections = sum(flightAB_infections),
#             daily_flightAB_prevalence = 100 * daily_flightAB_infections/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB),
#             daily_community_prevalence = 100 * mean(community_prevalence)/fixed_params$population_size) %>%
#   ungroup(time2)
# 
# set.seed(10)
# mod2 <- stoch_seir_dust$new(
# 
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma,
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
# 
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, num_flights = fixed_params$num_flights, num_flightsAB = fixed_params$num_flightsAB,
# 
#   # Sequencing Parameters
#   shedding_freq = 10, shedding_prop = 1,
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus,
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias,
#   seq_tot = fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
# 
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output2 <- mod2$run(1:(fixed_params$end/fixed_params$dt))
# output2 <- mod2$transform_variables(output2)
# ttd_fun(mod2, output2, fixed_params$num_reads)
# df2 <- data.frame(time = output2$time, 
#                   reads_det = output2$seq_reads_virus_aggFlight_det_Out,
#                   reads_stoch = output2$seq_reads_virus_aggFlight_stoch_Out,
#                   flightAB_infections = output2$n_inf_flightABOut,
#                   community_prevalence = output2$I)
# daily_df2 <- df2 %>%
#   dplyr::mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   dplyr::mutate(time3 = midpoints(time2)) %>%
#   group_by(time2, time3) %>%
#   summarise(daily_reads_det = sum(reads_det),
#             daily_reads_stoch = sum(reads_stoch),
#             daily_flightAB_infections = sum(flightAB_infections),
#             daily_flightAB_prevalence = 100 * daily_flightAB_infections/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB),
#             daily_community_prevalence = 100 * mean(community_prevalence)/fixed_params$population_size) %>%
#   ungroup(time2)
# 
# plot(daily_df1$daily_reads_det, ylim = c(0, 6))
# lines(daily_df2$daily_reads_det)
# 
# par(mfrow = c(2, 2))
# plot(output1$n_inf_flightABOut, type = "l")
# points(output2$n_inf_flightABOut, col = "red")
# 
# plot(output1$sample_amount_virus_aggFlight_det_Out * 10/0.05, type = "l")
# lines(output2$sample_amount_virus_aggFlight_det_Out, col = "red")
# 
# plot(output2$infected_indiv_shedding_events_Out, col = "red")
# lines(output1$infected_indiv_shedding_events_Out * 10/0.05, type = "l")
# 
# plot(output2$infected_indiv_shedding_events_Out/output2$uninfected_indiv_shedding_events_Out, col = "red")
# lines(output1$infected_indiv_shedding_events_Out/output1$uninfected_indiv_shedding_events_Out)
# 
# ### Runing and Exploring WHY We See the Results We Do Whilst Varying PropAB
# 
# # Running Model W/ Low Prop
# j <- 1000
# set.seed(j)
# mod1 <- stoch_seir_dust$new(
#   
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
#   
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, 
#   num_flights = 500, 
#   num_flightsAB = 1, 
#   
#   # Sequencing Parameters
#   shedding_prop = fixed_params$shedding_prop, shedding_freq = fixed_params$shedding_freq, 
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
#   seq_tot = 10 * fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output1 <- mod1$run(1:(fixed_params$end/fixed_params$dt))
# output1 <- mod1$transform_variables(output1)
# df1 <- data.frame(time = output1$time, 
#                   reads_det = output1$seq_reads_virus_aggFlight_det_Out,
#                   flightAB_infections = output1$n_inf_flightABOut)
# daily_df1 <- df1 %>%
#   dplyr::mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   dplyr::mutate(time3 = midpoints(time2)) %>%
#   group_by(time2, time3) %>%
#   summarise(daily_reads_det = sum(reads_det),
#             daily_flightAB_infections = sum(flightAB_infections)) %>%
#   ungroup(time2)
# 
# # Running Model W/ High Prop
# set.seed(j)
# mod2 <- stoch_seir_dust$new(
#   
#   # Epidemiological Parameters
#   beta = fixed_params$beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
#   population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
#   
#   # Flight Parameters
#   capacity_per_flight = fixed_params$capacity_per_flight, 
#   num_flights = 500, 
#   num_flightsAB = 10, 
#   
#   # Sequencing Parameters
#   shedding_prop = fixed_params$shedding_prop, shedding_freq = fixed_params$shedding_freq, 
#   virus_shed = fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus, 
#   non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
#   seq_tot = 10 * fixed_params$seq_tot, samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * fixed_params$num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = fixed_params$dt)
# output2 <- mod2$run(1:(fixed_params$end/fixed_params$dt))
# output2 <- mod2$transform_variables(output2)
# df2 <- data.frame(time = output2$time, 
#                   reads_det = output2$seq_reads_virus_aggFlight_det_Out,
#                   flightAB_infections = output2$n_inf_flightABOut)
# daily_df2 <- df2 %>%
#   dplyr::mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   dplyr::mutate(time3 = midpoints(time2)) %>%
#   group_by(time2, time3) %>%
#   summarise(daily_reads_det = sum(reads_det),
#             daily_flightAB_infections = sum(flightAB_infections)) %>%
#   ungroup(time2)
# 
# # Checking TTD - see that TTD for 5 reads is shorter for mod1 (lower prop, unexpected) but time
# # to 5 infections is shorter for mod2 (as expected)
# ttd_fun(mod1, output1, fixed_params$num_reads) # 1  AB flight
# ttd_fun(mod2, output2, fixed_params$num_reads) # 25 AB flights
# 
# par(mfrow = c(4, 2), oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 1))
# plot(daily_df1$daily_flightAB_infections, type = "l", xlab = "time (days)", ylab = "# infections A->B")
# lines(daily_df2$daily_flightAB_infections, type = "l", col = "red")
# legend(160, 20, legend=c("High Prop A->B", "Low Prop A->B"),
#        col=c("red", "black"), lty = 1, cex=1)
# plot(daily_df1$daily_reads_det, type = "l", xlab = "time (days)", ylab = "# reads")
# lines(daily_df2$daily_reads_det, type = "l", col = "red")
# legend(160, 450, legend=c("High Prop A->B", "Low Prop A->B"),
#        col=c("red", "black"), lty = 1, cex=1)
# 
# plot(daily_df1$daily_flightAB_infections, type = "l", xlim = c(50, 100), xlab = "time (days)", ylab = "# infections A->B")
# lines(daily_df2$daily_flightAB_infections, type = "l", col = "red")
# legend(48, 20, legend=c("High Prop A->B", "Low Prop A->B"),
#        col=c("red", "black"), lty = 1, cex=1)
# plot(daily_df1$daily_reads_det, type = "l", xlim = c(50, 100), xlab = "time (days)", ylab = "# reads")
# lines(daily_df2$daily_reads_det, type = "l", col = "red")
# legend(48, 450, legend=c("High Prop A->B ~2%", "Low Prop A->B ~0.2%"),
#        col=c("red", "black"), lty = 1, cex=1)
# 
# plot(output1$seq_reads_virus_aggFlight_det_Out/output1$infected_indiv_shedding_events_Out, ylab = "reads per shedding event", ylim = c(0, 15))
# points(output2$seq_reads_virus_aggFlight_det_Out/output2$infected_indiv_shedding_events_Out, col = "red")
# plot(output1$seq_reads_virus_aggFlight_det_Out/output1$n_inf_flightABOut, ylab = "reads per flight infection", ylim = c(0, 15))
# points(output2$seq_reads_virus_aggFlight_det_Out/output2$n_inf_flightABOut, col = "red")
# 
# plot(100 * output1$n_inf_flightABOut/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB), ylab = "Low Prop Flight Infection Prevalence", pch = 20, col = adjustcolor("black", alpha.f = 0.2))
# lines(100 * (output1$n_inf_flightOut * 1/500)/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB), ylab = "Flight Infection Prevalence", col = "black")
# plot(100 * output2$n_inf_flightABOut/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB), ylab = "High Prop Flight Infection Prevalence", col = adjustcolor("red", alpha.f = 0.2))
# lines(100 * (output2$n_inf_flightOut * 1/50)/(fixed_params$capacity_per_flight * fixed_params$num_flightsAB), ylab = "Flight Infection Prevalence", col = "red")
# 
# plot(100 * output1$n_inf_flightABOut)
# 
# plot(100 * output1$n_inf_flightOut)


# plot(output2$n_inf_flightABOut, output2$infected_indiv_shedding_events_Out, col = "red")
# points(output1$n_inf_flightABOut, output1$infected_indiv_shedding_events_Out)

####
# The problem is broadly as follows:
#   As expected, when we increase the proportion of all flights that are going from Location A to Location B, 
#   we get an increase in the total number of infections travelling from Location A to Location. 
#   However, we're also increasing the total number of individuals flying from Location A to Location B, and so 
#   the proportion of all those travelling from location A to location B stays approximately the same.  
#   I.e. when we increase proportion of all flights A -> B:
#     Number of Infections A -> B: Increases
#     Number of Uninfected A -> B: Increases
#     Prop Travellers A -> B Who Are Infected: Stays ~ The Same
# Within the model, it's the relative number of infected individuals to uninfected individuals that dictates the
# relative abundance of NA of interest to NA not of interest in wastewater, and determines the number of reads
# generated for a given sequencing amount. So that explains (broadly) why we don't see a decrease in time to detection 
# within the current model framework as we increase the proportion/number of flights going from Location A -> Location B.
#
# The next question is why we see (a slightly) SHORTER time-to-detection when we have the lowest number of flights 
# going from Location A -> Location B. When we have very few (e.g. only 1 flight) travelling from Location A -> Location B,
# a single infected individual represents a larger fraction of the A -> B travelling population. As a result, their 
# proportional contribution to the wastewater is larger and for a given fixed sequencing amount, results in a higher number of reads. 
# E.g. compare 1 infected person on one flight of 200 (prev = 0.5%), vs 2 infected people across 10 flights of 200 (prev = 0.1%).
# For a fixed sequencing volume, the former individual's shedding event generates more reads relating to NA of interest than the latter (because the latter
# has a higher proportion of uninfected people shedding NA not of interest into the system). In the first example, a single infection 
# contributes more, and so when you randomly get 1 or 2 infections travelling early on in the epidemic, these contribute 
# enough reads to count as "detection" under our current choice of threshold (5). Relatedly the time-series of reads for the 
# smaller number of flights is way more noisy but because we've set this as a binary yes/no detection based on 5 reads,
# noise is beneficial i.e. our definition of success is happy if randomly you have a spike and then several days of 0s. The
# larger number of flights version wouldn't have as much noise - it wouldn't (and doesn't) spike as high; but it also doesn't dip
# as low afterwards. The overall number of reads generated is the same across the entirety of the model run - which makes
# sense; in the larger flight volume example, we have more infected people travelling - they contribute fewer reads per shedding
# event, but there are more them. Whilst overall number of read generated across entire model run is similar though,
# there are important differences in the early stages of the epidemic; which have material consequences for detectability (as defined
# here) under different flight volume scenarios. 

# This doesn't capture our (or at least my) intuition, and I think it should, but I don't know how to change
# things so that it does. Currently the equation that converts infections -> reads only takes into account the 
# relative number of # infected individuals and uninfected individuals (or more precisely, their respective 
# shedding events and the # amount of nucleic acid released per shedding event). But my sense is that for the relative 
# proportions, # we should want more infections contributing to our system (and concomitantly, more uninfected 
# contributions) - # this is because we only sample a tiny fraction of the overall system (i.e. wastewater) and more 
# individuals # contributing to it means a higher concentration of viral NA, and a reduced chance that we don't randomly get 
# nothing/below our LOD in a sample just by chance (the chance of this is higher when we have fewer individuals
# contributing/shedding into the system. 