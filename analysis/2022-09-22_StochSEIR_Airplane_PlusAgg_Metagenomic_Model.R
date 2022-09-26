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
  capacity_per_flight <- user()    # Capacity of a single flight 
  num_flights <- user()            # Number of flights per day
  num_flightsAB <- user()          # Number of flights per day from Location A to NAO Location (Location B)

  ## Calculating the number of infected people taking a flight based on number of infections and airport/airplane parameters
  ### Note for rhyper calls: 1st param = # white balls; 2nd param # black balls; 3rd param # balls drawn; the output is the number of white balls that get drawn in this context
  n_inf_flight <- rhyper(I, population_size - I, dt * num_flights * capacity_per_flight) 
  
  n_inf_flightAB <- rhyper(n_inf_flight, num_flights * capacity_per_flight * dt - n_inf_flight, num_flightsAB * capacity_per_flight * dt) 
  # n_inf_flight has an implicit dt based on how it's defined above. So not needed in first term. 
  # 2nd arg divided by 3rd arg is ~equivalent to p for binomial, so dt in both (or neither) doesn't make any huge practical diff (especially when 1st arg << 2nd & 3rd args)

  # Distributing (detectable) infections that take a flight (n_inf_take_flight) across all the possible flights they could get on
  capacity_individual_flight[] <- capacity_per_flight  
  n_inf_specific_flightAB[] <- rmhyper(n_inf_flightAB, capacity_individual_flight)  

  ### Stochastic Model Updates for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  update(n_inf_flightOut) <- n_inf_flight
  update(n_inf_flightABOut) <- n_inf_flightAB
  update(n_inf_specific_flightABOut[]) <- n_inf_specific_flightAB[i]

  ### Initial Values for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  initial(n_inf_flightOut) <- 0
  initial(n_inf_flightABOut) <- 0
  initial(n_inf_specific_flightABOut[]) <- 0
  
  ### Dimension Assignments for Outputs Relevant to Surveillance
  dim(capacity_individual_flight) <- num_flightsAB
  dim(n_inf_specific_flightAB) <- num_flightsAB
  dim(n_inf_specific_flightABOut) <- num_flightsAB
  
  ################## Metagenomic Sequencing Model ##################
  
  ## Metagenomic and Sequencing Parameters
  shedding_freq <- user()              # Average number of defecation events per person per flight
  virus_shed <- user()                 # Average amount of material shed per event for our virus of interest (defecation)
  non_virus_shed <- user()             # Average amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  met_bias <- user()                   # Bias term for the metagenomic model
  seq_tot <- user()                    # Total amount of sequencing that is done
  samp_frac_indivFlight <- user()     # Fraction of each flight's wastewater that would be sampled if individual flights surveilled 
  samp_frac_aggFlight <- user()      # Fraction of all flight's wastewater that would be sampled if individual flight wastewater pooled and then sampled (via triturator) 
                                
  # Calculating the number of shedding events from infected and uninfected individuals on each airplane
  infected_indiv_shedding_events[] <- rpois(n_inf_specific_flightAB[i] * shedding_freq)
  uninfected_indiv_shedding_events[] <- rpois((capacity_individual_flight[i] - n_inf_specific_flightAB[i]) * shedding_freq)
  
  ## PER FLIGHT CALCULATIONS
  ### Calculating amount and concentration of nucleic acid shed into wastewater on each flight
  amount_virus_indivFlight[] <- infected_indiv_shedding_events[i] * virus_shed 
  amount_non_virus_indivFlight[] <- (uninfected_indiv_shedding_events[i] + infected_indiv_shedding_events[i]) * non_virus_shed
  sample_amount_virus_indivFlight[] <- rbinom(amount_virus_indivFlight[i], samp_frac_indivFlight)
  sample_amount_non_virus_indivFlight[] <- rbinom(amount_non_virus_indivFlight[i], samp_frac_aggFlight)

  ### NOTE: CONSIDER REVERTING BACK TO PREVIOUS FORMULATION
  ##        I.E. DON'T HAVE THE FRAC THING GOING ON FOR SAMPLING;
  ##        AND JUST HAVE THE DETERMINISTIC VERSION OF THE PROCESS,
  ##        AND THEN THE POISSON ON TOP. 
  
  ### Converting this nucleic acid abundance into sequencing reads for each flight
  seq_reads_virus_indivFlight[] <- if(sample_amount_virus_indivFlight[i] == 0 && sample_amount_non_virus_indivFlight[i] == 0) 0 else seq_tot * (sample_amount_virus_indivFlight[i] * met_bias)/((sample_amount_virus_indivFlight[i] * met_bias) + sample_amount_non_virus_indivFlight[i])
  seq_reads_non_virus_indivFlight[] <- seq_tot - seq_reads_virus_indivFlight[i]
  # stoch_seq_reads_virus_indivFlight[] <- rpois(seq_reads_virus_indivFlight[i]) 
  # stoch_seq_reads_non_virus_indivFlight[] <- seq_tot - stoch_seq_reads_virus_indivFlight[i]

  ## AGGREGATED FLIGHT CALCULATIONS (E.G. SAMPLING FROM TRITURATOR)
  ### Note that the binomial call accounts for fact that we're sampling a very small volume relative to total volume of the wastewater system, and therefore
  ### account for the fact that when there's very few virla genome copies, we might just not get any in a particular sample just by chance.
  
  ### Calculating amount and concentration of nucleic acid shed into aggregated flight wastewater
  amount_virus_aggFlight <- sum(amount_virus_indivFlight)
  amount_non_virus_aggFlight <- sum(amount_non_virus_indivFlight)
  sample_amount_virus_aggFlight <- rbinom(amount_virus_aggFlight, samp_frac_aggFlight)
  sample_amount_non_virus_aggFlight <- rbinom(amount_non_virus_aggFlight, samp_frac_aggFlight)
  
  ### Converting this nucleic acid abundance into sequencing reads for aggregated flight wastewater
  seq_reads_virus_aggFlight <- if(sample_amount_virus_aggFlight == 0 && sample_amount_non_virus_aggFlight == 0) 0 else seq_tot * (sample_amount_virus_aggFlight * met_bias)/((sample_amount_virus_aggFlight * met_bias) + sample_amount_non_virus_aggFlight)
  seq_reads_non_virus_aggFlight <- seq_tot - seq_reads_virus_aggFlight
  # stoch_seq_reads_virus_aggFlight <- rpois(seq_reads_virus_aggFlight) 
  # stoch_seq_reads_non_virus_aggFlight <- seq_tot - stoch_seq_reads_virus_aggFlight
  
  ### Stochastic Model Updates for Outputs Relevant to Metagenomic Sequencing
  update(amount_virus_indivFlight_Out[]) <- amount_virus_indivFlight[i]
  update(amount_non_virus_indivFlight_Out[]) <- amount_non_virus_indivFlight[i]
  update(sample_amount_virus_indivFlight_Out[]) <- sample_amount_virus_indivFlight[i]
  update(sample_amount_non_virus_indivFlight_Out[]) <- sample_amount_non_virus_indivFlight[i]
  update(seq_reads_virus_indivFlight_Out[]) <- seq_reads_virus_indivFlight[i]
  update(seq_reads_non_virus_indivFlight_Out[]) <- seq_reads_non_virus_indivFlight[i] 
  update(amount_virus_aggFlight_Out) <- amount_virus_aggFlight
  update(amount_non_virus_aggFlight_Out) <- amount_non_virus_aggFlight
  update(sample_amount_virus_aggFlight_Out) <- sample_amount_virus_aggFlight
  update(sample_amount_non_virus_aggFlight_Out) <- sample_amount_non_virus_aggFlight
  update(seq_reads_virus_aggFlight_Out) <- seq_reads_virus_aggFlight
  update(seq_reads_non_virus_aggFlight_Out) <- seq_reads_non_virus_aggFlight
  
  ### Initial Values for Outputs Relevant to Metagenomic Sequencing
  initial(amount_virus_indivFlight_Out[]) <- 0
  initial(amount_non_virus_indivFlight_Out[]) <- 0
  initial(sample_amount_virus_indivFlight_Out[]) <- 0
  initial(sample_amount_non_virus_indivFlight_Out[]) <- 0
  initial(seq_reads_virus_indivFlight_Out[]) <- 0
  initial(seq_reads_non_virus_indivFlight_Out[]) <- 0
  initial(amount_virus_aggFlight_Out) <- 0
  initial(amount_non_virus_aggFlight_Out) <- 0
  initial(sample_amount_virus_aggFlight_Out) <- 0
  initial(sample_amount_non_virus_aggFlight_Out) <- 0
  initial(seq_reads_virus_aggFlight_Out) <- 0
  initial(seq_reads_non_virus_aggFlight_Out) <- 0
  
  ### Dimension Assignments for Outputs Relevant to Metagenomic Sequencing
  dim(infected_indiv_shedding_events) <- num_flightsAB
  dim(uninfected_indiv_shedding_events) <- num_flightsAB
  dim(amount_virus_indivFlight) <- num_flightsAB
  dim(amount_virus_indivFlight_Out) <- num_flightsAB
  dim(amount_non_virus_indivFlight) <- num_flightsAB
  dim(amount_non_virus_indivFlight_Out) <- num_flightsAB
  dim(sample_amount_virus_indivFlight) <- num_flightsAB
  dim(sample_amount_virus_indivFlight_Out) <- num_flightsAB
  dim(sample_amount_non_virus_indivFlight) <- num_flightsAB
  dim(sample_amount_non_virus_indivFlight_Out) <- num_flightsAB
  dim(seq_reads_virus_indivFlight) <- num_flightsAB
  dim(seq_reads_virus_indivFlight_Out) <- num_flightsAB
  dim(seq_reads_non_virus_indivFlight) <- num_flightsAB
  dim(seq_reads_non_virus_indivFlight_Out) <- num_flightsAB

  ################# Miscellaneous Model Stuff ##################
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt 
  
})

# Specifying and Running the Model - Checking how long we need to run the model for R0 = 1.5
tic()
mod <- stoch_seir_dust$new(# Epidemiological Parameters
  beta = 0.6, gamma = 0.25, sigma = 0.2, population_size = 10^6, start_infections = 10,

  # Flight Parameters
  capacity_per_flight = 250, num_flights = 40, num_flightsAB = 20,

  # Sequencing Parameters
  shedding_freq = 1, virus_shed = 10000, non_virus_shed = 10000 * 10^6, met_bias = 1,
  seq_tot = 10^9,
  samp_frac_indivFlight = 1/800,
  samp_frac_aggFlight = 1/(20*800),

  # Miscellaenous Parameters
  dt = 0.2)

output <- mod$run(1:(1500/0.2))
output2 <- mod$transform_variables(output)
toc()

# # Plotting Overall Epidemic Dynamics
state_vars_df <- data.frame(time = output2$time, S = output2$S, E = output2$E, I = output2$I, R = output2$R) %>%
  pivot_longer(-time, names_to = "state_var", values_to = "value")
ggplot(state_vars_df, aes(x = time, y = value, colour = state_var)) +
  geom_line() +
  labs(x = "Time (Days)", y = "Number of Individuals in Each State") +
  lims(x = c(0, 200))

# Metric for Successful Detection
num_reads <- 5             # Number of reads required for successful detection

## Model Run Related Parameters
dt <- 0.2
seed <- rpois(100, 200) * rpois(100, 200) * rpois(100, 20) # Seed for each stochastic realisation
stochastic_sim <- 25                                       # # of stochastic realisations for each parameter value

## Fixed Parameters 
vol_flight_ww <- 800       # Assumption that requires updating
sample_flight_ww <- 1      # Assumption that requires updating
population_size <- 10^6    # Assumption that requires updating
capacity_per_flight <- 250 # Assumption that requires updating
non_virus_shed <- 2000000  # Assumption that requires updating
start_infections <- 10
met_bias <- 1

## Parameters To Vary

### Epidemiological
R0 <- seq(1.5, 3, 0.25)         # Range 1.50 to 3, increments of 0.25, 7 values
gamma <- 1/4
sigma <- 1/5
beta <- R0 * sigma 
end <- c(250, 200, 150, 125, 125, 100, 100) # Time to run the model for each of the beta values (computational cost saving trick, don't need to run model
                                            # as long when R0 is high and epidemic is over quickly

### Shedding 
shedding_freq <- 1                   # Range 0.2 to 3, increments of 0.4 - 8 values
ratio_virus_to_non_virus <- 1/10^7   # Range 1/10^6 to 1/10^8 - 8 values

### Flight
proportion_flying <- seq(0.002, 0.010, 0.002)       # Range 0.002 to 0.010 (i.e. 0.2% to 1% population flying daily), increments of 0.002 - 5 values
num_flights <- round(proportion_flying * population_size / capacity_per_flight, 0)          
proportion_AB <- seq(0.005, 0.025, 0.005)           # Range 0.005 to 0.025 (i.e. 1 in every 200 flights to 1 in every 40 flights), increment of 0.005 - 5 values
num_flightsAB <- num_flights * proportion_AB
function_pop <- proportion_flying * proportion_AB

### Sequencing
seq_tot <- round(lseq(10^7, 10^9, 8)) # Range TBD but 8 values

7 * 8 * 8 * 5 * 5 * 8

7 * # R0 values
  8 * # Viral Ratio Values
    5 * # Proportion Flying Values 
      8  # Sequencing Total Values


##### Sensitivity Analysis - Epidemiological Parameters (Beta) 
tic()
beta_sens <- beta
R0_df <- data.frame(beta = beta_sens, R0 = beta_sens/sigma)
beta_ttd_output <- data.frame(beta = rep(beta_sens, each = stochastic_sim), 
                              stochastic_realisation = 1:stochastic_sim, 
                              num_reads = num_reads, 
                              metric = "time_to_detection",
                              first_flight_reads = NA, 
                              avg_flight_reads = NA, 
                              agg_flight_reads = NA,  
                              first_flight_inf = NA, 
                              avg_flight_inf = NA)
beta_cuminf_output <- data.frame(beta = rep(beta_sens, each = stochastic_sim), 
                                 stochastic_realisation = 1:stochastic_sim, 
                                 num_reads = num_reads, 
                                 metric = "cumulative_incidence",
                                 first_flight_reads = NA,
                                 avg_flight_reads = NA, 
                                 agg_flight_reads = NA, 
                                 first_flight_inf = NA, 
                                 avg_flight_inf = NA)
counter <- 1
for (i in 1:length(beta_sens)) {
  
  for (j in 1:stochastic_sim) {
    
    # Running the Model
    set.seed(seed[j])
    mod <- stoch_seir_dust$new(# Epidemiological Parameters
      beta = beta_sens[i], gamma = gamma, sigma = sigma, population_size = population_size, start_infections = start_infections,
      
      # Flight Parameters
      capacity_per_flight = capacity_per_flight, num_flights = num_flights, num_flightsAB = num_flightsAB, 
      
      # Sequencing Parameters
      shedding_freq = shedding_freq, virus_shed = non_virus_shed * ratio_virus_to_non_virus, non_virus_shed = non_virus_shed, met_bias = met_bias, 
      seq_tot = seq_tot[8], 
      samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
      samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
      
      # Miscellaenous Parameters
      dt = dt)
    
    # Extracting Output
    output <- mod$run(1:(end[i]/dt))
    output2 <- mod$transform_variables(output)
    
    # Calculating Time to Detection
    beta_ttd_output$first_flight_reads[counter] <- min(time_to_detection_fun(output2, num_reads, "reads", "individual"))
    beta_ttd_output$avg_flight_reads[counter] <- mean(time_to_detection_fun(output2, num_reads, "reads", "individual")) # check whether this is correct or need to find when average line passes 5 for first time 
    beta_ttd_output$agg_flight_reads[counter] <- time_to_detection_fun(output2, num_reads, "reads", "aggregated")
    beta_ttd_output$first_flight_inf[counter] <- min(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    beta_ttd_output$avg_flight_inf[counter] <- mean(time_to_detection_fun(output2, num_reads, "infections", "individual")) # check whether this is correct or need to find when average line passes 5 for first time 
    
    # Calculating Cumulative Incidence at Time to Detection
    beta_cuminf_output$first_flight_reads[counter] <- if(!is.na(beta_ttd_output$first_flight_reads[counter])) sum(output2$n_SE_Output[1:(beta_ttd_output$first_flight_reads[counter]/dt)])/population_size else NA
    beta_cuminf_output$avg_flight_reads[counter] <- if(!is.na(beta_ttd_output$avg_flight_reads[counter])) sum(output2$n_SE_Output[1:(beta_ttd_output$avg_flight_reads[counter]/dt)])/population_size else NA
    beta_cuminf_output$agg_flight_reads[counter] <- if(!is.na(beta_ttd_output$agg_flight_reads[counter])) sum(output2$n_SE_Output[1:(beta_ttd_output$agg_flight_reads[counter]/dt)])/population_size else NA
    beta_cuminf_output$first_flight_inf[counter] <- if(!is.na(beta_ttd_output$first_flight_inf[counter])) sum(output2$n_SE_Output[1:(beta_ttd_output$first_flight_inf[counter]/dt)])/population_size else NA
    beta_cuminf_output$avg_flight_inf[counter] <- if(!is.na(beta_ttd_output$avg_flight_inf[counter])) sum(output2$n_SE_Output[1:(beta_ttd_output$avg_flight_inf[counter]/dt)])/population_size else NA
    
    counter <- counter + 1   
  }
  print(i)
}
toc()
saveRDS(list(ttd = beta_ttd_output, cuminf = beta_cuminf_output), "outputs/beta_sensitivity_analysis.rds")

## Summarising and Plotting the Output
beta_output <- readRDS("outputs/beta_sensitivity_analysis.rds")
beta_df_ttd <- beta_output$ttd %>%
  pivot_longer(-c(beta, num_reads, stochastic_realisation, metric), values_to = "time_to_detection", names_to = "method_calc") %>%
  left_join(R0_df, by = "beta") %>%
  select(-metric)
beta_df_cuminf <- beta_output$cuminf %>%
  pivot_longer(-c(beta, num_reads, stochastic_realisation, metric), values_to = "cumulative_incidence", names_to = "method_calc") %>%
  left_join(R0_df, by = "beta") %>%
  select(-metric)
beta_df_overall <- beta_df_ttd %>%
  left_join(beta_df_cuminf, by = c("beta", "R0", "stochastic_realisation", "num_reads", "method_calc")) %>%
  mutate(cumulative_incidence = 100 * cumulative_incidence) %>%
  filter(!grepl("*inf", method_calc))

beta_df_summary <- beta_df_overall %>%
  group_by(beta, R0, method_calc) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = min(time_to_detection, na.rm = TRUE),
            upper_ttd = max(time_to_detection, na.rm = TRUE),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = min(cumulative_incidence, na.rm = TRUE),
            upper_cuminf = max(cumulative_incidence, na.rm = TRUE),
            num_reached = paste0(stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (stochastic_sim - sum(is.na(time_to_detection)))/stochastic_sim) %>%
  group_by(method_calc) %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         scaling_factor_ttd = max(upper_ttd, na.rm = TRUE)/max(perc_reached, na.rm = TRUE), 
         scaled_ttd = perc_reached * max(upper_ttd, na.rm = TRUE)/max(perc_reached, na.rm = TRUE),
         scaling_factor_cuminf = max(upper_cuminf, na.rm = TRUE)/max(perc_reached, na.rm = TRUE), 
         scaled_cuminf = perc_reached * max(upper_cuminf, na.rm = TRUE)/max(perc_reached, na.rm = TRUE))

df_split <- beta_df_summary %>% 
  mutate(var = R0) %>%
  split(.$method_calc)

make_plot <- function(df, scaling_factor1, scaling_factor2, title, variable, transform) {
  
  if (transform == "log") {
    temp1 <- ggplot(data = df) +
      geom_line(aes(x = var, y = avg_ttd)) +
      geom_line(aes(x = var, y = scaled_ttd)) +
      geom_ribbon(aes(x = var, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), alpha = 0.2) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(limits = c(0, NA), 
                         sec.axis = sec_axis(~./scaling_factor1,
                                             breaks = c(0, 20, 40, 60, 80, 100),
                                             name = "% Simulations With\nSuccessful Detection")) +
      labs(x = paste0("Parameter Value - ", variable), y = "Time to Detection") +
      theme_bw() +
      labs(title = title)
    
    temp2 <- ggplot(data = df) +
      geom_line(aes(x = var, y = avg_cuminf)) +
      geom_line(aes(x = var, y = scaled_cuminf)) +
      geom_ribbon(aes(x = var, y = avg_cuminf, ymin = lower_cuminf, ymax = upper_cuminf), alpha = 0.2) +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(limits = c(0, NA), 
                         sec.axis = sec_axis(~./scaling_factor2,
                                             breaks = c(0, 20, 40, 60, 80, 100),
                                             name = "% Simulations With\nSuccessful Detection")) +
      labs(x = paste0("Parameter Value - ", variable), y = "Cumulative Incidence @\nTime of Detection") +
      theme_bw()
  } else {
    temp1 <- ggplot(data = df) +
      geom_line(aes(x = var, y = avg_ttd)) +
      geom_line(aes(x = var, y = scaled_ttd)) +
      geom_ribbon(aes(x = var, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), alpha = 0.2) +
      scale_y_continuous(limits = c(0, NA), 
                         sec.axis = sec_axis(~./scaling_factor1,
                                             breaks = c(0, 20, 40, 60, 80, 100),
                                             name = "% Simulations With\nSuccessful Detection")) +
      labs(x = paste0("Parameter Value - ", variable), y = "Time to Detection") +
      theme_bw() +
      labs(title = title)
    
    temp2 <- ggplot(data = df) +
      geom_line(aes(x = var, y = avg_cuminf)) +
      geom_line(aes(x = var, y = scaled_cuminf)) +
      geom_ribbon(aes(x = var, y = avg_cuminf, ymin = lower_cuminf, ymax = upper_cuminf), alpha = 0.2) +
      scale_y_continuous(limits = c(0, NA), 
                         sec.axis = sec_axis(~./scaling_factor2,
                                             breaks = c(0, 20, 40, 60, 80, 100),
                                             name = "% Simulations With\nSuccessful Detection")) +
      labs(x = paste0("Parameter Value - ", variable), y = "Cumulative Incidence @\nTime of Detection") +
      theme_bw()
  }
  
  temp3 <- temp1 + temp2 
  
  return(temp3)
}

scaling_factor_list_ttd <- list(length = length(df_split))
scaling_factor_list_cuminf <- list(length = length(df_split))
for (i in 1:length(df_split)) {
  scaling_factor_list_ttd[[i]] <- unique(beta_df_summary$scaling_factor_ttd)[i]  
  scaling_factor_list_cuminf[[i]] <- unique(beta_df_summary$scaling_factor_cuminf)[i]  
}

### CHECK THIS IS THE RIGHT ORDERING
list_titles <- list("Aggregated Flight Waste\n5 Reads for Detection",
                    "Individual Flight Waste\n1st Flight With 5 Reads",
                    "Individual Flight Waste\nAverage Flight With 5 Reads")

x <- purrr::pmap(list(df = df_split, 
                      scaling_factor1 = scaling_factor_list_ttd, 
                      scaling_factor2 = scaling_factor_list_cuminf,
                      title = list_titles,
                      variable = "R0",
                      transform = "none"), make_plot) %>% 
  patchwork::wrap_plots(ncol = 1, tag_level = "keep") 
x
ggsave("test.pdf", x, height = 8, width = 10)

test <- beta_df_overall %>% 
  pivot_longer(-c(beta, stochastic_realisation, num_reads, method_calc, R0), 
               values_to = "value", names_to = "metric") %>%
  group_by(R0, beta, method_calc, metric) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            lower = min(value, na.rm = TRUE),
            upper = max(value, na.rm = TRUE),
            num_reached = paste0(stochastic_sim - sum(is.na(value))),
            perc_reached = 100 * (stochastic_sim - sum(is.na(value)))/stochastic_sim) %>%
  group_by(method_calc) %>%
  mutate(lower = ifelse(is.infinite(lower), NaN, lower),
         upper = ifelse(is.infinite(upper), NaN, upper),
         scaling_factor = max(upper, na.rm = TRUE)/100) %>%
  ungroup() %>%
  mutate(scaled = perc_reached * scaling_factor)

a <- ggplot(data = test[test$metric == "cumulative_incidence", ]) +
  geom_line(aes(x = R0, y = avg, col = method_calc)) +
  geom_ribbon(aes(x = R0, y = avg, ymin = lower, ymax = upper, fill = method_calc), alpha = 0.2) +
  labs(x = "Parameter Value - R0", y = "Cumulative Incidence (%)") +
  theme_bw() +
  theme(plot.margin = margin(0, 1, 1, 1),
        legend.position = "none")

a2 <- ggplot(data = test[test$metric == "cumulative_incidence", ]) +
  geom_line(aes(x = R0, y = perc_reached, col = method_calc), position=position_jitter(w=0.025, h=0)) +
  labs(x = "", y = "% Detect\nSuccess") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(1, 1, 0, 1))

a3 <- a2 + a +
  plot_layout(nrow = 2, heights = c(1, 5))


b <- ggplot(data = test[test$metric == "time_to_detection", ]) +
  geom_line(aes(x = R0, y = avg, col = method_calc)) +
  geom_ribbon(aes(x = R0, y = avg, ymin = lower, ymax = upper, fill = method_calc), alpha = 0.2) +
  labs(x = "Parameter Value - R0", y = "Time to Detection (Days)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "none")
b3 <- a2 + b +
  plot_layout(nrow = 2, heights = c(1, 5))

legend <- cowplot::get_legend(b + theme(legend.position = "right"))

t <- cowplot::plot_grid(b3, a3, nrow = 1)
u <- cowplot::plot_grid(t, legend, ncol = 2, rel_widths = c(10, 1.5))


d3 <- b3 + a3 +
  plot_layout(nrow = 2, ncol = 2)



d3

layout <- "AABB
           AABB"

guides = 'collect',


##### Sensitivity Analysis - Sequencing Parameters (Seq_Total)
seq_tot_sens <- round(lseq(10^8, 10^10, 20))
seq_ttd_output <- data.frame(seq_total = rep(seq_tot_sens, each = stochastic_sim), 
                             stochastic_realisation = 1:stochastic_sim, 
                             num_reads = num_reads, 
                             metric = "time_to_detection",
                             first_flight_reads = NA, 
                             avg_flight_reads = NA, 
                             agg_flight_reads = NA,  
                             first_flight_inf = NA, 
                             avg_flight_inf = NA)
seq_cuminf_output <- data.frame(seq_total = rep(seq_tot_sens, each = stochastic_sim), 
                                stochastic_realisation = 1:stochastic_sim, 
                                num_reads = num_reads, 
                                metric = "cumulative_incidence",
                                first_flight_reads = NA,
                                avg_flight_reads = NA, 
                                agg_flight_reads = NA, 
                                first_flight_inf = NA, 
                                avg_flight_inf = NA)
counter <- 1
beta <- 0.4
for (i in 1:length(seq_tot_sens)) {
  
  for (j in 1:stochastic_sim) {
    
    # Running the Model
    set.seed(seed[j])
    mod <- stoch_seir_dust$new(# Epidemiological Parameters
      beta = beta, gamma = gamma, sigma = sigma, population_size = population_size, start_infections = start_infections,
      
      # Flight Parameters
      capacity_per_flight = capacity_per_flight, num_flights = num_flights, num_flightsAB = num_flightsAB, 
      
      # Sequencing Parameters
      shedding_freq = shedding_freq, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = met_bias, 
      seq_tot = seq_tot_sens[i], 
      samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
      samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
      
      # Miscellaenous Parameters
      dt = dt)
    
    # Extracting Output
    output <- mod$run(1:(end/dt))
    output2 <- mod$transform_variables(output)
    
    # Calculating Time to Detection
    seq_ttd_output$first_flight_reads[counter] <- min(time_to_detection_fun(output2, num_reads, "reads", "individual"))
    seq_ttd_output$avg_flight_reads[counter] <- mean(time_to_detection_fun(output2, num_reads, "reads", "individual")) # check whether this is correct or need to find when average line passes 5 for first time 
    seq_ttd_output$agg_flight_reads[counter] <- time_to_detection_fun(output2, num_reads, "reads", "aggregated")
    seq_ttd_output$first_flight_inf[counter] <- min(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    seq_ttd_output$avg_flight_inf[counter] <- mean(time_to_detection_fun(output2, num_reads, "infections", "individual")) # check whether this is correct or need to find when average line passes 5 for first time 
    
    # Calculating Cumulative Incidence at Time to Detection
    seq_cuminf_output$first_flight_reads[counter] <- if(!is.na(seq_ttd_output$first_flight_reads[counter])) sum(output2$n_SE_Output[1:(seq_ttd_output$first_flight_reads[counter]/dt)])/population_size else NA
    seq_cuminf_output$avg_flight_reads[counter] <- if(!is.na(seq_ttd_output$avg_flight_reads[counter])) sum(output2$n_SE_Output[1:(seq_ttd_output$avg_flight_reads[counter]/dt)])/population_size else NA
    seq_cuminf_output$agg_flight_reads[counter] <- if(!is.na(seq_ttd_output$agg_flight_reads[counter])) sum(output2$n_SE_Output[1:(seq_ttd_output$agg_flight_reads[counter]/dt)])/population_size else NA
    seq_cuminf_output$first_flight_inf[counter] <- if(!is.na(seq_ttd_output$first_flight_inf[counter])) sum(output2$n_SE_Output[1:(seq_ttd_output$first_flight_inf[counter]/dt)])/population_size else NA
    seq_cuminf_output$avg_flight_inf[counter] <- if(!is.na(seq_ttd_output$avg_flight_inf[counter])) sum(output2$n_SE_Output[1:(seq_ttd_output$avg_flight_inf[counter]/dt)])/population_size else NA
    
    counter <- counter + 1   
  }
  print(i)
}
saveRDS(list(ttd = seq_ttd_output, cuminf = seq_cuminf_output), "outputs/seq_sensitivity_analysis.rds")

## Summarising and Plotting the Output
seq_output <- readRDS("outputs/seq_sensitivity_analysis.rds")
seq_df_ttd <- seq_output$ttd %>%
  pivot_longer(-c(seq_total, num_reads, stochastic_realisation, metric), values_to = "time_to_detection", names_to = "method_calc") %>%
  select(-metric)
seq_df_cuminf <- seq_output$cuminf %>%
  pivot_longer(-c(seq_total, num_reads, stochastic_realisation, metric), values_to = "cumulative_incidence", names_to = "method_calc") %>%
  select(-metric)
seq_df_overall <- seq_df_ttd %>%
  left_join(seq_df_cuminf, by = c("seq_total", "stochastic_realisation", "num_reads", "method_calc")) %>%
  mutate(cumulative_incidence = 100 * cumulative_incidence) %>%
  filter(!grepl("*inf", method_calc))

seq_df_summary <- seq_df_overall %>%
  group_by(seq_total, method_calc) %>%
  summarise(avg_ttd = mean(time_to_detection, na.rm = TRUE),
            lower_ttd = min(time_to_detection, na.rm = TRUE),
            upper_ttd = max(time_to_detection, na.rm = TRUE),
            avg_cuminf = mean(cumulative_incidence, na.rm = TRUE),
            lower_cuminf = min(cumulative_incidence, na.rm = TRUE),
            upper_cuminf = max(cumulative_incidence, na.rm = TRUE),
            num_reached = paste0(stochastic_sim - sum(is.na(time_to_detection))),
            perc_reached = 100 * (stochastic_sim - sum(is.na(time_to_detection)))/stochastic_sim) %>%
  group_by(method_calc) %>%
  mutate(lower_ttd = ifelse(is.infinite(lower_ttd), NaN, lower_ttd),
         upper_ttd = ifelse(is.infinite(upper_ttd), NaN, upper_ttd),
         lower_cuminf = ifelse(is.infinite(lower_cuminf), NaN, lower_cuminf),
         upper_cuminf = ifelse(is.infinite(upper_cuminf), NaN, upper_cuminf),
         scaling_factor_ttd = max(upper_ttd, na.rm = TRUE)/max(perc_reached, na.rm = TRUE), 
         scaled_ttd = perc_reached * max(upper_ttd, na.rm = TRUE)/max(perc_reached, na.rm = TRUE),
         scaling_factor_cuminf = max(upper_cuminf, na.rm = TRUE)/max(perc_reached, na.rm = TRUE), 
         scaled_cuminf = perc_reached * max(upper_cuminf, na.rm = TRUE)/max(perc_reached, na.rm = TRUE))

df_split <- seq_df_summary %>% 
  mutate(var = seq_total) %>%
  split(.$method_calc)

scaling_factor_list_ttd <- list(length = length(df_split))
scaling_factor_list_cuminf <- list(length = length(df_split))
for (i in 1:length(df_split)) {
  scaling_factor_list_ttd[[i]] <- unique(seq_df_summary$scaling_factor_ttd)[i]  
  scaling_factor_list_cuminf[[i]] <- unique(seq_df_summary$scaling_factor_cuminf)[i]  
}

list_titles <- list("Aggregated Flight Waste\n5 Reads for Detection",
                    "Individual Flight Waste\nAverage Flight With 5 Reads",
                    "Individual Flight Waste\n1st Flight With 5 Reads")

x <- purrr::pmap(list(df = df_split, 
                      scaling_factor1 = scaling_factor_list_ttd, 
                      scaling_factor2 = scaling_factor_list_cuminf,
                      title = list_titles,
                      variable = "Total Seq",
                      transform = "log"), make_plot) %>% 
  patchwork::wrap_plots(ncol = 1, tag_level = "keep") 
x
ggsave("test.pdf", x, height = 8, width = 10)


test <- seq_df_overall %>% 
  pivot_longer(-c(seq_total, stochastic_realisation, num_reads, method_calc), 
               values_to = "value", names_to = "metric") %>%
  group_by(seq_total, method_calc, metric) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            lower = min(value, na.rm = TRUE),
            upper = max(value, na.rm = TRUE),
            num_reached = paste0(stochastic_sim - sum(is.na(value))),
            perc_reached = 100 * (stochastic_sim - sum(is.na(value)))/stochastic_sim) %>%
  group_by(method_calc) %>%
  mutate(lower = ifelse(is.infinite(lower), NaN, lower),
         upper = ifelse(is.infinite(upper), NaN, upper),
         scaling_factor = max(upper, na.rm = TRUE)/100) %>%
  ungroup() %>%
  mutate(scaled = perc_reached * scaling_factor)

a <- ggplot(data = test[test$metric == "cumulative_incidence", ]) +
  geom_line(aes(x = seq_total, y = avg, col = method_calc)) +
  geom_ribbon(aes(x = seq_total, y = avg, ymin = lower, ymax = upper, fill = method_calc), alpha = 0.2) +
  labs(x = "Parameter Value - Seq_Total", y = "Cumulative Incidence (%)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") 
b <- ggplot(data = test[test$metric == "time_to_detection", ]) +
  geom_line(aes(x = seq_total, y = avg, col = method_calc)) +
  geom_ribbon(aes(x = seq_total, y = avg, ymin = lower, ymax = upper, fill = method_calc), alpha = 0.2) +
  labs(x = "Parameter Value - Seq_Total", y = "Time to Detection (Days)") +
  theme_bw() +
  scale_x_continuous(trans = "log10") 

b + a +
  plot_layout(guides = 'collect')

#### OLD AND CHECKING INDIVIDUAL RUNS #### 
# 
# set.seed(2)
# mod <- stoch_seir_dust$new(# Epidemiological Parameters
#                            beta = 0.6, gamma = 1/4, sigma = 1/5, population_size = 10^6, start_infections = 10,
#                            
#                            # Flight Parameters
#                            capacity_per_flight = 250, num_flights = 50, num_flightsAB = num_flightsAB, 
#                            
#                            # Sequencing Parameters
#                            shedding_freq = 1, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = 1, 
#                            seq_tot = 10^9, 
#                            samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
#                            samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
#                            
#                            # Miscellaenous Parameters
#                            dt = dt)
# output <- mod$run(1:(end/dt))
# output2 <- mod$transform_variables(output)
# 
# # Plotting Overall Epidemic Dynamics
# state_vars_df <- data.frame(time = output2$time, S = output2$S, E = output2$E, I = output2$I, R = output2$R) %>%
#   pivot_longer(-time, names_to = "state_var", values_to = "value") 
# ggplot(state_vars_df, aes(x = time, y = value, colour = state_var)) +
#   geom_line() +
#   labs(x = "Time (Days)", y = "Number of Individuals in Each State")
#   
# # Plotting the Number of Infections on Each Flight Per Day
# 
# ## Wrangling data into right format
# num_flights <- mod$contents()[["num_flightsAB"]]
# indiv_flight_infections <- data.frame(time = output2$time, output2$n_inf_specific_flightABOut)
# colnames(indiv_flight_infections) <- c("time", paste0("flight", 1:num_flights))
# airplane_infections_df <- indiv_flight_infections %>%
#   pivot_longer(-time, names_to = "flight_num", values_to = "infections") %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(flight_num, time2) %>%
#   summarise(daily_flight_infections = sum(infections)) %>%
#   mutate(time3 = midpoints(time2))
# 
# agg_flight_infs <- data.frame(time = output2$time, output2$n_inf_flightABOut)
# colnames(agg_flight_infs) <- c("time", "infections")
# agg_inf_df <- agg_flight_infs %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(time2) %>%
#   summarise(daily_infections = sum(infections)/dim(output2$n_inf_specific_flightABOut)[2]) %>%
#   mutate(time3 = midpoints(time2))  
# 
# ### Boxplot 
# ggplot(airplane_infections_df, aes(x = time3, y = daily_flight_infections, group = time3)) +
#   geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
#   geom_jitter(size = 0.5, width = 0.15, alpha = 0.5, height = 0) +
#   theme(legend.position = "none") +
#   lims(x = c(0, end), y = c(0, max(airplane_infections_df$daily_flight_infections))) +
#   labs(x = "Time (Days)", y = "Number of Infections On Each Daily Flight")
# 
# ### Line
# ggplot() +
#   geom_line(data = airplane_infections_df, aes(x = time3, y = daily_flight_infections, colour = flight_num), 
#             fill = NA, outlier.colour = NA) +
#   geom_line(data = agg_inf_df, aes(x = time3, y = daily_infections), 
#             fill = NA) +
#   theme(legend.position = "none") +
#   lims(x = c(0, end), y = c(0, max(airplane_infections_df$daily_flight_infections))) +
#   labs(x = "Time (Days)", y = "Number of Infections On Each Daily Flight")
# 
# # Plotting the Number of Mapped Reads 
# 
# ## Wrangling data into right format
# indiv_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_indivFlight_Out)
# colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))
# airplane_reads_df <- indiv_flight_reads %>%
#   pivot_longer(-time, names_to = "flight_num", values_to = "reads") %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(flight_num, time2) %>%
#   summarise(daily_flight_reads = sum(reads)) %>%
#   mutate(time3 = midpoints(time2))  
# 
# agg_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_aggFlight_Out)
# colnames(agg_flight_reads) <- c("time", "reads")
# agg_reads_df <- agg_flight_reads %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(time2) %>%
#   summarise(daily_reads = sum(reads)) %>%
#   mutate(time3 = midpoints(time2))  
# 
# #### Boxplot
# ggplot() +
#   geom_boxplot(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3), 
#                fill = NA, outlier.colour = NA, alpha = 0.1, col = "grey") +
#   geom_jitter(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3), 
#               size = 0.5, width = 0.15, alpha = 0.5, height = 0, col = "grey") +
#   geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
#   theme(legend.position = "none") +
#   lims(x = c(0, end/2), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
#   labs(x = "Time (Days)", y = "Number of Reads In Each Flight's Wastewater")
# 
# ### Line
# ggplot() +
#   geom_line(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, colour = flight_num), 
#             fill = NA, outlier.colour = NA) +
#   geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
#   theme(legend.position = "none") +
#   lims(x = c(0, 45), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
#   labs(x = "Time (Days)", y = "Number of Reads In Flight Wastewater")


# # Runninng the Model and Checking
# set.seed(1)
# mod <- stoch_seir_dust$new(# Epidemiological Parameters
#   beta = 0.7, gamma = 1/4, sigma = 1/7, population_size = population_size, start_infections = 10,
#   
#   # Flight Parameters
#   capacity_per_flight = capacity_per_flight, num_flights = num_flights, num_flightsAB = num_flightsAB, 
#   
#   # Sequencing Parameters
#   shedding_freq = 1, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = 1, 
#   seq_tot = 10^9, 
#   samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
#   samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = dt)
# 
# output <- mod$run(1:(end/dt))
# output2 <- mod$transform_variables(output)
# 
# sum(output2$n_inf_flightABOut)
# sum(output2$n_inf_flightOut)
# sum(output2$n_inf_flightOut)/sum(output2$n_inf_flightABOut)
# sum(output2$n_inf_flightOut)/max(output2$R)
# (1/mod$contents()$sigma) * num_flights * capacity_per_flight / population_size
# 
# plot(output2$n_inf_flightOut/output2$n_inf_flightABOut, output2$n_inf_flightOut)
# x <- output2$n_inf_flightOut/output2$n_inf_flightABOut
# mean(x[!is.infinite(x) & !is.nan(x)])
# num_flights/num_flightsAB
# plot(output2$n_inf_flightOut/output2$n_inf_flightABOut, output2$n_inf_flightABOut)
# metrics <- unique(beta_df_summary$metric)
# letters <- c("a", "b", "c", "d", "d", "f", "g", "h", "j", "k")
# plot_list <- vector(mode = "list", length = length(letters))
# for (i in 1:length(metrics)) {
#   temp <- metrics[i]
#   new_df <- beta_df_summary %>%
#     filter(metric == temp) 
#   temp_plot <- ggplot(data = new_df) +
#     geom_line(aes(x = R0, y = avg)) +
#     geom_line(aes(x = R0, y = scaled)) +
#     geom_ribbon(aes(x = R0, y = avg, ymin = lower, ymax = upper), alpha = 0.2) +
#     lims(x = c(min(R0_df$R0), max(R0_df$R0))) + 
#     scale_y_continuous(limits = c(0, NA), 
#                        sec.axis = sec_axis(~./scaling_factor,
#                                            #~./unique(new_df$scaling_factor),
#                                            breaks = c(0, 20, 40, 60, 80, 100),
#                                            name = "% Simulations With\nSuccessful Detection")) +
#     labs(x = "Parameter Value - R0", y = "Time to Detection") +
#     theme_bw()
#   eval(parse(text = paste0(letters[i], " <- temp_plot")))
#   plot_list[[i]] <- temp_plot 
#   print(i)
# }

# ggplot(data = beta_df_summary) +
#   geom_line(aes(x = R0, y = avg_ttd)) +
#   geom_line(aes(x = R0, y = scaled)) +
#   geom_ribbon(aes(x = R0, y = avg_ttd, ymin = lower_ttd, ymax = upper_ttd), alpha = 0.2) +
#   facet_wrap(~method_calc, scales = "free_y") +
#   lims(x = c(min(R0_df$R0), max(R0_df$R0))) + 
#   scale_y_continuous(limits = c(0, NA), sec.axis=~.) +
#   geom_text(aes(x = R0, y = 0, label = num_reached, family = num_reached)) +
#   labs(x = "Parameter Value - R0", y = "Time to Detection") +
#   theme_bw()

# mod <- stoch_seir_dust$new(# Epidemiological Parameters
#   beta = 0.6, gamma = gamma, sigma = sigma, population_size = population_size, start_infections = start_infections,
#   
#   # Flight Parameters
#   capacity_per_flight = capacity_per_flight, num_flights = num_flights, num_flightsAB = num_flightsAB, 
#   
#   # Sequencing Parameters
#   shedding_freq = shedding_freq, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = met_bias, 
#   seq_tot = seq_tot, 
#   samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
#   samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
#   
#   # Miscellaenous Parameters
#   dt = dt)
# 
# output <- mod$run(1:(end/dt))
# output2 <- mod$transform_variables(output)
# 
# # Plotting the Number of Mapped Reads
# 
# ## Wrangling data into right format
# indiv_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_indivFlight_Out)
# colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flightsAB))
# airplane_reads_df <- indiv_flight_reads %>%
#   pivot_longer(-time, names_to = "flight_num", values_to = "reads") %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(flight_num, time2) %>%
#   summarise(daily_flight_reads = sum(reads)) %>%
#   mutate(time3 = midpoints(time2))
# 
# agg_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_aggFlight_Out)
# colnames(agg_flight_reads) <- c("time", "reads")
# agg_reads_df <- agg_flight_reads %>%
#   mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
#   group_by(time2) %>%
#   summarise(daily_reads = sum(reads)) %>%
#   mutate(time3 = midpoints(time2))
# 
# ggplot() +
#   geom_boxplot(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3),
#                fill = NA, outlier.colour = NA, alpha = 0.1, col = "grey") +
#   geom_jitter(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3),
#               size = 0.5, width = 0.15, alpha = 0.5, height = 0, col = "grey") +
#   geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
#   theme(legend.position = "none") +
#   lims(x = c(0, end/2), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
#   labs(x = "Time (Days)", y = "Number of Reads In Each Flight's Wastewater")
# 
# ggplot() +
#   geom_line(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, colour = flight_num),
#             fill = NA, outlier.colour = NA) +
#   geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
#   theme(legend.position = "none") +
#   lims(x = c(0, end/2), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
#   labs(x = "Time (Days)", y = "Number of Reads In Flight Wastewater")

# 
# sum(output2$n_SE_Output[1:75])/population_size
# sum(output2$n_SE_Output[1:(60/dt)])/population_size