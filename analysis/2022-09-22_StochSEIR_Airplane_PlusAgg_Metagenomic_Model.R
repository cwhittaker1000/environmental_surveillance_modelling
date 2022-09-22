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
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel)
source("functions/helper_functions.R")

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
  update(I) <- I + n_EI - n_IR - n_inf_flight
  update(R) <- R + n_IR         # NOTE: need to eventually return n_inf_flight to here - still TODO
  update(N) <- S + E + I + R
  update(new_infectious) <- n_EI
  update(n_EI_Output) <- n_EI   # odin.dust doesn't let you output calculated quantities like n_EI & n_IR without making 
  update(n_IR_Output) <- n_IR   # them tracked variables, which requires having initial() and update() calls for them

  ## Initial Values for Epidemiological States & Quantities of Interest
  initial(S) <- population_size - start_infections  
  initial(E) <- start_infections                    
  initial(I) <- 0
  initial(R) <- 0
  initial(N) <- S + E + I + R                       
  initial(new_infectious) <- 0
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
  n_inf_flightAB <- rhyper(n_inf_flight, num_flights * capacity_per_flight * dt, num_flightsAB * capacity_per_flight * dt) # Still figuring out whether dt is required here: 2nd arg divided by 3rd arg is ~equivalent to p for binomial, so doesn't make any huge practical diff (especially when 1st arg << 2nd & 3rd args)

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

# Specifying and Running the Model

# General Parameters
dt <- 0.05
end <- 250
num_flightsAB <- 25
virus_shed <- 1000000
virus_ratio <- 10^7  # assume 10^7 non-viral material shed for each virus shed
vol_flight_ww <- 800 # 800 litres assumed
sample_flight_ww <- 1 # assume we take 1 litre

# Functions to Extract Relevant Quantities
num_reads <- 5
agg_or_indiv <- "aggregated"

time_to_detection <- function(mod_output, num, inf_or_reads, agg_or_indiv) {

  # Checking inputs
  if (!(inf_or_reads %in% c("infections", "reads"))) {
    stop("agg_or_indiv must be either infections or reads")
  } 
  if (!(agg_or_indiv %in% c("individual", "aggregated"))) {
    stop("agg_or_indiv must be either individual or aggregated")
  } 
  
  # For aggregated flight outputs # NEED TO CONVERT PER TIMEPOINT TO PER DAY - IMPORTANT!!!!
  if (agg_or_indiv == "aggregated") {
    if (inf_or_reads == "infections") {
      index_greater_than_num_for_total <- which(mod_output$n_inf_flightABOut >= num)
      index_greater_than_num_for_average <- which(mod_output$n_inf_flightABOut/dim(mod_output$n_inf_specific_flightABOut)[2] >= num)
      if (identical(index_greater_than_num, integer(0))) {
        time <- NA
      } else {
        time_for_total <- mod_output$time[min(index_greater_than_num)]
      }
    } else if (inf_or_reads == "reads") {
      index_greater_than_num <- which(mod_output$seq_reads_virus_aggFlight_Out >= num)
      if (identical(index_greater_than_num, integer(0))) {
        time <- NA
      } else {
        time <- mod_output$time[min(index_greater_than_num)]
      }      
    }
  }
  
  # For individual flight outputs
  if (agg_or_indiv == "individual") {
    if (inf_or_reads == "infections") {
      indiv_flight_infections <- data.frame(time = mod_output$time, mod_output$n_inf_specific_flightABOut)
      colnames(indiv_flight_infections) <- c("time", paste0("flight", 1:num_flights))
      airplane_infections_df <- indiv_flight_infections %>%
        pivot_longer(-time, names_to = "flight_num", values_to = "infections") %>%
        mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
        group_by(flight_num, time2) %>%
        summarise(daily_flight_infections = sum(infections)) %>%
        mutate(time3 = midpoints(time2)) %>%
        group_by(flight_num) %>%
        summarise(time_num = time3[min(which(daily_flight_infections > num))])
      time <- unname(airplane_infections_df$time_num)
    } else if (inf_or_reads == "reads") {
      indiv_flight_reads <- data.frame(time = mod_output$time, mod_output$seq_reads_virus_indivFlight_Out)
      colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))
      airplane_reads_df <- indiv_flight_reads %>%
        pivot_longer(-time, names_to = "flight_num", values_to = "reads") %>%
        mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
        group_by(flight_num, time2) %>%
        summarise(daily_flight_reads = sum(reads)) %>%
        mutate(time3 = midpoints(time2)) %>%
        group_by(flight_num) %>%
        summarise(time_num = time3[min(which(daily_flight_reads > num))])
      time <- unname(airplane_reads_df$time_num)
    }
  }
  return(time)
}

time_to_detection(output2, 5, "infections", "aggregated")
time_to_detection(output2, 5, "infections", "individual")
time_to_detection(output2, 5, "reads", "aggregated")
time_to_detection(output2, 5, "reads", "individual")


# Looping Over Beta Parameters 
beta <- seq(0.3, 3, 0.1)
for (i in 1:length(beta)) {
  
  # Running the Model
  set.seed(2)
  mod <- stoch_seir_dust$new(# Epidemiological Parameters
    beta = beta[i], gamma = 1/4, sigma = 1/5, population_size = 10^6, start_infections = 10,
    
    # Flight Parameters
    capacity_per_flight = 250, num_flights = 50, num_flightsAB = num_flightsAB, 
    
    # Sequencing Parameters
    shedding_freq = 1, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = 1, 
    seq_tot = 10^9, 
    samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
    samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
    
    # Miscellaenous Parameters
    dt = dt)
  
  # Extracting Output
  output <- mod$run(1:(end/dt))
  output2 <- mod$transform_variables(output)
  
}


set.seed(2)
mod <- stoch_seir_dust$new(# Epidemiological Parameters
                           beta = 0.6, gamma = 1/4, sigma = 1/5, population_size = 10^6, start_infections = 10,
                           
                           # Flight Parameters
                           capacity_per_flight = 250, num_flights = 50, num_flightsAB = num_flightsAB, 
                           
                           # Sequencing Parameters
                           shedding_freq = 1, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = 1, 
                           seq_tot = 10^9, 
                           samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
                           samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
                           
                           # Miscellaenous Parameters
                           dt = dt)
output <- mod$run(1:(end/dt))
output2 <- mod$transform_variables(output)

# Plotting Overall Epidemic Dynamics
state_vars_df <- data.frame(time = output2$time, S = output2$S, E = output2$E, I = output2$I, R = output2$R) %>%
  pivot_longer(-time, names_to = "state_var", values_to = "value") 
ggplot(state_vars_df, aes(x = time, y = value, colour = state_var)) +
  geom_line() +
  labs(x = "Time (Days)", y = "Number of Individuals in Each State")
  
# Plotting the Number of Infections on Each Flight Per Day

## Wrangling data into right format
num_flights <- mod$contents()[["num_flightsAB"]]
indiv_flight_infections <- data.frame(time = output2$time, output2$n_inf_specific_flightABOut)
colnames(indiv_flight_infections) <- c("time", paste0("flight", 1:num_flights))
airplane_infections_df <- indiv_flight_infections %>%
  pivot_longer(-time, names_to = "flight_num", values_to = "infections") %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(flight_num, time2) %>%
  summarise(daily_flight_infections = sum(infections)) %>%
  mutate(time3 = midpoints(time2))

### Boxplot 
ggplot(airplane_infections_df, aes(x = time3, y = daily_flight_infections, group = time3)) +
  geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.5, height = 0) +
  theme(legend.position = "none") +
  lims(x = c(0, end), y = c(0, max(airplane_infections_df$daily_flight_infections))) +
  labs(x = "Time (Days)", y = "Number of Infections On Each Daily Flight")

### Line
ggplot(airplane_infections_df, aes(x = time3, y = daily_flight_infections, colour = flight_num)) +
  geom_line(fill = NA, outlier.colour = NA) +
  theme(legend.position = "none") +
  lims(x = c(0, end), y = c(0, max(airplane_infections_df$daily_flight_infections))) +
  labs(x = "Time (Days)", y = "Number of Infections On Each Daily Flight")

# Plotting the Number of Mapped Reads 

## Wrangling data into right format
indiv_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_indivFlight_Out)
colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))
airplane_reads_df <- indiv_flight_reads %>%
  pivot_longer(-time, names_to = "flight_num", values_to = "reads") %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(flight_num, time2) %>%
  summarise(daily_flight_reads = sum(reads)) %>%
  mutate(time3 = midpoints(time2))  

agg_flight_reads <- data.frame(time = output2$time, output2$seq_reads_virus_aggFlight_Out)
colnames(agg_flight_reads) <- c("time", "reads")
agg_reads_df <- agg_flight_reads %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(time2) %>%
  summarise(daily_reads = sum(reads)) %>%
  mutate(time3 = midpoints(time2))  

#### Boxplot
ggplot() +
  geom_boxplot(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3), 
               fill = NA, outlier.colour = NA, alpha = 0.1, col = "grey") +
  geom_jitter(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3), 
              size = 0.5, width = 0.15, alpha = 0.5, height = 0, col = "grey") +
  geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
  theme(legend.position = "none") +
  lims(x = c(0, end/2), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
  labs(x = "Time (Days)", y = "Number of Reads In Each Flight's Wastewater")

### Line
ggplot() +
  geom_line(data = airplane_reads_df, aes(x = time3, y = daily_flight_reads, colour = flight_num), 
            fill = NA, outlier.colour = NA) +
  geom_line(data = agg_reads_df, aes(x = time3, y = daily_reads)) +
  theme(legend.position = "none") +
  lims(x = c(0, 45), y = c(0, max(airplane_reads_df$daily_flight_reads))) +
  labs(x = "Time (Days)", y = "Number of Reads In Flight Wastewater")

