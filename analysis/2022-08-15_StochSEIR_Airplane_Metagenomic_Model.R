# Installing Required Libraries (If Required)
install_packages <- FALSE 
if(install_packages) { 
  install.packages("drat")
  drat:::add("mrc-ide")
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
  update(R) <- R + n_IR # NOTE: need to eventually return n_inf_flight to here 
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
  capacity_per_flight <- user()    # Capacity of a single flight - note does this need to be adjusted for time step e.g. capacity/timestep = actual capacity where timestep <1
  num_flights <- user()            # Number of flights per day
  num_flightsAB <- user()          # Number of flights per day from Location A to NAO Location (Location B)

  ## Calculating the number of infected people taking a flight based on number of infections and airport/airplane parameters
  ### Note for rhyper calls: 1st param = # white balls; 2nd param # black balls; 3rd param # balls drawn; the draw is the number of white balls that get drawn in this context
  n_inf_flight <- rhyper(I, population_size - I, dt * num_flights * capacity_per_flight) 
  n_inf_flightAB <- rhyper(n_inf_flight, num_flights * capacity_per_flight * dt, num_flightsAB * capacity_per_flight * dt) # not sure whether dt is required here, but 2nd arg/3rd arg is effectively p for binomial, so doesn't make any practical diff (I don't think)
  
  ### Note: Previously used rbinom calls in place of hypergeometric - probably fine as binomial approxs hypergeometric well
  ###       when number successes is << total population size. Keeping them here in case we decide to revert 
  #p_flight <- if ((num_flights * capacity_per_flight) > population_size) 1 else (num_flights * capacity_per_flight)/population_size
  #prop_flightsAB <- num_flightsAB/num_flights
  #n_inf_flight <- rbinom(I, dt * p_flight) # number of infected people taking flight - binomial approx of hypergeometric is fine here
  #n_inf_flightAB <- rbinom(n_inf_flight, prop_flightsAB) # number of infected people taking flight
  
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

  ################## Metagenomic Sequencing Model ##################
  
  ## Metagenomic and Sequencing Parameters
  virus_shed <- user()          # Amount of material shed per event for our virus of interest (defecation)
  non_virus_shed <- user()      # Amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  shedding_freq <- user()       # Rate at which individuals defecate on a flight
  volume_wastewater <- user()   # Volume of wastewater into which individuals shed (to convert abundance to concentration)
  bias <- user()                # Bas term for the metagenomic model
  m_tot <- user()               # Total amount of sequencing that is done
  
  # Calculating the amount of nucleic acid shed into wastewater on each flight (i.e. abundance)
  amount_virus_indiv_flight[] <- n_inf_specific_flightAB[i] * virus_shed * shedding_freq
  amount_non_virus_indiv_flight[] <- capacity_individual_flight[i] * non_virus_shed * shedding_freq
  
  # Converting the abundance of nucleic acid on each flight into the concentration
  conc_virus_indiv_flight[]  <- amount_virus_indiv_flight[i] / volume_wastewater
  conc_non_virus_indiv_flight[]  <- amount_non_virus_indiv_flight[i] / volume_wastewater
  
  # Converting this nucleic acid concentration into sequencing reads
  m_i[] <- m_tot * (conc_virus_indiv_flight[i] * bias)/((conc_virus_indiv_flight[i] * bias) + conc_non_virus_indiv_flight[i])
  stoch_m_i_pois[] <- rpois(m_i[i]) # NBinom to be implemented shortly
  
  ### Stochastic Model Updates for Outputs Relevant to Metagenomic Sequencing
  
  ### Initial Values for Outputs Relevant to Metagenomic Sequencing
  initial(amount_virus_indiv_flight_Out[]) <- 0
  initial(amount_non_virus_indiv_flight_Out[]) <- 0
  initial(conc_virus_indiv_flight_Out[]) <- 0
  initial(conc_non_virus_indiv_flight_Out[]) <- 0
  initial(m_i_Out[]) <- 0
  initial(stoch_m_i_pois_Out[]) <- 0
  
  update(amount_virus_indiv_flight_Out[]) <- amount_virus_indiv_flight[i]
  update(amount_non_virus_indiv_flight_Out[]) <- amount_non_virus_indiv_flight[i]
  update(conc_virus_indiv_flight_Out[]) <- conc_virus_indiv_flight[i]
  update(conc_non_virus_indiv_flight_Out[]) <- conc_non_virus_indiv_flight[i]
  update(m_i_Out[]) <- m_i[i]
  update(stoch_m_i_pois_Out[]) <- stoch_m_i_pois[i] # replace with NegBinom when implemented
  
  ################# Miscellaneous Model Stuff ##################
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # Specifying the dimensions of the different vectors in the model (required for C compilation)
  dim(capacity_individual_flight) <- num_flightsAB
  dim(n_inf_specific_flightAB) <- num_flightsAB
  dim(n_inf_specific_flightABOut) <- num_flightsAB
  dim(amount_virus_indiv_flight) <- num_flightsAB
  dim(amount_non_virus_indiv_flight) <- num_flightsAB
  dim(conc_virus_indiv_flight) <- num_flightsAB
  dim(conc_non_virus_indiv_flight) <- num_flightsAB
  dim(m_i) <- num_flightsAB
  dim(stoch_m_i_pois) <- num_flightsAB
  dim(amount_virus_indiv_flight_Out) <- num_flightsAB
  dim(amount_non_virus_indiv_flight_Out) <- num_flightsAB
  dim(conc_virus_indiv_flight_Out) <- num_flightsAB
  dim(conc_non_virus_indiv_flight_Out) <- num_flightsAB
  dim(m_i_Out) <- num_flightsAB
  dim(stoch_m_i_pois_Out) <- num_flightsAB

})

# Specifying and Running the Model
dt <- 0.05
end <- 50
mod <- stoch_seir_dust$new(beta = 2, gamma = 1, sigma = 1, population_size = 10^6, start_infections = 50,
                           capacity_per_flight = 2000, num_flights = 50, num_flightsAB = 25, 
                           dt = dt,
                           virus_shed = 10, non_virus_shed = 1000, volume_wastewater = 50000,
                           shedding_freq = 1, bias = 1,
                           m_tot = 1000000)
output <- mod$run(1:(end/dt))
output2 <- mod$transform_variables(output)

# Plotting Overall Epidemic Dynamics
state_vars_df <- data.frame(time = output2$time, S = output2$S, E = output2$E, I = output2$I, R = output2$R) %>%
  pivot_longer(-time, names_to = "state_var", values_to = "value") 
ggplot(state_vars_df, aes(x = time, y = value, colour = state_var)) +
  geom_line() +
  labs(x = "Time (Days)", y = "Number of Individuals in Each State")
  
# Plotting the Number of Infections on Each Flight Per Day
num_flights <- mod$contents()[["num_flightsAB"]]
indiv_flight_infections <- data.frame(time = output2$time, output2$n_inf_specific_flightABOut)
colnames(indiv_flight_infections) <- c("time", paste0("flight", 1:num_flights))

airplane_infections_df <- indiv_flight_infections %>%
  pivot_longer(-time, names_to = "flight_num", values_to = "infections") %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(flight_num, time2) %>%
  summarise(daily_flight_infections = sum(infections)) %>%
  mutate(time3 = midpoints(time2))

## Boxplot showing the variation in number of infections across the flights on a given day
ggplot(airplane_infections_df, aes(x = time3, y = daily_flight_infections, group = time3)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  lims(x = c(0, 20)) +
  labs(x = "Time (Days)", y = "Daily Number of Detectable Infections On Flight Route")

# Plotting the Number of Mapped Reads on Each Flight Per Day (Deterministic Model)
indiv_flight_reads <- data.frame(time = output2$time, output2$m_i_Out)
colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))

airplane_reads_df <- indiv_flight_reads %>%
  pivot_longer(-time, names_to = "flight_num", values_to = "reads") %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(flight_num, time2) %>%
  summarise(daily_flight_reads = sum(reads)) %>%
  mutate(time3 = midpoints(time2))  

ggplot(airplane_reads_df, aes(x = time3, y = daily_flight_reads, group = time3)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  lims(x = c(0, 20)) +
  labs(x = "Time (Days)", y = "Daily Number of Reads On Flight Route")
