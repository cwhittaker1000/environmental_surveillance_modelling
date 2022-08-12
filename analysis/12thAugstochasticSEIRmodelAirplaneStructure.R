# If running for the first time, run this line
#install.packages("drat"); 
#drat:::add("ncov-ic")
#drat:::add("mrc-ide")
#install.packages("odin")
#install.packages("mcstate")
#install.packages("remotes")
#install.packages("odin")
#remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE)
#pkgbuild::check_build_tools()

# Load required libraries & source helper function
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel)
source("functions/helper_functions.R")

# Changing memory limit to handle large outputs
memory.limit(size = 30000)

# Stochastic SEIR
stoch_seir_dust <- odin::odin({
  
  ################## Model Parameters ##################
  
  ## Epidemiological Parameters
  beta <- user()              # probability of a contact successfully transmitting the disease
  gamma <- user()             # rate of transition from Exposed -> Infectious (incubation period)
  sigma <- user()             # rate of transition from Infectious -> Recovered (rate of recovery)
  population_size <- user()   # overall size of population
  start_infections <- user()  # starting number of infections (in Exposed compartment)
  
  ## NAO/Surveillance Parameters
  capacity_per_flight <- user()    # Capacity of a single flight - note does this need to be adjusted for time step e.g. capacity/timestep = actual capacity where timestep <1
  num_flights <- user()            # Number of flights per day
  num_flightsAB <- user()          # Number of flights per day from Location A to NAO Location (Location B)
  p_detection <- user()            # probability of detecting the virus the infected person shed into the system
  
  ########### Rate -> Probability Conversion ############
  
  ## Converting Epidemiological Rates to Probabilities of Leaving Each Compartment
  lambda <- ((beta * I) / N) 
  p_SE <- 1 - exp(-lambda* dt) # S to E 
  p_EI <- 1 - exp(-gamma * dt) # E to I 
  p_IR <- 1 - exp(-sigma * dt) # I to R
  
  ############# Stochastic Model Updates ###############
  
  ## Epidemiological Parameters
  
  ### Binomial Draws for Number of People Moving Between Compartments
  n_SE <- rbinom(S, p_SE)    # number of individuals infected at each timestep (S->E)
  n_EI <- rbinom(E, p_EI)    # number of individuals becoming infectious at each timestep (E->I)
  n_IR <- rbinom(I, p_IR)    # number of individuals recovering at each timestep (I->R)
  
  ### Stochastic Model Updates for Epidemiological Parameters
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR - n_inf_flight
  update(R) <- R + n_IR
  update(N) <- S + E + I + R
  update(new_exposed) <- n_SE
  update(new_infectious) <- n_EI
  
  initial(n_EI_out) <- 0
  update(n_EI_out) <- n_EI
  initial(n_IR_out) <- 0
  update(n_IR_out) <- n_IR
  
  ## NAO/Surveillance Parameters
  
  # Calculating the number of infected people taking a flight
  p_flight <- if ((num_flights * capacity_per_flight) > population_size) 1 else (num_flights * capacity_per_flight)/population_size
  n_inf_flight <- rbinom(I, dt * p_flight) # number of infected people (E+I) taking flight
  prop_flightsAB <- num_flightsAB/num_flights
  n_inf_flightAB <- rbinom(n_inf_flight, prop_flightsAB) # number of infected people (E+I) taking flight
  
  # Distributing (detectable) infections that take a flight (n_inf_take_flight) across all the possible flights they could get on
  capacity_individual_flight[] <- capacity_per_flight
  n_inf_specific_flightAB[] <- rmhyper(n_inf_flightAB, capacity_individual_flight) # be careful about when n_inf_take_flight gets >>> capacity specific flight
  n_inf_specific_flight_detect[] <- rbinom(n_inf_specific_flightAB[i], p_detection)

  # Tallying the Total Number of Infected (and detected) Individuals Flying from A to B (NAO)
  total_inf_NAO_site <- sum(n_inf_specific_flightAB)
  total_inf_NAO_site_detect <- sum(n_inf_specific_flight_detect)
  
  ### Stochastic Model Updates for Epidemiological Parameters
  update(n_inf_flightOut) <- n_inf_flight
  update(n_inf_flightABOut) <- n_inf_flightAB
  update(n_inf_specific_flightABOut[]) <- n_inf_specific_flightAB[i]
  update(n_inf_specific_flight_detectOut[]) <- n_inf_specific_flight_detect[i]
  update(n_inf_NAOOut) <- total_inf_NAO_site
  update(n_inf_NAO_detectOut) <- total_inf_NAO_site_detect
  
  ## Metagenomic and Sequencing Parameters
  
  virus_per_shed <- user()
  non_virus_per_shed <- user()
  shedding_freq <- user()
  volume_wastewater <- user()
  bias <- user()
  m_tot <- user()
  
  amount_virus_indiv_flight[] <- n_inf_specific_flightAB[i] * virus_per_shed * shedding_freq
  amount_non_virus_indiv_flight[] <- capacity_individual_flight[i] * non_virus_per_shed * shedding_freq
  
  conc_virus_indiv_flight[]  <- amount_virus_indiv_flight[i] / volume_wastewater
  conc_non_virus_indiv_flight[]  <- amount_non_virus_indiv_flight[i] / volume_wastewater
  
  m_i[] <- m_tot * (conc_virus_indiv_flight[i] * bias)/((conc_virus_indiv_flight[i] * bias) + conc_non_virus_indiv_flight[i])
  stoch_m_i_pois[] <- rpois(m_i[i]) # ask Rich about exact parameterisation of nbinom when you ask him create it
  
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
  update(stoch_m_i_pois_Out[]) <- stoch_m_i_pois[i] # replace with NegBinom - ask Rich to implement and for details on how to do mean/dispersion parameterisation
  
  ################# Initial Values ##################
  
  ## Epidemiological Parameters
  initial(N) <- S + E + I + R                       # total population size
  initial(S) <- population_size - start_infections  
  initial(E) <- start_infections                    # starting number of infections
  initial(I) <- 0
  initial(R) <- 0
  initial(new_exposed) <- 0
  initial(new_infectious) <- 0
  
  ## NAO/Surveillance Parameters
  initial(n_inf_flightOut) <- 0
  initial(n_inf_flightABOut) <- 0
  initial(n_inf_specific_flightABOut[]) <- 0
  initial(n_inf_specific_flight_detectOut[]) <- 0
  initial(n_inf_NAOOut) <- 0
  initial(n_inf_NAO_detectOut) <- 0
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # Specifying the dimensions of the vectors
  dim(capacity_individual_flight) <- num_flightsAB
  dim(n_inf_specific_flightAB) <- num_flightsAB
  dim(n_inf_specific_flight_detect) <- num_flightsAB
  dim(n_inf_specific_flightABOut) <- num_flightsAB
  dim(n_inf_specific_flight_detectOut) <- num_flightsAB
  
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


# Specifying the Model
dt <- 0.05
end <- 50
mod <- stoch_seir_dust$new(beta = 2, gamma = 1, sigma = 1, population_size = 10^6, start_infections = 5,
                           capacity_per_flight = 100, num_flights = 200, num_flightsAB = 150, 
                           p_detection = 0.5, dt = dt,
                           virus_per_shed = 10, non_virus_per_shed = 1000, volume_wastewater = 50000,
                           shedding_freq = 1, bias = 1,
                           m_tot = 1000000)
output <- mod$run(1:(end/dt))

initial(amount_virus_indiv_flight_Out[]) <- 0
initial(amount_non_virus_indiv_flight_Out[]) <- 0
initial(conc_virus_indiv_flight_Out[]) <- 0
initial(conc_non_virus_indiv_flight_Out[]) <- 0
initial(m_i_Out[]) <- 0
initial(stoch_m_i_pois_Out[]) <- 0

plot(output[, "time"], output[, "n_inf_specific_flightABOut[1]"], type = "l")
plot(output[, "time"], output[, "amount_virus_indiv_flight_Out[1]"], type = "l")
plot(output[, "time"], output[, "conc_virus_indiv_flight_Out[1]"], type = "l")
plot(output[, "time"], output[, "m_i_Out[1]"], type = "l")
plot(output[, "time"], output[, "stoch_m_i_pois_Out[1]"], type = "l")

# Plotting Overall Epidemic Dynamics
plot(output[, "time"], output[, "S"], ylim = c(0, 10^6), type = "l")
lines(output[, "time"], output[, "I"], ylim = c(0, 10^6), col = "blue")
lines(output[, "time"], output[, "R"], ylim = c(0, 10^6), col = "red")

colnames(output)

sum(output[, "n_inf_flightOut"])/sum(output[, "new_infectious"])

plot(output[, "time"], output[, "m_i_Out[1]"], type = "l")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(25)

for (i in 1:25) {
  if (i == 1) {
    plot(get_daily_outputs(end, dt, paste0("n_inf_specific_flight_detectOut.", i, "."), data.frame(output)), type = "l", col = adjustcolor(cols[i], alpha.f = 0.5), ylim = c(0, 20))
  } else {
    lines(get_daily_outputs(end, dt, paste0("n_inf_specific_flight_detectOut.", i, "."), data.frame(output)), type = "l", col = adjustcolor(cols[i], alpha.f = 0.5))
  }
}

index <- grep("time|n_inf_specific_flight_detectOut*", colnames(output))
indiv_airplanes <- data.frame(output[, index])
colnames(indiv_airplanes) <- gsub("n_inf_specific_flight_detectOut.", "flight", colnames(indiv_airplanes))
colnames(indiv_airplanes) <- gsub("\\.", "", colnames(indiv_airplanes))




airplanes_df <- indiv_airplanes %>%
  pivot_longer(-time, names_to = "flight_num", values_to = "infections") %>%
  mutate(time2 = cut(time, breaks = max(time))) %>%
  group_by(flight_num, time2) %>%
  summarise(daily_flight_infections = sum(infections)) %>%
  mutate(time3 = midpoints(time2))

ggplot(airplanes_df, aes(x = time3, y = daily_flight_infections, col = flight_num)) +
  geom_line(alpha = 0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time (Days)", y = "Daily Number of Detectable Infections On Flight Route")

time_to_detection <- airplanes_df %>%
  group_by(flight_num) %>%
  summarise(time_to_first_inf = time3[min(which(daily_flight_infections != 0))],
            time_to_first_inf2 = time2[min(which(daily_flight_infections != 0))])

set.seed(4200)
ggplot(time_to_detection, aes(x = 1, y = time_to_first_inf)) +
  geom_boxplot(col = "dark blue", outlier.shape = NA) +
  geom_jitter(width = 0.2, 
              colour="black", pch = 21, size = 4, fill = "dark blue") +
  scale_y_continuous(breaks = seq(16, 28, 2)) +
  labs(x = "", y = "Time to First Detectable Infection For Each Flight") +
  theme_bw() +
  coord_flip()

