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
dt <- 0.2
end <- 100
num_flights <- 50
num_flightsAB <- 25
virus_shed <- 1000000
virus_ratio <- 10^7  # assume 10^7 non-viral material shed for each virus shed
vol_flight_ww <- 800 # 800 litres assumed
sample_flight_ww <- 1 # assume we take 1 litre

# Functions to Extract Relevant Quantities

# Looping Over Beta Parameters 
options(dplyr.summarise.inform = FALSE)
tic()
num_reads <- 5
stochastic_sim <- 50
beta <- seq(0.3, 2, 0.1)
output_df <- data.frame(beta = rep(beta, each = stochastic_sim), stochastic_realisation = 1:stochastic_sim, num_reads = num_reads, 
                        first_flight_reads = NA, avg_flight_reads = NA, agg_flight_reads = NA, 
                        first_flight_inf = NA, avg_flight_inf = NA)
seed <- rpois(100, lambda = 200) * rpois(100, lambda = 200) * rpois(100, lambda = 20) 
counter <- 1
for (i in 1:length(beta)) {
  
  for (j in 1:stochastic_sim) {
    
    # Running the Model
    set.seed(seed[j])
    mod <- stoch_seir_dust$new(# Epidemiological Parameters
      beta = beta[i], gamma = 1/4, sigma = 1/5, population_size = 10^6, start_infections = 10,
      
      # Flight Parameters
      capacity_per_flight = 250, num_flights = num_flights, num_flightsAB = num_flightsAB, 
      
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
    
    # Summarising Output as Time to Detection
    
    output_df$first_flight_reads[counter] <- min(time_to_detection_fun(output2, num_reads, "reads", "individual"))
    output_df$avg_flight_reads[counter] <- mean(time_to_detection_fun(output2, num_reads, "reads", "individual")) # check whether this is correct or need to find when average line passes 5 for first time
    output_df$agg_flight_reads[counter] <- time_to_detection_fun(output2, num_reads, "reads", "aggregated")
    output_df$first_flight_inf[counter] <- min(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    output_df$avg_flight_inf[counter] <- mean(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    
    counter <- counter + 1 
    #print(j)
  }
  print(i)
}
toc()






## Summarising the Output
x <- output_df %>%
  pivot_longer(-c(beta, num_reads, stochastic_realisation), values_to = "am_time_to_detection", names_to = "metric")
ggplot(x, aes(x = beta, y = am_time_to_detection, group = stochastic_realisation)) +
  geom_line() +
  facet_wrap(~metric)

y <- x %>%
  group_by(beta, metric) %>%
  summarise(avg = mean(am_time_to_detection, na.rm = TRUE),
            lower = min(am_time_to_detection, na.rm = TRUE),
            upper = max(am_time_to_detection, na.rm = TRUE),
            num_not_reached = paste0(sum(is.na(am_time_to_detection))))
  
metric_names <- list(
  "agg_flight_reads"=paste0("Time to >", num_reads, " Reads, Agg. Wastewater"),
  "avg_flight_inf"=paste0("Time to >", num_reads, " Average Infections Per Flight"),
  "avg_flight_reads"=paste0("Time to >", num_reads, " Average Reads Per Flight"),
  "first_flight_inf"=paste0("Time to First Flight With >", num_reads, " Infections"),
  "first_flight_reads"=paste0("Time to First Flight With >", num_reads, " SReads"))
metric_labeller <- function(variable,value){
  return(metric_names[value])
}

ggplot() +
  geom_line(data = y, aes(x = beta, y = avg)) +
  geom_ribbon(data = y, aes(x = beta, y = avg, ymin = lower, ymax = upper), alpha = 0.2) +
  lims(x = c(0, max(beta)), y = c(0, max(y$upper))) +
  facet_wrap(~metric, labeller = metric_labeller) +
  geom_text(data = y, aes(x = beta, y = max(upper), label = num_not_reached, family = num_not_reached)) +
  labs(x = "Parameter Value - Beta", y = "Time to Detection") +
  theme_bw()



# Looping Over Total Sequencing
options(dplyr.summarise.inform = FALSE)
tic()
num_reads <- 5
stochastic_sim <- 50
seq_tot <- round(lseq(10^9, 10^12, 20))
output_df <- data.frame(seq_tot = rep(seq_tot, each = stochastic_sim), stochastic_realisation = 1:stochastic_sim, num_reads = num_reads, 
                        first_flight_reads = NA, avg_flight_reads = NA, agg_flight_reads = NA, 
                        first_flight_inf = NA, avg_flight_inf = NA)
seed <- rpois(100, lambda = 200) * rpois(100, lambda = 200) * rpois(100, lambda = 20) 
counter <- 1
for (i in 1:length(seq_tot)) {
  
  for (j in 1:stochastic_sim) {
    
    # Running the Model
    set.seed(seed[j])
    mod <- stoch_seir_dust$new(# Epidemiological Parameters
      beta = 0.6, gamma = 1/4, sigma = 1/5, population_size = 10^6, start_infections = 10,
      
      # Flight Parameters
      capacity_per_flight = 250, num_flights = num_flights, num_flightsAB = num_flightsAB, 
      
      # Sequencing Parameters
      shedding_freq = 1, virus_shed = virus_shed, non_virus_shed = virus_shed * virus_ratio, met_bias = 1, 
      seq_tot = seq_tot[i], 
      samp_frac_indivFlight = (sample_flight_ww/num_flightsAB)/vol_flight_ww, 
      samp_frac_aggFlight = sample_flight_ww/(vol_flight_ww * num_flightsAB),
      
      # Miscellaenous Parameters
      dt = dt)
    
    # Extracting Output
    output <- mod$run(1:(end/dt))
    output2 <- mod$transform_variables(output)
    
    # Summarising Output as Time to Detection
    
    output_df$first_flight_reads[counter] <- min(time_to_detection_fun(output2, num_reads, "reads", "individual"))
    output_df$avg_flight_reads[counter] <- mean(time_to_detection_fun(output2, num_reads, "reads", "individual")) # check whether this is correct or need to find when average line passes 5 for first time
    output_df$agg_flight_reads[counter] <- time_to_detection_fun(output2, num_reads, "reads", "aggregated")
    output_df$first_flight_inf[counter] <- min(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    output_df$avg_flight_inf[counter] <- mean(time_to_detection_fun(output2, num_reads, "infections", "individual"))
    
    counter <- counter + 1 
    #print(j)
  }
  print(i)
}
toc()

## Summarising the Output
x <- output_df %>%
  pivot_longer(-c(seq_tot, num_reads, stochastic_realisation), values_to = "am_time_to_detection", names_to = "metric")
ggplot(x, aes(x = seq_tot, y = am_time_to_detection, group = stochastic_realisation)) +
  geom_line() +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~metric)

y <- x %>%
  group_by(seq_tot, metric) %>%
  filter(metric %in% c("agg_flight_reads", "avg_flight_reads", "first_flight_reads")) %>%
  summarise(avg = mean(am_time_to_detection, na.rm = TRUE),
            lower = min(am_time_to_detection, na.rm = TRUE),
            upper = max(am_time_to_detection, na.rm = TRUE),
            num_not_reached = paste0(sum(is.na(am_time_to_detection))))

metric_names <- list(
  "agg_flight_reads"=paste0("Time to >", num_reads, " Reads, Agg. Wastewater"),
  "avg_flight_inf"=paste0("Time to >", num_reads, " Average Infections Per Flight"),
  "avg_flight_reads"=paste0("Time to >", num_reads, " Average Reads Per Flight"),
  "first_flight_inf"=paste0("Time to First Flight With >", num_reads, " Infections"),
  "first_flight_reads"=paste0("Time to First Flight With >", num_reads, " SReads"))
metric_labeller <- function(variable,value){
  return(metric_names[value])
}

ggplot() +
  geom_line(data = y, aes(x = seq_tot, y = avg)) +
  geom_ribbon(data = y, aes(x = seq_tot, y = avg, ymin = lower, ymax = upper), alpha = 0.2) +
  lims(x = c(0, max(seq_tot)), y = c(0, max(y$upper))) +
  facet_wrap(~metric, labeller = metric_labeller) +
  geom_text(data = y, aes(x = seq_tot, y = max(upper), label = num_not_reached, family = num_not_reached)) +
  labs(x = "Parameter Value - Total Sequenced", y = "Time to Detection") +
  scale_x_continuous(trans = "log10") +
  theme_bw()


#### OLD AND CHECKING INDIVIDUAL RUNS #### 

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

agg_flight_infs <- data.frame(time = output2$time, output2$n_inf_flightABOut)
colnames(agg_flight_infs) <- c("time", "infections")
agg_inf_df <- agg_flight_infs %>%
  mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
  group_by(time2) %>%
  summarise(daily_infections = sum(infections)/dim(output2$n_inf_specific_flightABOut)[2]) %>%
  mutate(time3 = midpoints(time2))  

### Boxplot 
ggplot(airplane_infections_df, aes(x = time3, y = daily_flight_infections, group = time3)) +
  geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.5, height = 0) +
  theme(legend.position = "none") +
  lims(x = c(0, end), y = c(0, max(airplane_infections_df$daily_flight_infections))) +
  labs(x = "Time (Days)", y = "Number of Infections On Each Daily Flight")

### Line
ggplot() +
  geom_line(data = airplane_infections_df, aes(x = time3, y = daily_flight_infections, colour = flight_num), 
            fill = NA, outlier.colour = NA) +
  geom_line(data = agg_inf_df, aes(x = time3, y = daily_infections), 
            fill = NA) +
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

