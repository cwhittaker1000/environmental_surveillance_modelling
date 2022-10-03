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
  shedding_prop <- user()              # Proportion of infected people who are shedding
  virus_shed <- user()                 # Average amount of material shed per event for our virus of interest (defecation)
  non_virus_shed <- user()             # Average amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  met_bias <- user()                   # Bias term for the metagenomic model
  seq_tot <- user()                    # Total amount of sequencing that is done
  samp_frac_aggFlight <- user()      # Fraction of all flight's wastewater that would be sampled if individual flight wastewater pooled and then sampled (via triturator) 
  
  ## AGGREGATED FLIGHT CALCULATIONS (E.G. SAMPLING FROM TRITURATOR)
  
  # Calculating the number of shedding events from infected and uninfected individuals on each airplane
  infected_indiv_shedding_events <- rpois(n_inf_flightAB * shedding_freq * shedding_prop)
  uninfected_indivs_flightsAB <- (capacity_per_flight * num_flightsAB) - (n_inf_flightAB * shedding_prop)
  uninfected_indiv_shedding_events <- rpois(uninfected_indivs_flightsAB * shedding_freq) 
  
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
  update(infected_indiv_shedding_events_Out) <- infected_indiv_shedding_events
  update(uninfected_indiv_shedding_events_Out) <- uninfected_indiv_shedding_events
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
  initial(infected_indiv_shedding_events_Out) <- 0
  initial(uninfected_indiv_shedding_events_Out) <- 0
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


#### even more det tester
# infected_indiv_shedding_events_det <- n_inf_flightAB * shedding_freq * shedding_prop
# uninfected_indivs_flightsAB_det <- (capacity_per_flight * num_flightsAB) - (n_inf_flightAB * shedding_prop)
# uninfected_indiv_shedding_events_det <- uninfected_indivs_flightsAB * shedding_freq 
# amount_virus_aggFlight_det <- infected_indiv_shedding_events_det * virus_shed 
# amount_non_virus_aggFlight_det <- (uninfected_indiv_shedding_events_det + infected_indiv_shedding_events_det) * non_virus_shed
# 
# sample_amount_virus_aggFlight_detdet <- amount_virus_aggFlight_det * samp_frac_aggFlight
# sample_amount_non_virus_aggFlight_detdet <- amount_non_virus_aggFlight_det * samp_frac_aggFlight
# seq_reads_virus_aggFlight_detdet <- if(sample_amount_virus_aggFlight_detdet == 0 && sample_amount_non_virus_aggFlight_detdet == 0) 0 else seq_tot * (sample_amount_virus_aggFlight_detdet * met_bias)/((sample_amount_virus_aggFlight_detdet * met_bias) + sample_amount_non_virus_aggFlight_detdet)
# seq_reads_non_virus_aggFlight_detdet <- seq_tot - seq_reads_virus_aggFlight_det
# 
# update(infected_indiv_shedding_events_det_Out) <- infected_indiv_shedding_events_det
# update(uninfected_indiv_shedding_events_det_Out) <- uninfected_indiv_shedding_events_det
# update(amount_virus_aggFlight_det_Out) <- amount_virus_aggFlight_det
# update(amount_non_virus_aggFlight_det_Out) <- amount_non_virus_aggFlight_det
# update(sample_amount_virus_aggFlight_detdet_Out) <- sample_amount_virus_aggFlight_detdet
# update(sample_amount_non_virus_aggFlight_detdet_Out) <- sample_amount_non_virus_aggFlight_detdet
# update(seq_reads_virus_aggFlight_detdet_Out) <- seq_reads_virus_aggFlight_detdet
# update(seq_reads_non_virus_aggFlight_detdet_Out) <- seq_reads_non_virus_aggFlight_detdet
# 
# initial(infected_indiv_shedding_events_det_Out) <- 0
# initial(uninfected_indiv_shedding_events_det_Out) <- 0
# initial(amount_virus_aggFlight_det_Out) <- 0
# initial(amount_non_virus_aggFlight_det_Out) <- 0
# initial(sample_amount_virus_aggFlight_detdet_Out) <- 0
# initial(sample_amount_non_virus_aggFlight_detdet_Out) <- 0
# initial(seq_reads_virus_aggFlight_detdet_Out) <- 0
# initial(seq_reads_non_virus_aggFlight_detdet_Out) <- 0