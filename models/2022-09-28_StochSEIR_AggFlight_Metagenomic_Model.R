# Stochastic SEIR
stoch_seir <- odin::odin({
  
  ################## Epidemiological Model ##################
  
  ## Epidemiological Parameters
  beta <- user()              # transmission rate - rate of successful infectious contact
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
  update(I) <- I + n_EI - n_IR - n_inf_flight_leaveA 
  update(R) <- R + n_IR         # NOTE: need to eventually return n_inf_flightA to here - still TODO
  update(N) <- S + E + I + R
  update(n_SE_Output) <- n_SE   
  update(n_EI_Output) <- n_EI   
  update(n_IR_Output) <- n_IR
  
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
  num_flights_leaveA <- user()                                    # Total number of flights per day leaving Location A
  num_flights_AtoB <- user()                                      # Number of flights per day from Location A to NAO Location (Location B)
  num_flights_arriveB <- user()                                   # Total number of flights per day arriving at Location B from ALL locations. 
  num_flights_OtoB <- num_flights_arriveB - num_flights_AtoB      # Total number of flights per day arriving at Location B from OTHER Locations (i.e. NOT A)
  capacity_per_flight <- user()                                   # Capacity of single flight (assumed to be same across all flights)
  
  ## Calculating the number of infected people taking a flight based on number of infections and airport/airplane parameters
  ### rhyper takes 3 params: m, n and k. m = # white balls, n = # black balls, k = # balls drawn from urn (without replacement), output of draw is # white balls drawn
  n_inf_flight_leaveA <- rhyper(I, population_size - I, num_flights_leaveA * capacity_per_flight * dt) 
  n_inf_flight_AtoB <- rhyper(n_inf_flight_leaveA, (num_flights_leaveA * capacity_per_flight * dt) - n_inf_flight_leaveA, num_flights_AtoB * capacity_per_flight * dt) 
  
  ### Stochastic Model Updates for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  update(n_inf_flight_leaveA_Out) <- n_inf_flight_leaveA
  update(n_inf_flight_AtoB_Out) <- n_inf_flight_AtoB
  
  ### Initial Values for Outputs Relevant to Surveillance (Number Infected On Flights Etc)
  initial(n_inf_flight_leaveA_Out) <- 0
  initial(n_inf_flight_AtoB_Out) <- 0
  
  ################## Metagenomic Sequencing Model ##################
  ## Note: this aggregates wastewater from all flights arriving in Location B and assumes we're sampling from triturator
  
  ## Metagenomic and Sequencing Parameters
  shedding_freq <- user()              # Average number of defecation events per person per flight
  shedding_prop <- user()              # Proportion of infected people who are shedding
  virus_shed <- user()                 # Average amount of material shed per event for our virus of interest (defecation)
  non_virus_shed <- user()             # Average amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  met_bias <- user()                   # Bias term for the metagenomic model
  seq_tot <- user()                    # Total amount of sequencing that is done
  samp_frac_aggFlight <- user()        # Fraction of all flight's wastewater that would be sampled if individual flight wastewater pooled and then sampled (via triturator) 
  
  # Calculating the number of shedding events from infected and uninfected individuals across flights
  ## Note there are 4 types of individual we're tracking - 1) infected & shedding individuals travelling A->B 
  ##                                                       2) infected and not-shedding individuals travelling A->B
  ##                                                       3) uninfected individuals travelling A->B
  ##                                                       4) uninfected individuals travelling O->B (where O is all other locations apart from A)
  ##                                                      For our purposes, 2) and 3) are the same (not shedding NA of interest, travelling A->B)
  n_inf_AtoB_shedding <- rbinom(n_inf_flight_AtoB, shedding_prop)                               # Calculating number of infected people travelling A->B who actually shed NA of interest
  n_inf_AtoB_shedding_events <- rpois(n_inf_AtoB_shedding * shedding_freq)                      # Multiply number of infected people on A->B flights who are shedding by a rate of shedding (i.e. defectations per flight)
  
  n_uninf_AtoB_shedding <- (capacity_per_flight * num_flights_AtoB * dt) - n_inf_AtoB_shedding  # Calculating number of people travelling A->B who do not shed NA of interest - includes uninfected people and those infected but not shedding
  n_uninf_AtoB_shedding_events <- rpois(n_uninf_AtoB_shedding * shedding_freq)                  # Multiply number of people on A->B flights not shedding NA of interest by a rate of shedding (i.e. defectations per flight)
  
  n_uninf_OtoB_shedding <- (capacity_per_flight * num_flights_OtoB  * dt)                       # Calculating number of people travelling O->B who are not shedding NA of interest
  n_uninf_OtoB_shedding_events <- rpois(n_uninf_OtoB_shedding * shedding_freq)                  # Multiply this by rate of shedding (i.e. defectations per flight)
  
  ### Calculating amount and concentration of nucleic acid shed into aggregated flight wastewater
  amount_virus_aggFlight <- n_inf_AtoB_shedding_events * virus_shed 
  amount_non_virus_aggFlight <- (n_inf_AtoB_shedding_events + n_uninf_AtoB_shedding_events + n_uninf_OtoB_shedding_events) * non_virus_shed
  
  ### Converting this nucleic acid abundance into sequencing reads for aggregated flight wastewater
  sample_amount_virus_aggFlight <- amount_virus_aggFlight * samp_frac_aggFlight
  sample_amount_non_virus_aggFlight <- amount_non_virus_aggFlight * samp_frac_aggFlight
  
  # Note that I don't think this is quite right currently, as this is basically calculating relative abundance and reads produced
  # from each timestep, rather than reads produced from the overall relative abundance across a single day (which results from aggregating
  # across timesteps). Probably best to do outside the model, downstream following model running, using
  # sample_amount_virus_aggFlight and sample_amount_non_virus_aggFlight as the model outputs
  seq_reads_virus_aggFlight <- if(sample_amount_virus_aggFlight == 0 && sample_amount_non_virus_aggFlight == 0) 0 else seq_tot * (sample_amount_virus_aggFlight * met_bias)/((sample_amount_virus_aggFlight * met_bias) + sample_amount_non_virus_aggFlight)
  seq_reads_non_virus_aggFlight <- seq_tot - seq_reads_virus_aggFlight
  
  ### Stochastic Model Updates for Outputs Relevant to Metagenomic Sequencing
  update(n_inf_AtoB_shedding_events_Out) <- n_inf_AtoB_shedding_events
  update(n_uninf_AtoB_shedding_events_Out) <- n_uninf_AtoB_shedding_events
  update(n_uninf_OtoB_shedding_events_Out) <- n_uninf_OtoB_shedding_events
  update(amount_virus_aggFlight_Out) <- amount_virus_aggFlight
  update(amount_non_virus_aggFlight_Out) <- amount_non_virus_aggFlight
  update(sample_amount_virus_aggFlight_Out) <- sample_amount_virus_aggFlight
  update(sample_amount_non_virus_aggFlight_Out) <- sample_amount_non_virus_aggFlight
  update(seq_reads_virus_aggFlight_Out) <- seq_reads_virus_aggFlight
  update(seq_reads_non_virus_aggFlight_Out) <- seq_reads_non_virus_aggFlight
  
  ### Initial Values for Outputs Relevant to Metagenomic Sequencing
  initial(n_inf_AtoB_shedding_events_Out) <- 0
  initial(n_uninf_AtoB_shedding_events_Out) <- 0
  initial(n_uninf_OtoB_shedding_events_Out) <- 0
  initial(amount_virus_aggFlight_Out) <- 0
  initial(amount_non_virus_aggFlight_Out) <- 0
  initial(sample_amount_virus_aggFlight_Out) <- 0
  initial(sample_amount_non_virus_aggFlight_Out) <- 0
  initial(seq_reads_virus_aggFlight_Out) <- 0
  initial(seq_reads_non_virus_aggFlight_Out) <- 0
  
  ################# Miscellaneous Model Stuff ##################
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt 
  
})


## For Exploration at Later Date

## EXPERIMENTAL - THINKING ABOUT WHETHER ADDED STOCHASTICITY FROM SAMPLING SMALL VOLUMER IS RELEVANT CONSIDERATION
# sample_amount_virus_aggFlight_stoch <- rbinom(amount_virus_aggFlight, samp_frac_aggFlight)
# sample_amount_non_virus_aggFlight_stoch <- rbinom(amount_non_virus_aggFlight, samp_frac_aggFlight)
# seq_reads_virus_aggFlight_stoch <- if(sample_amount_virus_aggFlight_stoch == 0 && sample_amount_non_virus_aggFlight_stoch == 0) 0 else seq_tot * (sample_amount_virus_aggFlight_stoch * met_bias)/((sample_amount_virus_aggFlight_stoch * met_bias) + sample_amount_non_virus_aggFlight_stoch)
# seq_reads_non_virus_aggFlight_stoch <- seq_tot - seq_reads_virus_aggFlight_stoch
