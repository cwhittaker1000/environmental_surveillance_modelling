# Helper Functions


# Function to extract outputs from incidence outputs from a single model run 
# and convert them from "per timestep" to "daily"
# Takes the following arguments:
#     time_period  - the overall time the model was run for
#     dt           - the step-size used in the model
#     variable     - which variable you want to generate daily outputs for
#     model_output - the model output, processed into a named dataframe.
#                    Note must be a 2D array with nrow = number of timesteps
#                    abd each column a different model output.
get_daily_outputs <- function(time_period, dt, variable, model_output) {
  
  if (!(variable %in% colnames(model_output))) {
    stop("Variable must be in colnames(model_output)")
  }
  if (!(is.data.frame(model_output))) {
    stop("Model output must be a data.frame")
  }
  
  daily_var <- sapply(1:time_period, function(i) {
    if (i == 1) {
      pos <- seq(i, (i/dt))
    } else if (i == time_period) {
      pos <- seq((i-1)/dt + 1, i/dt + 1)
    }
    else {
      pos <- seq((i-1)/dt + 1, i/dt)
    }
    sum(model_output[pos, variable])
  })
  return(daily_var)
  
}

# get midpoints as numerical value after using "cut"
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}

time_to_detection_fun <- function(mod_output, num, inf_or_reads, agg_or_indiv) {
  
  # Checking inputs
  if (!(inf_or_reads %in% c("infections", "reads"))) {
    stop("agg_or_indiv must be either infections or reads")
  } 
  if (!(agg_or_indiv %in% c("individual", "aggregated"))) {
    stop("agg_or_indiv must be either individual or aggregated")
  } 
  
  num_flights <- dim(mod_output$n_inf_specific_flightABOut)[2]
  
  # For aggregated flight outputs
  if (agg_or_indiv == "aggregated") {
    if (inf_or_reads == "infections") {
      agg_flight_infs <- data.frame(time = mod_output$time, mod_output$n_inf_flightABOut)
      colnames(agg_flight_infs) <- c("time", "infections")
      agg_inf_df <- agg_flight_infs %>%
        mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
        group_by(time2) %>%
        summarise(daily_infections = sum(infections)/num_flights) %>%
        mutate(time3 = midpoints(time2)) %>%
        summarise(time_num = time3[min(which(daily_infections > num))])
      time <- unname(agg_inf_df$time_num)
    } else if (inf_or_reads == "reads") {
      agg_flight_reads <- data.frame(time = mod_output$time, mod_output$seq_reads_virus_aggFlight_Out)
      colnames(agg_flight_reads) <- c("time", "reads")
      agg_reads_df <- agg_flight_reads %>%
        mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
        group_by(time2) %>%
        summarise(daily_reads = sum(reads)) %>%
        mutate(time3 = midpoints(time2)) %>%
        summarise(time_num = time3[min(which(daily_reads > num))])
      time <- unname(agg_reads_df$time_num)
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

ttd_fun <- function(mod, mod_output, num_reads) {
  
  num_flightsAB <- mod$contents()$num_flightsAB
  capacity_per_flight <- mod$contents()$capacity_per_flight
  population_size <- mod$contents()$population_size
  
  df <- data.frame(time = mod_output$time, 
                   reads_det = mod_output$seq_reads_virus_aggFlight_det_Out,
                   reads_stoch = mod_output$seq_reads_virus_aggFlight_stoch_Out,
                   flightAB_infections = mod_output$n_inf_flightABOut,
                   community_prevalence = mod_output$I)
  daily_df <- df %>%
    mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
    mutate(time3 = midpoints(time2)) %>%
    group_by(time2, time3) %>%
    summarise(daily_reads_det = sum(reads_det),
              daily_reads_stoch = sum(reads_stoch),
              daily_flightAB_infections = sum(flightAB_infections),
              daily_flightAB_prevalence = 100 * daily_flightAB_infections/(capacity_per_flight * num_flightsAB),
              daily_community_prevalence = 100 * mean(community_prevalence)/population_size) %>%
    ungroup(time2)
  
  time_to_detection <- daily_df  %>%
    select(-time2) %>%
    summarise(time_to_detection = time3[min(which(daily_reads_det > num_reads))])
  time_to_detection <- time_to_detection$time_to_detection
  infection_measures <- daily_df %>%
    filter(time3 == time_to_detection)

  return(list(time = time_to_detection, 
              daily_reads = if(!is.na(time_to_detection)) infection_measures$daily_reads_det else NA,
              daily_flightAB_infections = if(!is.na(time_to_detection)) infection_measures$daily_flightAB_infections else NA,
              daily_flight_AB_prevalence = if(!is.na(time_to_detection)) infection_measures$daily_flightAB_prevalence else NA,
              daily_community_prevalence = if(!is.na(time_to_detection)) infection_measures$daily_community_prevalence else NA))
}

lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}


## Function to be Run on the Cluster
cluster_model_running <- function(fixed_params, variable_params, generator) {
  
  # Setting Up Storage for Model Outputs
  variable_params$num_reads <- fixed_params$num_reads
  variable_params$ttd_lower <- NA
  variable_params$ttd_mean <- NA
  variable_params$ttd_upper <- NA
  variable_params$cuminf_lower <- NA
  variable_params$cuminf_mean <- NA
  variable_params$cuminf_upper <- NA
  variable_params$perc_success <- NA
  
  # Setting Up the Model
  for (i in 1:dim(variable_params)[1]) {
    
    # Loading in set of parameters
    beta <- variable_params$beta[i]
    shedding_freq <- variable_params$shedding_freq[i]
    ratio_virus_to_non_virus <- variable_params$ratio_virus_to_non_virus[i]
    num_flights <- variable_params$num_flights[i]
    num_flightsAB <- variable_params$num_flightsAB[i]
    seq_tot <- variable_params$seq_tot[i]
    
    # Generating the model instance
    mod <- generator$new(
      
      # Epidemiological Parameters
      beta = beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
      population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
      
      # Flight Parameters
      capacity_per_flight = fixed_params$capacity_per_flight, 
      num_flights = num_flights, num_flightsAB = num_flightsAB, 
      samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * num_flightsAB),
      
      # Sequencing Parameters
      shedding_freq = shedding_freq, virus_shed = fixed_params$non_virus_shed * ratio_virus_to_non_virus, 
      non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
      seq_tot = seq_tot, 
      
      # Miscellaenous Parameters
      dt = dt)
    
    # Setting up temporary storage for stochastic realisation results
    temp_storage <- matrix(data = NA, nrow = stochastic_sim, ncol = 3)
    colnames(temp_storage) <- c("ttd", "cuminf", "success")
    
    for (j in 1:fixed_params$stochastic_sim) {
      
      # Setting Seed
      set.seed(fixed_params$seed[j])
      
      # Running the Model
      end <- fixed_params$end$end[which(fixed_params$end$beta == beta)]
      output <- mod$run(1:(end/fixed_params$dt))
      output2 <- mod$transform_variables(output)
      
      # Extracting Outputs
      temp_storage[j, "ttd"] <- time_to_detection_fun(output2, fixed_params$num_reads, "reads", "aggregated")
      temp_storage[j, "cuminf"] <- if(!is.na(temp_storage[j, "ttd"])) sum(output2$n_SE_Output[1:(temp_storage[j, "ttd"]/dt)])/fixed_params$population_size else NA
      temp_storage[j, "success"] <- if(!is.na(temp_storage[j, "cuminf"])) 1 else 0
    }
    
    variable_params$ttd_lower[i] <- min(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$ttd_mean[i] <- mean(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$ttd_upper[i] <- max(temp_storage[, "ttd"], na.rm = TRUE)
    variable_params$cuminf_lower[i] <- min(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$cuminf_mean[i] <- mean(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$cuminf_upper[i] <- max(temp_storage[, "cuminf"], na.rm = TRUE)
    variable_params$perc_success[i] <- 100 * sum(temp_storage[, "success"]) / stochastic_sim
    
    rm(temp_storage)
    print(i) 
  }
  
  # Saving and Writing the Output for the MCMC Chains
  saveRDS(variable_params, file = "cluster_test_output.rds")
  
  # Returning the Output
  return(variable_params)
  
}
