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
    dplyr::mutate(time2 = cut(time, breaks = max(time))) %>% # aggregating by day
    dplyr::mutate(time3 = midpoints(time2)) %>%
    group_by(time2, time3) %>%
    summarise(daily_reads_det = sum(reads_det),
              daily_reads_stoch = sum(reads_stoch),
              daily_flightAB_infections = sum(flightAB_infections),
              daily_flightAB_prevalence = 100 * daily_flightAB_infections/(capacity_per_flight * num_flightsAB),
              daily_community_prevalence = 100 * mean(community_prevalence)/population_size) %>%
    ungroup(time2)
  
  time_to_detection <- daily_df  %>%
    dplyr::select(-time2) %>%
    dplyr::summarise(time_to_detection = time3[min(which(daily_reads_det > num_reads))])
  time_to_detection <- time_to_detection$time_to_detection
  infection_measures <- daily_df %>%
    filter(time3 == time_to_detection)
  time_to_infection <- daily_df  %>%
    dplyr::select(-time2) %>%
    dplyr::summarise(time_to_infection = time3[min(which(daily_flightAB_infections > num_reads))])

  return(list(time = time_to_detection, 
              daily_reads = if(!is.na(time_to_detection)) infection_measures$daily_reads_det else NA,
              daily_flightAB_infections = if(!is.na(time_to_detection)) infection_measures$daily_flightAB_infections else NA,
              daily_flight_AB_prevalence = if(!is.na(time_to_detection)) infection_measures$daily_flightAB_prevalence else NA,
              daily_community_prevalence = if(!is.na(time_to_detection)) infection_measures$daily_community_prevalence else NA,
              time_to_infection = time_to_infection$time_to_infection))
}

lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}


# Function to Run the Model on the Cluster
parallel_model_running <- function(variable_params_list, fixed_params, generator) {
  
  # Loading in set of parameters
  beta <- variable_params_list$beta
  shedding_freq <- variable_params_list$shedding_freq
  ratio_virus_to_non_virus <- variable_params_list$ratio_virus_to_non_virus
  num_flights <- variable_params_list$num_flights
  seq_tot <- variable_params_list$seq_tot
  
  # Figuring out how long to run the model for 
  if (length(fixed_params$end == 1)) {
    temp_end <- fixed_params$end
  } else {
    tol <- 1e-4
    end_index <- ((fixed_params$end$beta - tol) <= beta) & ((fixed_params$end$beta + tol) >= beta)
    temp_end <- fixed_params$end$end[end_index] 
  }
  
  # Generating the model instance
  mod <- generator$new(
    
    # Epidemiological Parameters
    beta = beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
    population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
    
    # Flight Parameters
    capacity_per_flight = fixed_params$capacity_per_flight, 
    num_flights = num_flights, num_flightsAB = num_flights * fixed_params$proportion_AB, 
    samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * num_flights * fixed_params$proportion_AB),
    
    # Sequencing Parameters
    shedding_freq = shedding_freq, virus_shed = fixed_params$non_virus_shed * ratio_virus_to_non_virus, 
    non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
    seq_tot = seq_tot, 
    
    # Miscellaenous Parameters
    dt = fixed_params$dt)
  
  # Setting up temporary storage for stochastic realisation results
  temp_storage <- matrix(data = NA, nrow = fixed_params$stochastic_sim, ncol = 5)
  colnames(temp_storage) <- c("ttd", "num_infs", "flight_prev", "cuminf", "success")
  
  # Running Multiple Stochastic Realisations For One Parameter Set and Summarising the Results 
  for (j in 1:fixed_params$stochastic_sim) {
    
    # Setting Seed
    set.seed(fixed_params$seed[j])
    
    # Running the Model
    output <- mod$run(1:(temp_end/fixed_params$dt))
    output2 <- mod$transform_variables(output)
    
    # Extracting Outputs
    ttd_metrics <- ttd_fun(mod, output2, fixed_params$num_reads)
    
    temp_storage[j, "ttd"] <- ttd_metrics$time
    temp_storage[j, "num_infs"] <- ttd_metrics$daily_flightAB_infections
    temp_storage[j, "flight_prev"] <- ttd_metrics$daily_flight_AB_prevalence
    temp_storage[j, "cuminf"] <- if(!is.na(ttd_metrics$time)) sum(output2$n_SE_Output[1:(ttd_metrics$time/fixed_params$dt)])/fixed_params$population_size else NA
    temp_storage[j, "success"] <- if(!is.na(ttd_metrics$time)) 1 else 0 
    
  }
  
  temp_output <- data.frame(# Params
                            beta = beta,
                            shedding_freq = shedding_freq,
                            ratio_virus_to_non_virus = ratio_virus_to_non_virus,
                            num_flights = num_flights,
                            seq_tot = seq_tot,
    
                            # Time to Detection
                            ttd_lower = if(is.infinite(min(temp_storage[, "ttd"], na.rm = TRUE))) NA else min(temp_storage[, "ttd"], na.rm = TRUE), 
                            ttd_mean = mean(temp_storage[, "ttd"], na.rm = TRUE),
                            ttd_upper = if(is.infinite(max(temp_storage[, "ttd"], na.rm = TRUE))) NA else max(temp_storage[, "ttd"], na.rm = TRUE),
                            
                            # Number of Infections @ Time of Detection
                            lower_inf = if(is.infinite(min(temp_storage[, "num_infs"], na.rm = TRUE))) NA else min(temp_storage[, "num_infs"], na.rm = TRUE),
                            avg_inf = mean(temp_storage[, "num_infs"], na.rm = TRUE),
                            upper_inf = if(is.infinite(max(temp_storage[, "num_infs"], na.rm = TRUE))) NA else max(temp_storage[, "num_infs"], na.rm = TRUE), 
                            
                            # Flight Prevalence @ Time of Detection
                            lower_flight_prev = if(is.infinite(min(temp_storage[, "flight_prev"], na.rm = TRUE))) NA else min(temp_storage[, "flight_prev"], na.rm = TRUE),
                            avg_flight_prev = mean(temp_storage[, "flight_prev"], na.rm = TRUE), 
                            upper_flight_prev = if(is.infinite(max(temp_storage[, "flight_prev"], na.rm = TRUE))) NA else max(temp_storage[, "flight_prev"], na.rm = TRUE),
                            
                            # Cumulative Incidence in Location A @ Time of Detection
                            cuminf_lower =  if(is.infinite(min(temp_storage[, "cuminf"], na.rm = TRUE))) NA else min(temp_storage[, "cuminf"], na.rm = TRUE),
                            cuminf_mean = mean(temp_storage[, "cuminf"], na.rm = TRUE),
                            cuminf_upper = if(is.infinite(max(temp_storage[, "cuminf"], na.rm = TRUE))) NA else max(temp_storage[, "cuminf"], na.rm = TRUE),
                            
                            # Percentage of Stochastic Sims Where Detection Successfully Achieved
                            perc_success = 100 * sum(temp_storage[, "success"]) / fixed_params$stochastic_sim)
  
  # Returning the Output
  return(temp_output)
  
}

wrapped_parallel <- function(variable_params_list, fixed_params, generator, filename, cluster) {
  output <- parLapply(cluster, variable_params_list, parallel_model_running, fixed_params, generator)
  output_list <- list(variable_params = bind_rows(variable_params_list),
                      fixed_params = fixed_params, 
                      model_output = output)
  saveRDS(output_list, file = filename)
  return(output_list)
}


# Old Function to Run the Model on the Cluster
# old_cluster_model_running <- function(fixed_params, variable_params, generator, save_output) {
#   
#   # Setting Up Storage for Model Outputs
#   variable_params$num_reads <- fixed_params$num_reads
#   variable_params$ttd_lower <- NA
#   variable_params$ttd_mean <- NA
#   variable_params$ttd_upper <- NA
#   variable_params$cuminf_lower <- NA
#   variable_params$cuminf_mean <- NA
#   variable_params$cuminf_upper <- NA
#   variable_params$avg_inf <- NA
#   variable_params$lower_inf <- NA
#   variable_params$upper_inf <- NA
#   variable_params$avg_flight_prev <- NA
#   variable_params$lower_flight_prev <- NA
#   variable_params$upper_flight_prev <- NA
#   variable_params$perc_success <- NA
#   
#   # Setting Up the Model
#   for (i in 1:dim(variable_params)[1]) {
#     
#     # Loading in set of parameters
#     beta <- variable_params$beta_sens[i]
#     shedding_freq <- variable_params$shedding_sens[i]
#     ratio_virus_to_non_virus <- variable_params$ratio_sens[i]
#     num_flights <- variable_params$num_flights_sens[i]
#     seq_tot <- variable_params$seq_tot_sens[i]
#     
#     # Figuring out how long to run the model for 
#     tol <- 1e-4
#     end_index <- ((fixed_params$end$beta - tol) <= beta) & ((fixed_params$end$beta + tol) >= beta)
#     temp_end <- fixed_params$end$end[end_index]
#     
#     # Generating the model instance
#     mod <- generator$new(
#       
#       # Epidemiological Parameters
#       beta = beta, gamma = fixed_params$gamma, sigma = fixed_params$sigma, 
#       population_size = fixed_params$population_size, start_infections = fixed_params$start_infections,
#       
#       # Flight Parameters
#       capacity_per_flight = fixed_params$capacity_per_flight, 
#       num_flights = num_flights, num_flightsAB = num_flights * fixed_params$proportion_AB, 
#       samp_frac_aggFlight = fixed_params$sample_flight_ww/(fixed_params$vol_flight_ww * num_flights * fixed_params$proportion_AB),
#       
#       # Sequencing Parameters
#       shedding_freq = shedding_freq, virus_shed = fixed_params$non_virus_shed * ratio_virus_to_non_virus, 
#       non_virus_shed = fixed_params$non_virus_shed, met_bias = fixed_params$met_bias, 
#       seq_tot = seq_tot, 
#       
#       # Miscellaenous Parameters
#       dt = fixed_params$dt)
#     
#     # Setting up temporary storage for stochastic realisation results
#     temp_storage <- matrix(data = NA, nrow = fixed_params$stochastic_sim, ncol = 5)
#     colnames(temp_storage) <- c("ttd", "num_infs", "flight_prev", "cuminf", "success")
#     
#     # Running Multiple Stochastic Realisations For One Parameter Set and Summarising the Results 
#     for (j in 1:fixed_params$stochastic_sim) {
#       
#       # Setting Seed
#       set.seed(fixed_params$seed[j])
#       
#       # Running the Model
#       end <- fixed_params$end$end[which(fixed_params$end$beta == beta)]
#       output <- mod$run(1:(temp_end/fixed_params$dt))
#       output2 <- mod$transform_variables(output)
#       
#       # Extracting Outputs
#       ttd_metrics <- ttd_fun(mod, output2, fixed_params$num_reads)
#       
#       temp_storage[j, "ttd"] <- ttd_metrics$time
#       temp_storage[j, "num_infs"] <- ttd_metrics$daily_flightAB_infections
#       temp_storage[j, "flight_prev"] <- ttd_metrics$daily_flight_AB_prevalence
#       temp_storage[j, "cuminf"] <- if(!is.na(ttd_metrics$time)) sum(output2$n_SE_Output[1:(ttd_metrics$time/fixed_params$dt)])/fixed_params$population_size else NA
#       temp_storage[j, "success"] <- if(!is.na(ttd_metrics$time)) 1 else 0 
#       
#     }
#     
#     variable_params$ttd_lower[i] <- min(temp_storage[, "ttd"], na.rm = TRUE)
#     variable_params$ttd_mean[i] <- mean(temp_storage[, "ttd"], na.rm = TRUE)
#     variable_params$ttd_upper[i] <- max(temp_storage[, "ttd"], na.rm = TRUE)
#     variable_params$lower_inf[i] <- min(temp_storage[, "num_infs"], na.rm = TRUE)
#     variable_params$avg_inf[i] <- mean(temp_storage[, "num_infs"], na.rm = TRUE)
#     variable_params$upper_inf[i] <- max(temp_storage[, "num_infs"], na.rm = TRUE)
#     variable_params$lower_flight_prev[i] <- min(temp_storage[, "flight_prev"], na.rm = TRUE)
#     variable_params$avg_flight_prev[i] <- mean(temp_storage[, "flight_prev"], na.rm = TRUE)
#     variable_params$upper_flight_prev[i] <- max(temp_storage[, "flight_prev"], na.rm = TRUE)
#     variable_params$cuminf_lower[i] <- min(temp_storage[, "cuminf"], na.rm = TRUE)
#     variable_params$cuminf_mean[i] <- mean(temp_storage[, "cuminf"], na.rm = TRUE)
#     variable_params$cuminf_upper[i] <- max(temp_storage[, "cuminf"], na.rm = TRUE)
#     variable_params$perc_success[i] <- 100 * sum(temp_storage[, "success"]) / stochastic_sim
#     
#     rm(temp_storage)
#     print(i) 
#   }
#   
#   # Saving and Writing the Output for the MCMC Chains
#   if (save_output) {
#     saveRDS(variable_params, file = "cluster_test_output.rds")
#   }
#   
#   # Returning the Output
#   return(variable_params)
#   
# }
