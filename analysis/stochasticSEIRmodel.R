# If running for the first time, run this line
#install.packages("drat"); #install.packages("odin.dust"); #install.packages("mcstate")

# Load required libraries
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel)

# Sourcing functions
source("functions/helper_functions.R")

# Changing memory limit to handle large outputs
memory.limit(size = 30000)

# Stochastic SEIR
stoch_seir_dust <- odin.dust::odin_dust({
  
  ################## Model Parameters ##################
  
  ## Epidemiological Parameters
  beta <- user()              # probability of a contact successfully transmitting the disease
  gamma <- user()             # rate of transition from Exposed -> Infectious (incubation period)
  sigma <- user()             # rate of transition from Infectious -> Recovered (rate of recovery)
  population_size <- user()   # overall size of population
  start_infections <- user()  # starting number of infections (in Exposed compartment)
  
  ## NAO/Surveillance Parameters
  p_flight <- user()          # probability of individual taking a flight
  p_flightAB <- user()        # conditional on individual taking a flight, what is the probability of that flight 
                              #   going from location A (epidemic start site) to location B (airport where NAO is operational)
  p_excretion <- user()       # probability of infected individual on flight A->B using the toilet in a way that sheds
                              #   virus into the system
  p_detection <- user()       # probability of detecting the virus the infected person shed into the system
  
  
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
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(N) <- S + E + I + R
  update(I_total) <- I_total + n_SE - n_IR
  update(new_exposed) <- n_SE
  update(new_infectious) <- n_EI
  
  ## NAO/Surveillance Parameters
  
  ### Binomial Draws for Number of People Moving Between Compartments
  average_length_infection <- (1/gamma) + (1/sigma)
  n_inf_take_flight <- rbinom(I_total, dt * p_flight / average_length_infection) # number of infected people (E+I) taking flight
  n_take_flight_AB <- rbinom(n_inf_take_flight, p_flightAB)                      # number of infected people (E+I) taking flight to NAO location
  n_excrete <- rbinom(n_take_flight_AB, p_excretion)                             # number flying to NAO location who shed virus en route
  n_detect <- rbinom(n_excrete, p_detection)                                     # number of individuals who shed on NAO bound flight who are detected
  # Note: need to check whether probabilities here need to be multiplied by dt (I think the answer's no except for p_flight)
  
  ### Stochastic Model Updates for Epidemiological Parameters
  update(num_infs_flight) <- n_inf_take_flight
  update(num_infs_flightAB) <- n_take_flight_AB
  update(num_infs_flightAB_excrete) <- n_excrete
  update(num_infs_flightAB_excrete_detect) <- n_detect
  
  ################# Initial Values ##################
  
  ## Epidemiological Parameters
  initial(N) <- S + E + I + R                       # total population size
  initial(S) <- population_size - start_infections  
  initial(E) <- start_infections                    # starting number of infections
  initial(I) <- 0
  initial(I_total) <- E + I
  initial(R) <- 0
  initial(new_exposed) <- 0
  initial(new_infectious) <- 0
  
  ## NAO/Surveillance Parameters
  initial(num_infs_flight) <- 0
  initial(num_infs_flightAB) <- 0
  initial(num_infs_flightAB_excrete) <- 0
  initial(num_infs_flightAB_excrete_detect) <- 0
  
  ## Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
})

# Note important to distinguish between the number of infected individuals and 
#      number of new infections at each timepoint
# Note that as epidemic spreads, important to take account of 2nd order routes that virus could reach
#      a location. I.e., in early stage of spread, A -> B will be dominant route; later on, when spreading in C,
#      A->C->B also represents a viable route. Think this will lead to systematic underestimation of detection
#      probability later on in epidemic. Implicit assumption for now is that in the earliest stages of spread,
#      dynamics are dominated by transmission explicitly from the initial source location.
# Note that we currently assume detection as a series of independent bernoulli trials for each infected individual.
#      This is line with Mike and Adam's initial document - however, there's probably some non-linear relationship
#      that links number of infected people on flight (note we don't represent flights explicitly currently) and 
#      p_detection (i.e. they're not independent and our focus is not on number detected per se but overall prob of detection).
# Note - does timestep matter for the binomia draws for n_detect etc? As long as I'm summing up to produce day as output
#        then I don't think so, but need to check this! 
# Note - assume if someone gets on the plane when infected, they don't cease to become infected before they've shed on flight. 
# Note - when someone is on the plane going elsewhere, they're not contributing to FOI in location they've come from.
#        Need to figure out whether this is important/relevant to the framework considered here. 
# Note - stochastic fadeout is a thing - I can remove the runs where fadeout occurs in postprocessing, but important to be aware of
# Note - generate canonical parameter sets for each "threat" - e.g. 500 draws from MVN for set of parameters with some covariance
#        Note light tails of MVN and think about whether we need a different distribution instead (with heavier tails and able to consider low prob high impact risks)

# Running and Checking Output of the Model For 1 Parameter Set
var_names <- c("N", "S", "E", "I", "I_total", "R", "new_exposed", "new_infectious", 
               "num_infs_flight", "num_infs_flight_AB", "num_infs_flightAB_excrete", "num_infs_flightAB_excrete_detect",
               "time")
dt <- 0.05
time_period <- 100
time <- time_period/dt
population_size <- 10^5
start_infections <- 5
epi_params <- list(beta = 0.5, 
                   gamma = 0.2, 
                   sigma = 0.4, 
                   population_size = population_size, 
                   start_infections = start_infections, 
                   dt = dt)
nao_params <- list(p_flight = 0.5,
                   p_flightAB = 0.5,
                   p_excretion = 0.9,
                   p_detection = 0.1)
params <- append(epi_params, nao_params)

# Running the Model For One Parameter Set
mod <- stoch_seir_dust$new(pars = params, step = 1, n_particles = 1L, n_threads = 2L, seed = 100L)
mod$update_state(step = 0)
raw_output <- mod$simulate(step_end = 0:time)
output2 <- t(matrix(raw_output, ncol = time+1))
output2 <- data.frame(output2)
colnames(output2) <-  var_names# this should match the ordering in the model above

# Plotting Overall Epidemic Dynamics
plot(output2$time, output2$S, ylim = c(0, population_size), type = "l")
lines(output2$time, output2$I_total, col = "blue")
lines(output2$time, output2$R, col = "purple")

# Plotting Incidence of New Infections/Infectious
time_span <- 50
new_infections_daily <- get_daily_outputs(time_period, dt, "new_exposed", output2)
new_infectious_daily <- get_daily_outputs(time_period, dt, "new_infectious", output2)
plot(output2$time[1:(time_span/dt)], output2$I[1:(time_span/dt)], col = "blue", pch = 20, type = "l") # number of infected people over time
points(seq(0.5, time_span-0.5), new_infections_daily[1:time_span], col = "black", pch = 20) # number of new infections (S -> E transitions)
points(seq(0.5, time_span-0.5), new_infectious_daily[1:time_span], col = "red", pch = 20) # number of newly infectious people (E -> I transitions)

# Plotting Infections Flying to Boston and Detected by NAO
new_infections_flying <- get_daily_outputs(time_period, dt, "num_infs_flight", output2)
new_infections_flyingAB <- get_daily_outputs(time_period, dt, "num_infs_flight_AB", output2)
new_infections_flyingAB_excrete <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete", output2)
new_infections_flyingAB_excrete_detect <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete_detect", output2)

plot(1:time_span, new_infections_flying[1:time_span], col = "blue", pch = 20, type = "l") 
points(1:time_span, new_infections_flyingAB[1:time_span], col = "black", pch = 20, type = "l") 
points(1:time_span, new_infections_flyingAB_excrete[1:time_span], col = "red", pch = 20, type = "l") 
points(1:time_span, new_infections_flyingAB_excrete_detect[1:time_span], col = "purple", pch = 20, type = "l") 

plot(output2$time[1:(time_span/dt)], output2$I_total[1:(time_span/dt)], col = "blue", type = "l", ylim = c(0, max(c(output2$I_total[1:(time_span/dt)], new_infections_flying[1:time_span]))))
lines(1:time_span, new_infections_flying[1:time_span])

# Checking that the proportions for each of the downstrean NAO probabilities match up
sum(new_infections_flying)/sum(new_infections_daily)
sum(new_infections_flyingAB)/sum(new_infections_flying)
sum(new_infections_flyingAB_excrete)/sum(new_infections_flyingAB)
sum(new_infections_flyingAB_excrete_detect)/sum(new_infections_flyingAB_excrete)

# Running the Model For One Parameter Set But Multiple Times
tic()
mod <- stoch_seir_dust$new(pars = params, step = 1, n_particles = 250L, n_threads = 8L, seed = 1L) # need to set seed or will that run same seed for all particles?
mod$update_state(step = 0)
raw_output <- mod$simulate(step_end = 0:time)
toc()

tic()
output <- aperm(raw_output, c(3, 1, 2))
toc()

for (i in 1:250) {
  if (i == 1) {
    plot(output[, 2, i], type = "l")
  } else {
    lines(output[, 2, i], type = "l")
  }
}

# Note need to adjust and remove all the runs that went extinct by themselves
tic()
threshold <- epi_params$population_size - epi_params$start_infections - (epi_params$population_size * 0.005)
runs_retained <- apply(output, 3, function(x) {
  min(x[, 2]) < threshold
})
# sum(runs_retained)
retained_runs <- which(runs_retained)
toc()

for (i in retained_runs) {
  if (i == min(retained_runs)) {
    plot(output[, 2, i], type = "l")
  } else {
    lines(output[, 2, i], type = "l")
  }
}

# Generating Summary Outputs
tic()
summary_outputs <- apply(output[, , retained_runs], c(1, 2), quantile, c(0.5, 0.05, 0.95))
toc()

tic()
cl <- makeCluster(8)
summary_outputs2 <- parApply(cl, output[, , retained_runs], c(1, 2), quantile, c(0.5, 0.05, 0.95))
stopCluster(cl)
toc()

median <- summary_outputs[1, , ]
median <- data.frame(median)
lower <- summary_outputs[2, , ]
lower <- data.frame(lower)
upper <- summary_outputs[3, , ]
upper <- data.frame(upper)
colnames(median) <- colnames(lower) <- colnames(upper) <- var_names

median_new_infections_flying <- get_daily_outputs(time_period, dt, "num_infs_flight", median)
median_new_infections_flyingAB <- get_daily_outputs(time_period, dt, "num_infs_flight_AB", median)
median_new_infections_flyingAB_excrete <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete", median)
median_new_infections_flyingAB_excrete_detect <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete_detect", median)

lower_new_infections_flying <- get_daily_outputs(time_period, dt, "num_infs_flight", lower)
lower_new_infections_flyingAB <- get_daily_outputs(time_period, dt, "num_infs_flight_AB", lower)
lower_new_infections_flyingAB_excrete <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete", lower)
lower_new_infections_flyingAB_excrete_detect <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete_detect", lower)

upper_new_infections_flying <- get_daily_outputs(time_period, dt, "num_infs_flight", upper)
upper_new_infections_flyingAB <- get_daily_outputs(time_period, dt, "num_infs_flight_AB", upper)
upper_new_infections_flyingAB_excrete <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete", upper)
upper_new_infections_flyingAB_excrete_detect <- get_daily_outputs(time_period, dt, "num_infs_flightAB_excrete_detect", upper)

par(mfrow = c(1, 4), mar = c(2, 2, 2, 2))
plot(median_new_infections_flying, type = "l", ylim = c(0, max(upper_new_infections_flying)))
lines(lower_new_infections_flying, type = "l")
lines(upper_new_infections_flying, type = "l")

plot(median_new_infections_flyingAB, type = "l", ylim = c(0, max(upper_new_infections_flying)))
lines(lower_new_infections_flyingAB, type = "l")
lines(upper_new_infections_flyingAB, type = "l")

plot(median_new_infections_flyingAB_excrete, type = "l", ylim = c(0, max(upper_new_infections_flying)))
lines(lower_new_infections_flyingAB_excrete, type = "l")
lines(upper_new_infections_flyingAB_excrete, type = "l")

plot(median_new_infections_flyingAB_excrete_detect, type = "l", ylim = c(0, max(upper_new_infections_flying)))
lines(lower_new_infections_flyingAB_excrete_detect, type = "l")
lines(upper_new_infections_flyingAB_excrete_detect, type = "l")
toc()

# Running the Model With Multiple Parameter Sets

# Creating the parameter sets
tic()
var_names <- c("N", "S", "E", "I", "I_total", "R", 
               "new_exposed", "new_infectious", 
               "num_infs_flight", "num_infs_flight_AB", "num_infs_flightAB_excrete", "num_infs_flightAB_excrete_detect",
               "time")
dt <- 0.05
time_period <- 200
time <- time_period/dt
population_size <- 10^5
start_infections <- 5
nao_params <- list(p_flight = 0.5, p_flightAB = 0.5, p_excretion = 0.9, p_detection = 0.1)
beta_range <- seq(0.5, 1.5, 0.25)
gamma_range <- seq(0.2, 0.8, 0.1)
sigma_range <- seq(0.2, 0.8, 0.1)
population_size <- 10^5
start_infections <- 5
param_combos <- expand_grid(beta = beta_range,
                            gamma = gamma_range,
                            sigma = sigma_range,
                            population_size = population_size,
                            start_infections = start_infections,
                            dt = dt)
params_list <- vector(mode = "list", length = nrow(param_combos))
for (i in 1:nrow(param_combos)) {
  temp_list <- as.list(param_combos[i, ])
  params_list[[i]] <- append(temp_list, nao_params)
}

# Running the model
particles <- 250
mod <- stoch_seir_dust$new(pars = params_list, step = 1L, 
                           n_particles = particles, 
                           pars_multi = TRUE, 
                           n_threads = 6L,
                           seed = 10L) # need to set seed or will that run same seed for all particles?
mod$update_state(step = 0)
mod$set_index(c(S = 2, I_total = 5, new_inf = 7))
raw_output <- mod$simulate(step_end = 0:time) # state variables * particles per param set * param combos * timepoints
dim(raw_output)

# Removing stochastic fadeout
threshold <- epi_params$population_size - epi_params$start_infections - (epi_params$population_size * 0.005)
dims <- dim(raw_output)
runs_retained <- (raw_output["S", , , dims[4]]) < threshold
prop_runs_retained <- apply(runs_retained, 2, sum)/num_particles
runs_retained_vector <- as.vector(runs_retained)
runs_retained_vector2 <- rep(runs_retained_vector, each = num_state_vars)

# Creating overall dataset
x <- array_flatten(raw_output, c(1, 2, 3))
num_particles <- particles
num_param_combos<- dim(param_combos)[1]
state_var_names <- names(mod$index())
num_state_vars <- length(state_var_names)
rm(raw_output)
y <- data.frame(state_variable = rep(state_var_names, num_particles * num_param_combos),
                particle_number = rep(rep(1:num_particles, each = num_state_vars), num_param_combos),
                param_combo = rep(1:num_param_combos, each = num_state_vars * num_particles),
                retain_run = runs_retained_vector2,
                prop_runs_retained = rep(prop_runs_retained, each = num_state_vars * num_particles),
                num_runs_retained = rep(apply(runs_retained, 2, sum), each = num_state_vars * num_particles),
                beta = rep(param_combos$beta, each = num_state_vars * num_particles),
                gamma = rep(param_combos$gamma, each = num_state_vars * num_particles),
                sigma = rep(param_combos$sigma, each = num_state_vars * num_particles),
                p_flight = nao_params$p_flight,
                p_flightAB = nao_params$p_flightAB,
                p_excretion = nao_params$p_detection,
                p_detection = nao_params$p_excretion,
                x)
rm(x)
n <- which(runs_retained_vector2)
z <- y[n, ]
rm(y)
toc()

gc()

format(object.size(z), standard = "legacy", units = "Mb")

tail(colnames(z))
new <- z %>% 
  pivot_longer(cols = -c(state_variable:p_detection), names_to = "timestep", values_to = "val") %>%
  mutate(timestep = as.numeric(gsub("X", "", timestep)),
         time = (timestep - 1) * dt)

data.table::setDT(z)
ab <- data.table::melt(z, id.vars = c("state_variable", "particle_number", "param_combo",
                          "retain_run", "prop_runs_retained", "num_runs_retained", 
                          "beta", "gamma", "sigma", "p_flight", "p_flightAB", 
                          "p_excretion", "p_detection"))

rm(z)
gc()


# creating overall dataset
# output <- aperm(raw_output, c(4, 1, 2, 3)) # timepoints * variables * particles * param sets
# dim(output)
# runs_retained <- apply(output, c(3, 4), function(x) {
#   min(x[, 2]) < threshold
# })
# dim(runs_retained)
# sum(runs_retained)
# apply(runs_retained, 2, sum)

# temp <- raw_output[, , 1, ]
# dim(temp)
# for (i in 1:100) {
#   if (i == 1) {
#     plot(temp[2, i, ], type = "l", ylim = c(0, 10^5))
#   } else {
#     lines(temp[2, i, ], type = "l")
#   } 
# }
# param_combos[2, ]
# apply(runs_retained, 2, sum) # something not quite working here - index position 2 = 0, 
#                              # yet param_combos has R0 > 1 and manual plotting above suggests there are non-fadeout runs
#                              # *Ah* - I think they're growing they just haven't hit the threshold yet because R0 is marginal 
#                              # (approx 1.2) and so with low starting infections, takes ages to increase
# # NEED TO FIGURE OUT WHAT TO DO RE THIS
# 
# 
# 
# # 
# # runs_retained <- apply(raw_output, c(2, 3), function(x) {
# #   min(x[2, ]) < threshold
# # })
# # dim(runs_retained)
# # sum(runs_retained)
# # apply(runs_retained, 2, sum)
# retained_index <- which(runs_retained)
# toc()
# 
# dim(output)
# 
# dim(output)
# # timepoints x state vars x particles x param sets
# 
# dim(runs_retained)
# # particles x param sets
# 
# dim(output[, , , 1])
# length(runs_retained[, 1])
# 
# dim(output[, , runs_retained[, 1], 1]) # correctly retained particles for param set 1
# 
# temp <- output[, , runs_retained]
# 
# tic()
# cl <- makeCluster(8)
# 
# summary_outputs2 <- parApply(cl, output[, , retained_runs], c(1, 2), quantile, 0.5)
# stopCluster(cl)
# toc()
# 
# 
# # two types of summary - per param set and overall
# 
# tic()
# summary_outputs <- apply(output[, , retained_runs], c(1, 2), quantile, c(0.5, 0.05, 0.95))
# toc()
# 
# tic()
# cl <- makeCluster(8)
# summary_outputs2 <- parApply(cl, output[, , retained_runs], c(1, 2), quantile, c(0.5, 0.05, 0.95))
# stopCluster(cl)
# toc()
# 
# 
# dim(output)
# 
# mod$update_state(step = 0)
# mod$step()
# 
# tic()
# output <- mod$simulate(step_end = 1:time)
# toc()
# dim(output)
# 
# output2 <- aperm(output, c(4, 1, 3, 2))
# output2 <- output2[, , , 1]
# dim(output2)
# colnames(output2) <- c("N", "S", "E", "I", "I_total", "I_total2", "R", "time")
# head(output2)
# 
# plot(output2[, "time", 1], output2[, "S", 1], ylim = c(0, population_size), type = "l")
# lines(output2[, "time", 1], output2[, "I_total2", 1], col = "black")
# lines(output2[, "time", 1], output2[, "R", 1], col = "black")
# 
# lines(output2[, "time", 2], output2[, "S", 2], ylim = c(0, population_size), type = "l", col = "red")
# lines(output2[, "time", 2], output2[, "I_total2", 2], col = "red")
# lines(output2[, "time", 2], output2[, "R", 2], col = "red")
# 
# dim(output2)
# 







