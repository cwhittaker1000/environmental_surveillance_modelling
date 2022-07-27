# If running for the first time, run this line
#install.packages("drat"); #install.packages("odin.dust"); #install.packages("mcstate")

# Load required libraries
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc)

# Sourcing functions
source("functions/helper_functions.R")

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
  n_inf_take_flight <- rbinom(I_total, p_flight)             # number of infected people (E+I) taking flight
  n_take_flight_AB <- rbinom(n_inf_take_flight, p_flightAB)  # number of infected people (E+I) taking flight to NAO location
  n_excrete <- rbinom(n_take_flight_AB, p_excretion)         # number flying to NAO location who shed virus en route
  n_detect <- rbinom(n_excrete, p_detection)                 # number of individuals who shed on NAO bound flight who are detected
  # Note: need to check whether probabilities here need to be multiplied by dt (I think the answer's no)
  
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

# Running and Checking Output of the Model For 1 Parameter Set
dt <- 0.02
time_period <- 200
time <- time_period/dt
population_size <- 10^5
start_infections <- 2
epi_params <- list(beta = 1, 
                   gamma = 0.5, 
                   sigma = 0.5, 
                   population_size = population_size, 
                   start_infections = start_infections, 
                   dt = dt)
average_length_infection <- (1/epi_params$gamma) + (1/epi_params$sigma)
nao_params <- list(p_flight = dt * 0.5 / average_length_infection, 
                   # need to be a little careful of this one - needs adjusting for average length of infection
                   # also needs adjusting for something else
                   p_flightAB = 0.1,
                   p_excretion = 1,
                   p_detection = 0.1)
params <- append(epi_params, nao_params)

# Running the Model For One Parameter Set
mod <- stoch_seir_dust$new(pars = params, step = 1, n_particles = 1L, n_threads = 1L, seed = 10L)
mod$update_state(step = 0)
raw_output <- mod$simulate(step_end = 0:time)
output2 <- t(matrix(raw_output, ncol = time+1))
output2 <- data.frame(output2)
colnames(output2) <- c("N", "S", "E", "I", "I_total", "R", "new_exposed", "new_infectious", 
                       "num_infs_flight", "num_infs_flight_AB", "num_infs_flightAB_excrete", "num_infs_flightAB_excrete_detect",
                       "time") # this should match the ordering in the model above

# Plotting Overall Epidemic Dynamics
plot(output2$time, output2$S, ylim = c(0, population_size), type = "l")
lines(output2$time, output2$I_total, col = "blue")
lines(output2$time, output2$R, col = "purple")

# Plotting Incidence of New Infections/Infectious
new_infections_daily <- get_daily_outputs(time_period, dt, "new_exposed", output2)
new_infectious_daily <- get_daily_outputs(time_period, dt, "new_infectious", output2)
time_span <- 25
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

plot(output2$time[1:(time_span/dt)], output2$I_total[1:(time_span/dt)], col = "blue", type = "l",
     ylim = c(0, max(c(output2$I_total[1:(time_span/dt)],
                       new_infections_flying[1:time_span]))))
lines(1:time_span, new_infections_flying[1:time_span])

y <- (output2$num_infs_flight/output2$I_total)
mean(y, na.rm = TRUE) # currently p_flight is specifying the proportion of the infected population
                      # travelling at each timestep
                      # that's not what we want though - I think we want "probability infected individual
                      # makes trip at some point whilst infected
                      # need to multiply by timestep
# gives non-sensical results when e.g new_infections_flying is 50x greater than daily number of infections
sum(new_infections_flying)/sum(new_infections_daily)

sum(new_infections_flyingAB)/sum(new_infections_flying)


sum(new_infections_flyingAB_excrete)/sum(new_infections_flyingAB)

sum(new_infections_flyingAB_excrete_detect)/sum(new_infections_flyingAB_excrete)















# Running the Model With Multiple Parameter Sets
beta_range <- seq(3, 6, 0.5)
gamma_range <- seq(1, 3, 0.1)
sigma_range <- seq(1, 3, 0.1)
population_size <- 10^5
start_infections <- 100
dt <- 0.05
param_combos <- expand_grid(beta = beta_range,
                            gamma = gamma_range,
                            sigma = sigma_range,
                            population_size = population_size,
                            start_infections = start_infections,
                            dt = dt)
params_list <- vector(mode = "list", length = nrow(param_combos))
for (i in 1:nrow(param_combos)) {
  params_list[[i]] <- param_combos[i, ]
}

mod <- stoch_seir_dust$new(pars = params_list, step = 1, n_particles = 1L, n_threads = 4L, seed = 10L, pars_multi = TRUE)
mod$update_state(step = 0)
mod$step()

tic()
output <- mod$simulate(step_end = 1:time)
toc()
dim(output)

output2 <- aperm(output, c(4, 1, 3, 2))
output2 <- output2[, , , 1]
dim(output2)
colnames(output2) <- c("N", "S", "E", "I", "I_total", "I_total2", "R", "time")
head(output2)

plot(output2[, "time", 1], output2[, "S", 1], ylim = c(0, population_size), type = "l")
lines(output2[, "time", 1], output2[, "I_total2", 1], col = "black")
lines(output2[, "time", 1], output2[, "R", 1], col = "black")

lines(output2[, "time", 2], output2[, "S", 2], ylim = c(0, population_size), type = "l", col = "red")
lines(output2[, "time", 2], output2[, "I_total2", 2], col = "red")
lines(output2[, "time", 2], output2[, "R", 2], col = "red")

dim(output2)








