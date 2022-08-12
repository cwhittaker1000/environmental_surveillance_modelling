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
