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
library(odin.dust); library(mcstate); library(tidyverse); library(tictoc); library(parallel); library(patchwork)
library(hrbrthemes); library(viridis); library(emdbook); library(ecoflux); library(ggExtra); library(patchwork)
source("functions/helper_functions.R")
source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R")
options(dplyr.summarise.inform = FALSE)

# Fixed Parameters 
stochastic_sim <- 150
fixed_params <- list(# Misc Params
                     num_reads = 5,
                     dt = 0.2, 
                     end = 350,
                     seed = rpois(stochastic_sim, 200) * rpois(stochastic_sim, 200) * rpois(stochastic_sim, 20),
                     stochastic_sim = stochastic_sim, 
                     start_infections = 10, 
                     
                     # Epi Params
                     beta = 0.5,    # corresponding to R0 for 2.5 for fixed analyses
                     gamma = 0.25,
                     sigma = 0.2,
                     population_size = 10^7, 
                     
                     # Flight Params
                     num_flights = NA,
                     num_flightsAB = NA,
                     vol_flight_ww = 800, 
                     sample_flight_ww = 1, 
                     capacity_per_flight = 100, 
                     
                     # Shedding Params
                     non_virus_shed = 2*10^11, 
                     shedding_prop = 0.75, 
                     shedding_freq = 1.2,   # based off Janvi's spreadsheet
                     ratio_virus_to_non_virus = 1 * 10^-7,   # based off Janvi's spreadsheet - need to change to 10^-8 I think
                     
                     # Sequencing Params
                     met_bias = 1,
                     seq_tot = 1 * 10^9) # 10^9 reads per day; corresponding to approx $2000 per day, assuming 200bp read length and $10 per Gbp (Illumina-like)

# Infilling Flight Parameters
prop_pop_flying <- 0.005   # proportion of population taking flight every day - fixed, to 0.5% (based off Janvi's spreadsheet) i.e. 1 in every 200 people.
prop_pop_flying_to_AB <- 0.03 # proportion of flight taking population who take flights to location with NAO - fixed to 3% based off Janvi's spreadsheet
overall_prop_pop_AB <- prop_pop_flying * prop_pop_flying_to_AB # equivalent to 1 in every 10,000 of the population taking flight to the location
num_flights <- floor((prop_pop_flying * fixed_params$population_size) / fixed_params$capacity_per_flight)
num_flights_AB <- num_flights * prop_pop_flying_to_AB
fixed_params$num_flights <- num_flights
fixed_params$num_flightsAB <- num_flights_AB

# Sense Check Parameters
R0_sc <- fixed_params$beta/fixed_params$sigma  # R0
flight_sc <- 100 * fixed_params$num_flightsAB * fixed_params$capacity_per_flight / fixed_params$population_size  # % of Pop Travelling A->B every day
virus_shed_sc <- fixed_params$non_virus_shed * fixed_params$ratio_virus_to_non_virus # Amount of virus shed per event


# Create Overall Set of Parameter Values and Run Bivariate Sensitivity Analyses

## Generating parameter sets
seq_sens <- round(lseq(10^8, 10^10, 20))
ratio_sens <- lseq(10^-9, 10^-6, 20)
ratio_seqTotal_params <- expand.grid(ratio_virus_to_non_virus = ratio_sens, seq_tot = seq_sens)
ratio_seqTotal_params$shedding_prop <- fixed_params$shedding_prop
ratio_seqTotal_params$beta <- fixed_params$beta
ratio_seqTotal_params$num_flights <- fixed_params$num_flights
ratio_seqTotal_params_list <- vector(mode = "list", length = dim(ratio_seqTotal_params)[1])
for (i in 1:length(ratio_seqTotal_params_list)) {
  ratio_seqTotal_params_list[[i]]$beta <- ratio_seqTotal_params$beta[i]
  ratio_seqTotal_params_list[[i]]$seq_tot <- ratio_seqTotal_params$seq_tot[i]
  ratio_seqTotal_params_list[[i]]$shedding_prop <- ratio_seqTotal_params$shedding_prop[i]
  ratio_seqTotal_params_list[[i]]$ratio_virus_to_non_virus <- ratio_seqTotal_params$ratio_virus_to_non_virus[i]
  ratio_seqTotal_params_list[[i]]$num_flights <- ratio_seqTotal_params$num_flights[i]
}

## Running model in parallel across all parameter values
numCores <- 6
cluster <- makeCluster(numCores)
clusterEvalQ(cluster, {
  library(Rcpp)
  library(odin)
  library(dplyr)
})
clusterEvalQ(cluster, source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"))
clusterEvalQ(cluster, source("functions/helper_functions.R"))
ratio_seqTotal_output <- wrapped_parallel(variable_params_list = ratio_seqTotal_params_list,
                                          fixed_params = fixed_params, 
                                          generator = stoch_seir_dust,
                                          filename = "ratioseqTotal_bivariateSensitivity.rds",
                                          cluster = cluster)
stopCluster(cluster)

## Plotting the output
bivariate_output <- dplyr::bind_rows(ratio_seqTotal_output$model_output)
scaleFUN <- function(x) sprintf("%.0f", x)

# Give extreme colors:
bivar_heatmap <- ggplot(bivariate_output, aes(x = ratio_virus_to_non_virus , y = seq_tot, fill = 100 * cuminf_mean)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(100 * bivariate_output$cuminf_mean, na.rm = TRUE))) +
  scale_x_continuous(trans = "log10", 
                     #limits = c(min(ratio_sens), max(ratio_sens)),
                     labels = scientific_10x(lseq(min(ratio_sens), max(ratio_sens), 1 + log10(max(ratio_sens)/min(ratio_sens)))), 
                     breaks = lseq(min(ratio_sens), max(ratio_sens), 1 + log10(max(ratio_sens)/min(ratio_sens))), 
                     expand = c(min(ratio_sens), max(ratio_sens))) +
  scale_y_continuous(trans = "log10", 
                     #limits = c(min(seq_sens), max(seq_sens)),
                     labels = scientific_10x(lseq(min(seq_sens), max(seq_sens), 1 + log10(max(seq_sens)/min(seq_sens)))), 
                     breaks = lseq(min(seq_sens), max(seq_sens), 1 + log10(max(seq_sens)/min(seq_sens))), 
                     expand = c(0, 0)) +
  labs(x = "Parameter - Shedding Ratio Virus:Non-Virus", y = "Parameter - Total Sequencing") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = margin(0, 15, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(0.05, 'npc'), 
        legend.key.width = unit(0.15, 'npc'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(fill = guide_colourbar(title = "Cumulative Incidence @ Detection (%)", title.position = "top"))

avg_cumInf_ratio <- bivariate_output %>%
  group_by(ratio_virus_to_non_virus) %>%
  summarise(cuminf_mean2 = mean(c(cuminf_mean, rep(0.60, sum(is.na(cuminf_mean)))), na.rm = TRUE),
            cuminf_lower2 = mean(c(cuminf_lower, rep(0.60, sum(is.na(cuminf_lower)))), na.rm = TRUE),
            cuminf_upper2 = mean(c(cuminf_upper, rep(0.60, sum(is.na(cuminf_upper)))), na.rm = TRUE), 
            cuminf_mean = mean(cuminf_mean, na.rm = TRUE),
            cuminf_lower = mean(cuminf_lower, na.rm = TRUE),
            cuminf_upper = mean(cuminf_upper, na.rm = TRUE))
ratio_mar <- ggplot(avg_cumInf_ratio, aes(x = ratio_virus_to_non_virus, y = 100 * cuminf_mean2)) +
  geom_line() +
  geom_ribbon(aes(x = ratio_virus_to_non_virus, ymin = 100 * cuminf_lower2, ymax = 100 * cuminf_upper2), alpha = 0.5) +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0),
                     trans = "log10") +
  labs(x = "", y = "CumInf(%)") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,-25,0),
        axis.title.y = element_text(angle = 360, hjust = 1, vjust = 0.5)) 

avg_cumInf_seq <- bivariate_output %>%
  group_by(seq_tot) %>%
  summarise(cuminf_mean2 = mean(c(cuminf_mean, rep(0.60, sum(is.na(cuminf_mean)))), na.rm = TRUE),
            cuminf_lower2 = mean(c(cuminf_lower, rep(0.60, sum(is.na(cuminf_lower)))), na.rm = TRUE),
            cuminf_upper2 = mean(c(cuminf_upper, rep(0.60, sum(is.na(cuminf_upper)))), na.rm = TRUE),
            cuminf_mean = mean(cuminf_mean, na.rm = TRUE),
            cuminf_lower = mean(cuminf_lower, na.rm = TRUE),
            cuminf_upper = mean(cuminf_upper, na.rm = TRUE))
seq_mar <- ggplot(avg_cumInf_seq, aes(x = seq_tot, y = 100 * cuminf_mean2)) +
  geom_line() +
  geom_ribbon(aes(x = seq_tot, ymin = 100 * cuminf_lower2, ymax = 100 * cuminf_upper2), alpha = 0.5) +
  scale_x_continuous(trans = "log10", expand = c(0, 0),
                     labels = scientific_10x(c(10^5, 10^6, 10^7, 10^8, 10^9), digits = 0),
                     breaks = c(10^5, 10^6, 10^7, 10^8, 10^9),
                     position = "top") +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip() +
  labs(x = "", y = "CumInf(%)") +
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1, 1, 1, 5))

layout <- "
AA#
BB#
CCD
CCD
"
ratio_seq_bivariate_plot <- ratio_mar + plot_spacer() + bivar_heatmap + seq_mar + 
  plot_layout(design = layout, heights = c(1, -0.25, 5), widths = c(1, 5))
ggsave(filename = "figures/virusRatio_seqTotal_bivariate_sensitivity.pdf",
       plot = ratio_seq_bivariate_plot,
       width = 9,
       height = 9)

## Generating parameter sets
seq_sens <- round(lseq(10^8, 10^10, 20))
prop_sens <- seq(0.05, 1, 0.1)
prop_seqTotal_params <- expand.grid(shedding_prop = prop_sens, seq_tot = seq_sens)
prop_seqTotal_params$ratio_virus_to_non_virus <- fixed_params$ratio_virus_to_non_virus
prop_seqTotal_params$beta <- fixed_params$beta
prop_seqTotal_params$num_flights <- fixed_params$num_flights
prop_seqTotal_params_list <- vector(mode = "list", length = dim(prop_seqTotal_params)[1])
for (i in 1:length(prop_seqTotal_params_list)) {
  prop_seqTotal_params_list[[i]]$beta <- prop_seqTotal_params$beta[i]
  prop_seqTotal_params_list[[i]]$seq_tot <- prop_seqTotal_params$seq_tot[i]
  prop_seqTotal_params_list[[i]]$shedding_prop <- prop_seqTotal_params$shedding_prop[i]
  prop_seqTotal_params_list[[i]]$ratio_virus_to_non_virus <- prop_seqTotal_params$ratio_virus_to_non_virus[i]
  prop_seqTotal_params_list[[i]]$num_flights <- prop_seqTotal_params$num_flights[i]
}

## Running model in parallel across all parameter values
numCores <- 6
cluster <- makeCluster(numCores)
clusterEvalQ(cluster, {
  library(Rcpp)
  library(odin)
  library(dplyr)
})
clusterEvalQ(cluster, source("models/2022-09-28_StochSEIR_AggFlight_Metagenomic_Model.R"))
clusterEvalQ(cluster, source("functions/helper_functions.R"))
prop_seqTotal_output <- wrapped_parallel(variable_params_list = prop_seqTotal_params_list,
                                          fixed_params = fixed_params, 
                                          generator = stoch_seir_dust,
                                          filename = "shedProp_seqTotal_bivariateSensitivity.rds",
                                          cluster = cluster)
stopCluster(cluster)

## Plotting the output
bivariate_output <- dplyr::bind_rows(prop_seqTotal_output$model_output)
scaleFUN <- function(x) sprintf("%.0f", x)

# Give extreme colors:
bivar_heatmap <- ggplot(bivariate_output, aes(x = shedding_prop, y = seq_tot, fill = 100 * cuminf_mean)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(100 * bivariate_output$cuminf_mean, na.rm = TRUE))) +
  scale_x_continuous(#limits = c(min(ratio_sens), max(ratio_sens)),
                     labels = seq(0, 1, 0.2),
                     breaks = seq(0, 1, 0.2),
                     expand = c(min(ratio_sens), max(ratio_sens))) +
  scale_y_continuous(trans = "log10", 
                     #limits = c(min(seq_sens), max(seq_sens)),
                     labels = scientific_10x(lseq(min(seq_sens), max(seq_sens), 1 + log10(max(seq_sens)/min(seq_sens)))), 
                     breaks = lseq(min(seq_sens), max(seq_sens), 1 + log10(max(seq_sens)/min(seq_sens))), 
                     expand = c(0, 0)) +
  labs(x = "Parameter - Proportion of Infectious Shedding", y = "Parameter - Total Sequencing") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(0.05, 'npc'), 
        legend.key.width = unit(0.15, 'npc'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(fill = guide_colourbar(title = "Cumulative Incidence @ Detection (%)", title.position = "top"))

avg_cumInf_prop <- bivariate_output %>%
  group_by(shedding_prop) %>%
  summarise(cuminf_mean2 = mean(c(cuminf_mean, rep(0.60, sum(is.na(cuminf_mean)))), na.rm = TRUE),
            cuminf_lower2 = mean(c(cuminf_lower, rep(0.60, sum(is.na(cuminf_lower)))), na.rm = TRUE),
            cuminf_upper2 = mean(c(cuminf_upper, rep(0.60, sum(is.na(cuminf_upper)))), na.rm = TRUE), 
            cuminf_mean = mean(cuminf_mean, na.rm = TRUE),
            cuminf_lower = mean(cuminf_lower, na.rm = TRUE),
            cuminf_upper = mean(cuminf_upper, na.rm = TRUE))
prop_mar <- ggplot(avg_cumInf_prop, aes(x = shedding_prop, y = 100 * cuminf_mean2)) +
  geom_line() +
  geom_ribbon(aes(x = shedding_prop, ymin = 100 * cuminf_lower2, ymax = 100 * cuminf_upper2), alpha = 0.5) +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0),
                     trans = "log10") +
  labs(x = "", y = "CumInf (%)") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,-25,0),
        axis.title.y = element_text(angle = 360, hjust = 1, vjust = 0.5)) 

avg_cumInf_seq <- bivariate_output %>%
  group_by(seq_tot) %>%
  summarise(cuminf_mean2 = mean(c(cuminf_mean, rep(0.60, sum(is.na(cuminf_mean)))), na.rm = TRUE),
            cuminf_lower2 = mean(c(cuminf_lower, rep(0.60, sum(is.na(cuminf_lower)))), na.rm = TRUE),
            cuminf_upper2 = mean(c(cuminf_upper, rep(0.60, sum(is.na(cuminf_upper)))), na.rm = TRUE),
            cuminf_mean = mean(cuminf_mean, na.rm = TRUE),
            cuminf_lower = mean(cuminf_lower, na.rm = TRUE),
            cuminf_upper = mean(cuminf_upper, na.rm = TRUE))
seq_mar <- ggplot(avg_cumInf_seq, aes(x = seq_tot, y = 100 * cuminf_mean2)) +
  geom_line() +
  geom_ribbon(aes(x = seq_tot, ymin = 100 * cuminf_lower2, ymax = 100 * cuminf_upper2), alpha = 0.5) +
  scale_x_continuous(trans = "log10", expand = c(0, 0),
                     labels = scientific_10x(c(10^5, 10^6, 10^7, 10^8, 10^9), digits = 0),
                     breaks = c(10^5, 10^6, 10^7, 10^8, 10^9),
                     position = "top") +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip() +
  labs(x = "", y = "CumInf(%)") +
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1, 1, 1, 5))

layout <- "
AA#
BB#
CCD
CCD
"
prop_seq_bivariate_plot <- prop_mar + plot_spacer() + bivar_heatmap + seq_mar + 
  plot_layout(design = layout, heights = c(1, -0.25, 5), widths = c(1, 5))

ggsave(filename = "figures/shedProp_seqTotal_bivariate_sensitivity.pdf",
       plot = prop_seq_bivariate_plot,
       width = 9,
       height = 9)


