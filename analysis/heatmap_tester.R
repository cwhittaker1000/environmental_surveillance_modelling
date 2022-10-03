# Load required libraries
library(tidyverse); library(hrbrthemes); library(viridis); library(emdbook); library(ecoflux); library(ggExtra)
library(patchwork)

# Load in real data
clust_df <- readRDS("outputs/R0_seqTotal_bivariate_sensitivity_analysis.rds") 
bivariate_output <- dplyr::bind_rows(clust_df$model_output)
bivariate_output$R0 <- bivariate_output$beta/clust_df$fixed_params$sigma
scaleFUN <- function(x) sprintf("%.0f", x)

# Give extreme colors:
x <- ggplot(bivariate_output, aes(x = R0 , y = seq_tot, fill = ttd_mean)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option = "magma", 
                     limits = c(0, max(bivariate_output$ttd_mean, na.rm = TRUE))) +
                     #na.value="white") +
  # scale_fill_distiller(palette = "RdPu", limits = c(0, max(df$ttd_mean))) +
  # scale_fill_gradient(low = "white", high = "blue", limits = c(0, max(df$ttd_mean))) +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0)) +
  scale_y_continuous(trans = "log10", expand = c(0, 0),
                     labels = scientific_10x(c(10^5, 10^6, 10^7, 10^8, 10^9), digits = 0),
                     breaks = c(10^5, 10^6, 10^7, 10^8, 10^9)) +
  labs(x = "Parameter - R0", y = "Parameter - Total Sequencing") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(0.05, 'npc'), 
        legend.key.width = unit(0.15, 'npc'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  guides(fill = guide_colourbar(title = "Time to Detection (Days)", title.position = "top"))

avg_ttd_by_R0 <- bivariate_output %>%
  group_by(R0) %>%
  summarise(ttd_mean = mean(c(ttd_mean, rep(250, sum(is.na(ttd_mean)))), na.rm = TRUE),
            lower = mean(c(ttd_lower,
                           rep(250, sum(is.na(ttd_lower)))), na.rm = TRUE),
            upper = mean(c(ttd_upper,
                           rep(250, sum(is.na(ttd_upper)))), na.rm = TRUE),
            inf_mean = mean(c(avg_inf, rep(250, sum(is.na(ttd_mean)))), na.rm = TRUE))
y <- ggplot(avg_ttd_by_R0, aes(x = R0, y = ttd_mean)) +
  geom_line() +
  geom_ribbon(aes(x = R0, ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scaleFUN) +
  labs(x = "", y = "TTD") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,-25,0),
        axis.title.y = element_text(angle = 360, hjust = 1, vjust = 0.5)) 

yalt <- ggplot(avg_ttd_by_R0, aes(x = R0, y = ttd_mean, fill = ttd_mean)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = R0, ymin = lower, ymax = upper)) +
  scale_fill_viridis(discrete=FALSE, option = "magma", 
                     limits = c(0, max(bivariate_output$ttd_mean, na.rm = TRUE))) +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scaleFUN) +
  labs(x = "", y = "TTD") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,-25,0),
        axis.title.y = element_text(angle = 360, hjust = 1, vjust = 0.5),
        legend.position = "none") 


avg_ttd_by_seq <- bivariate_output %>%
  group_by(seq_tot) %>%
  summarise(ttd_mean = mean(c(ttd_mean, 
                              rep(250, sum(is.na(ttd_mean)))), na.rm = TRUE),
            lower = mean(c(ttd_lower,
                           rep(250, sum(is.na(ttd_lower)))), na.rm = TRUE),
            upper = mean(c(ttd_upper,
                           rep(250, sum(is.na(ttd_upper)))), na.rm = TRUE))
            
                                

z <- ggplot(avg_ttd_by_seq, aes(x = seq_tot, y = ttd_mean)) +
  geom_line() +
  geom_ribbon(aes(x = seq_tot, ymin = lower, ymax = upper), alpha = 0.5) +
  scale_fill_viridis(discrete=FALSE, option = "magma", 
                     limits = c(0, max(bivariate_output$ttd_mean, na.rm = TRUE))) +
  scale_x_continuous(trans = "log10", expand = c(0, 0),
                     labels = scientific_10x(c(10^5, 10^6, 10^7, 10^8, 10^9), digits = 0),
                     breaks = c(10^5, 10^6, 10^7, 10^8, 10^9),
                     position = "top") +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip() +
  labs(x = "", y = "TTD") +
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1, 1, 1, 5))

zalt <- ggplot(avg_ttd_by_seq, aes(x = seq_tot, y = ttd_mean, fill = ttd_mean)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = seq_tot, ymin = lower, ymax = upper)) +
  scale_fill_viridis(discrete=FALSE, option = "magma", 
                     limits = c(0, max(bivariate_output$ttd_mean, na.rm = TRUE))) +
  scale_x_continuous(trans = "log10", expand = c(0, 0),
                     labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)) +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip() +
  labs(x = "", y = "TTD") +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1, 1, 1, 5))

layout <- "
AA#
BB#
CCD
CCD
"
y + plot_spacer() + x + z + 
  plot_layout(design = layout, heights = c(1, -0.25, 5), widths = c(1, 5))


yalt + plot_spacer() + x + zalt + 
  plot_layout(design = layout, heights = c(1, -0.25, 5), widths = c(1, 5))

