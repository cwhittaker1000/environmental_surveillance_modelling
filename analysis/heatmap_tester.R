# Load required libraries
library(tidyverse); library(hrbrthemes); library(viridis); library(emdbook); library(ecoflux); library(ggExtra)
library(patchwork)

# Load in real data
table(clust_df$shedding_sens)
table(clust_df$num_flights_sens)
table(clust_df$ratio_sens)

clust_df <- readRDS("outputs/full_sensitivity_analysis_cluster.rds") %>%
  filter(round(shedding_sens, 1) == round(1.2,1),
         round(num_flights_sens, 0) == round(150, 0),
         ratio_sens > 2.31012970008316e-06 & ratio_sens < 8.11130830789687e-06)
table(clust_df$beta_sens)
table(clust_df$seq_tot_sens)
clust_df$R0 <- clust_df$beta_sens / 0.2
clust_df$seq_total <- clust_df$seq_tot_sens
clust_df$time_to_detection <- clust_df$ttd_mean

scaleFUN <- function(x) sprintf("%.0f", x)


# Give extreme colors:
x <- ggplot(clust_df, aes(x = R0, y = seq_total, fill = time_to_detection)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(clust_df$time_to_detection, na.rm = TRUE))) +
  # scale_fill_distiller(palette = "RdPu", limits = c(0, max(df$time_to_detection))) +
  # scale_fill_gradient(low = "white", high = "blue", limits = c(0, max(df$time_to_detection))) +
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

avg_ttd_by_R0 <- clust_df %>%
  group_by(R0) %>%
  summarise(time_to_detection = mean(c(time_to_detection, 
                                       rep(200, sum(is.na(time_to_detection)))),
                                     na.rm = TRUE),
            lower = mean(c(ttd_lower,
                           rep(200, sum(is.na(ttd_lower)))), na.rm = TRUE),
            upper = mean(c(ttd_upper,
                           rep(200, sum(is.na(ttd_upper)))), na.rm = TRUE))
y <- ggplot(avg_ttd_by_R0, aes(x = R0, y = time_to_detection)) +
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

yalt <- ggplot(avg_ttd_by_R0, aes(x = R0, y = time_to_detection, fill = time_to_detection)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(clust_df$time_to_detection, na.rm = TRUE))) +
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


avg_ttd_by_seq <- clust_df %>%
  mutate(max_time_to_detection = max(time_to_detection, na.rm = TRUE)) %>%
  group_by(seq_total) %>%
  summarise(time_to_detection = mean(c(time_to_detection, 
                                       rep(200, sum(is.na(time_to_detection)))),
                                     na.rm = TRUE),
            lower = mean(c(ttd_lower,
                           rep(200, sum(is.na(ttd_lower)))), na.rm = TRUE),
            upper = mean(c(ttd_upper,
                           rep(200, sum(is.na(ttd_upper)))), na.rm = TRUE))
            
                                

z <- ggplot(avg_ttd_by_seq, aes(x = seq_total, y = time_to_detection)) +
  geom_line() +
  geom_ribbon(aes(x = seq_total, ymin = lower, ymax = upper), alpha = 0.5) +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(clust_df$time_to_detection, na.rm = TRUE))) +
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

zalt <- ggplot(avg_ttd_by_seq, aes(x = seq_total, y = time_to_detection, fill = time_to_detection)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = seq_total, ymin = lower, ymax = upper)) +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(clust_df$time_to_detection, na.rm = TRUE))) +
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

