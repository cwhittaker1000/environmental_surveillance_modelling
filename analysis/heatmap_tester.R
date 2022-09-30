# Load required libraries
library(tidyverse); library(hrbrthemes); library(viridis); library(emdbook); library(ecoflux); library(ggExtra)
library(patchwork)

# Generate dummy data 
R0 <- seq(0.5, 4.4, 0.1)
seq_total <- round(lseq(10^5, 10^9, 40))
parameter_values <- expand.grid(R0, seq_total)
ttd <- rpois(nrow(parameter_values), 10)
df <- data.frame(R0 = parameter_values$Var1, seq_total = parameter_values$Var2, time_to_detection = ttd)

scaleFUN <- function(x) sprintf("%.0f", x)


# Give extreme colors:
x <- ggplot(df, aes(x = R0, y = seq_total, fill = time_to_detection)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, option = "magma", limits = c(0, max(df$time_to_detection))) +
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

avg_ttd_by_R0 <- df %>%
  group_by(R0) %>%
  summarise(time_to_detection = mean(ttd, na.rm = TRUE))
y <- ggplot(avg_ttd_by_R0, aes(x = R0, y = time_to_detection)) +
  geom_line() +
  scale_x_continuous(labels = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scaleFUN) +
  labs(x = "", y = "TTD") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,-15,0),
        axis.title.y = element_text(angle = 360, hjust = 1, vjust = 0.5)) 

avg_ttd_by_seq <- df %>%
  group_by(seq_total) %>%
  summarise(time_to_detection = mean(ttd, na.rm = TRUE))
z <- ggplot(avg_ttd_by_seq, aes(x = seq_total, y = time_to_detection)) +
  geom_line() +
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

layout <- "
AA#
BB#
CCD
CCD
"
y + plot_spacer() + x + z + 
  plot_layout(design = layout, heights = c(1, -0.25, 5), widths = c(1, 5))


