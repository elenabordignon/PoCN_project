# Premium Plotting Script for Schelling Model Phase Transition & Susceptibility
# Generates scientific quality plots with clean modern aesthetics.
# THIS CODE IS WRITTEN WITH AI TOOLS AND THEN CHECKED 

library(ggplot2)
library(dplyr)

# Output folder
out_dir <- "../output"
summary_file <- file.path(out_dir, "phase_transition_summary.csv")

if (!file.exists(summary_file)) {
  stop("Errore: Esegui prima shelling_phase_transition.r per generare i dati.")
}

# Load data
df <- read.csv(summary_file)

# Convert p_rewire to factor for clean plot legends
df$p_rewire_factor <- as.factor(df$p_rewire)

# Compute standard deviation for error bars
df$sd_seg_index <- sqrt(df$var_seg_index)

# Setup dynamic premium color palette based on unique p values
p_unique <- sort(unique(df$p_rewire))
premium_colors <- c(
  "#2c3e50", # Dark Slate
  "#2980b9", # Classic Blue
  "#27ae60", # Emerald Green
  "#f39c12", # Warm Orange
  "#e74c3c", # Vibrant Red
  "#9b59b6", # Amethyst Purple
  "#1abc9c"  # Turquoise
)
# Make sure we have enough colors, otherwise repeat
colors_to_use <- premium_colors[1:length(p_unique)]
color_map <- setNames(colors_to_use, as.character(p_unique))

# Define premium theme matching academic publication style
premium_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#2c3e50", margin = margin(b = 6), hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = "#7f8c8d", margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 10, color = "#2c3e50"),
    axis.text = element_text(size = 8, color = "#34495e"),
    legend.title = element_text(face = "bold", size = 9, color = "#2c3e50"),
    legend.text = element_text(size = 8, color = "#34495e"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "#ecf0f1", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# ==========================================================
# PLOT 1: SEGREGATION INDEX VS TOLERANCE (F)
# ==========================================================
cat("Generating Segregation Phase Transition Plot...\n")

p1 <- ggplot(df, aes(x = F_th, y = mean_seg_index, color = p_rewire_factor, group = p_rewire_factor)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_seg_index - sd_seg_index), 
                    ymax = pmin(1, mean_seg_index + sd_seg_index)), 
                width = 0.005, alpha = 0.4) +
  scale_color_manual(values = color_map, name = "Rewire Prob. (p)") +
  labs(
    title = "Segregation Phase Transition Curve",
    subtitle = "Steady-state segregation index (s = 1 - 2n) vs Tolerance (F). Mean ± SD",
    x = "Tolerance Threshold (F)",
    y = "Segregation Index (s)"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_segregation_index.png"), plot = p1, width = 7, height = 5, dpi = 300)

# ==========================================================
# PLOT 2: SUSCEPTIBILITY VS TOLERANCE (F)
# ==========================================================
cat("Generating Susceptibility Curve Plot...\n")

p2 <- ggplot(df, aes(x = F_th, y = susceptibility, color = p_rewire_factor, group = p_rewire_factor)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_map, name = "Rewire Prob. (p)") +
  labs(
    title = "Susceptibility",
    x = "Tolerance Threshold (F)",
    y = "Susceptibility"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_susceptibility.png"), plot = p2, width = 7, height = 5, dpi = 300)

cat("Plots generated successfully in phase_transition_output/!\n")
