# Plotting Script for the Schelling Model Sensitivity Analysis (v2)
# Generates high-quality comparative plots and heatmaps in English.
# SWEEP ANALYSIS AND PLOTS ARE GENERATED WITH AI AND THEN CHECKED
library(ggplot2)
library(dplyr)

# Set input and output directories
out_dir <- "../output"
summary_file <- file.path(out_dir, "sensitivity_analysis_summary.csv")
raw_file <- file.path(out_dir, "sensitivity_analysis_raw.csv")

sw_summary_file <- file.path(out_dir, "smallworld_sweep_summary.csv")
sw_raw_file <- file.path(out_dir, "smallworld_sweep_raw.csv")

if (!file.exists(summary_file) || !file.exists(sw_summary_file)) {
  stop("Errore: Esegui prima sensitivity_analysis.r per generare i dati.")
}

# Load datasets
df <- read.csv(summary_file)
raw_df <- read.csv(raw_file)

df_sw <- read.csv(sw_summary_file)
raw_df_sw <- read.csv(sw_raw_file)

# Color palettes
color_palette <- c(
  "lattice" = "#2c3e50",     # Sleek dark blue
  "smallworld" = "#e74c3c",  # Coral red
  "scalefree" = "#2ecc71",   # Emerald green
  "erdos" = "#f1c40f"        # Mustard yellow
)

# Colors for states in the heatmap
state_colors <- c(
  "Frozen" = "#2980b9",      # Soft Blue
  "Segregated" = "#27ae60",  # Soft Green
  "Mixed" = "#c0392b"        # Soft Red
)

# Premium Theme (GGPlot)
premium_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#2c3e50", margin = margin(b = 8)),
    plot.subtitle = element_text(size = 9, color = "#7f8c8d", margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 10, color = "#2c3e50"),
    axis.text = element_text(size = 8, color = "#34495e"),
    legend.title = element_text(face = "bold", size = 9, color = "#2c3e50"),
    legend.text = element_text(size = 8, color = "#34495e"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "#ecf0f1"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# ==========================================================
# PART 1: COMPARATIVE PLOTS FOR TOPOLOGIES (English, <= 4 words labels)
# ==========================================================

# 1.1 Convergence Time vs F_th
p1 <- ggplot(df, aes(x = F_th, y = mean_conv_time, color = topology, group = topology)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_conv_time - sd_conv_time), 
                    ymax = pmin(10000, mean_conv_time + sd_conv_time)), 
                width = 0.005, alpha = 0.6) +
  scale_color_manual(values = color_palette, name = "Topology") +
  labs(
    title = "Convergence Time Comparison",
    subtitle = "Mean ± SD over 10 runs per topology",
    x = "Tolerance Threshold (F)",
    y = "Convergence Time (Steps)"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_convergence_time.png"), plot = p1, width = 7, height = 5, dpi = 300)

# 1.2 Segregation Index vs F_th
p2 <- ggplot(df, aes(x = F_th, y = 1 - mean_interface_density, color = topology, group = topology)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, 1 - mean_interface_density - sd_interface_density), 
                    ymax = pmin(1, 1 - mean_interface_density + sd_interface_density)), 
                width = 0.005, alpha = 0.6) +
  scale_color_manual(values = color_palette, name = "Topology") +
  labs(
    title = "Phase Transition Curves (Segregation)",
    subtitle = "Fraction of same-type links at steady state",
    x = "Tolerance Threshold (F)",
    y = "Segregation Index (S)"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_interface_density.png"), plot = p2, width = 7, height = 5, dpi = 300)

# 1.3 Number of Clusters vs F_th
p3 <- ggplot(df, aes(x = F_th, y = mean_num_clusters, color = topology, group = topology)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_num_clusters - sd_num_clusters), 
                    ymax = mean_num_clusters + sd_num_clusters), 
                width = 0.005, alpha = 0.6) +
  scale_color_manual(values = color_palette, name = "Topology") +
  labs(
    title = "Number of Clusters",
    subtitle = "Count of same-type connected components",
    x = "Tolerance Threshold (F)",
    y = "Number of Clusters"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_number_of_clusters.png"), plot = p3, width = 7, height = 5, dpi = 300)

# 1.4 Normalized Cluster Size vs F_th
p4 <- ggplot(df, aes(x = F_th, y = mean_cluster_size, color = topology, group = topology)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_cluster_size - sd_cluster_size), 
                    ymax = pmin(1, mean_cluster_size + sd_cluster_size)), 
                width = 0.005, alpha = 0.6) +
  scale_color_manual(values = color_palette, name = "Topology") +
  labs(
    title = "Normalized Cluster Size",
    subtitle = "Average cluster size divided by network nodes (N=400)",
    x = "Tolerance Threshold (F)",
    y = "Normalized Cluster Size"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_average_cluster_size.png"), plot = p4, width = 7, height = 5, dpi = 300)


# ==========================================================
# PART 2: SMALL-WORLD SWEEP PLOTS (F_th vs. p_rewire)
# ==========================================================

# Convert p_rewire to factor for clean plot legends
df_sw$p_rewire_factor <- as.factor(df_sw$p_rewire)

# 2.1 SW Convergence Time
sw_p1 <- ggplot(df_sw, aes(x = F_th, y = mean_conv_time, color = p_rewire_factor, group = p_rewire_factor)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_conv_time - sd_conv_time), 
                    ymax = pmin(10000, mean_conv_time + sd_conv_time)), 
                width = 0.005, alpha = 0.5) +
  scale_color_brewer(palette = "YlOrRd", name = "Rewire Prob. (p)") +
  labs(
    title = "Small-World Convergence Time",
    subtitle = "Impact of rewiring probability (p) on convergence speed",
    x = "Tolerance Threshold (F)",
    y = "Convergence Time (Steps)"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_sw_convergence_time.png"), plot = sw_p1, width = 7, height = 5, dpi = 300)

# 2.2 SW Segregation Index
sw_p2 <- ggplot(df_sw, aes(x = F_th, y = 1 - mean_interface_density, color = p_rewire_factor, group = p_rewire_factor)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, 1 - mean_interface_density - sd_interface_density), 
                    ymax = pmin(1, 1 - mean_interface_density + sd_interface_density)), 
                width = 0.005, alpha = 0.5) +
  scale_color_brewer(palette = "YlOrRd", name = "Rewire Prob. (p)") +
  labs(
    title = "Small-World Segregation Index",
    subtitle = "Impact of rewiring probability (p) on segregation phase transitions",
    x = "Tolerance Threshold (F)",
    y = "Segregation Index (S)"
  ) +
  premium_theme

ggsave(file.path(out_dir, "plot_sw_interface_density.png"), plot = sw_p2, width = 7, height = 5, dpi = 300)


# ==========================================================
# PART 3: PHASE DIAGRAM HEATMAP (Small-World F vs. p)
# ==========================================================

# Helper function to find the mode (most common state)
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Group raw data to find the consensus state for each (F, p) grid cell
sw_states <- raw_df_sw %>%
  group_by(F_th, p_rewire) %>%
  summarise(state = get_mode(state), .groups = 'drop')

# Make the state factor with ordered levels
sw_states$state <- factor(sw_states$state, levels = c("Mixed", "Segregated", "Frozen"))

# Generate the Heatmap
p_heatmap <- ggplot(sw_states, aes(x = F_th, y = p_rewire, fill = state)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = state_colors, name = "Phase State") +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_y_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  labs(
    title = "Small-World Phase Diagram",
    subtitle = "Steady-state regimes as a function of tolerance (F) and rewiring (p)",
    x = "Tolerance (F)",
    y = "Rewiring Probability (p)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#2c3e50", margin = margin(b = 8)),
    plot.subtitle = element_text(size = 9, color = "#7f8c8d", margin = margin(b = 12)),
    axis.title = element_text(face = "bold", size = 10, color = "#2c3e50"),
    axis.text = element_text(size = 8, color = "#34495e"),
    legend.title = element_text(face = "bold", size = 9, color = "#2c3e50"),
    legend.text = element_text(size = 8, color = "#34495e"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(out_dir, "plot_sw_phase_heatmap.png"), plot = p_heatmap, width = 7, height = 5, dpi = 300)

# ==========================================================
# PART 4: VACANCY SWEEP PLOT (Segregation vs. Hole Density)
# ==========================================================
vac_summary_file <- file.path(out_dir, "vacancy_sweep_summary.csv")

if (file.exists(vac_summary_file)) {
  df_vac <- read.csv(vac_summary_file)
  
  p_vac <- ggplot(df_vac, aes(x = vacancy_density, y = 1 - mean_interface_density, color = topology, group = topology)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = pmax(0, 1 - mean_interface_density - sd_interface_density), 
                      ymax = pmin(1, 1 - mean_interface_density + sd_interface_density)), 
                  width = 0.005, alpha = 0.6) +
    scale_color_manual(values = color_palette, name = "Topology") +
    labs(
      title = "Segregation vs. Vacancy Density",
      subtitle = "Steady-state segregation index for F = 0.5 over 10 runs",
      x = "Vacancy Density (Hole Ratio)",
      y = "Segregation Index (S)"
    ) +
    premium_theme
  
  ggsave(file.path(out_dir, "plot_segregation_vs_holes.png"), plot = p_vac, width = 7, height = 5, dpi = 300)
}

cat("All plots generated successfully in shelling_sensitvity_v2/!\n")
