# High-resolution phase transition and susceptibility analysis
# for the Schelling model on Small-World networks.
# THIS CODE IS MODIFIED WITH AI 

source("shelling_parameters.r")
source("shelling_functions.r")
source("shelling_time_simulation.r")

# Output folder
out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE)

# Preview mode configuration (toggle to FALSE for final run)
PREVIEW_MODE <- FALSE 

if (PREVIEW_MODE) {
  p_values <- c(0.0, 0.1, 0.5, 1.0)
  thresholds <- seq(0.1, 0.8, by = 0.05)
  replications <- 3
  max_steps <- 5000
} else {
  p_values <- c(0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0)
  thresholds <- seq(0.1, 0.85, by = 0.02)
  replications <- 30
  max_steps <- 10000
}

# Turn off R plotting device during simulation
pdf(NULL)

cat(sprintf("Starting simulation sweep (PREVIEW_MODE = %s)...\n", PREVIEW_MODE))

grid_params <- expand.grid(
  p_rewire = p_values,
  F_th = thresholds,
  replica = 1:replications,
  stringsAsFactors = FALSE
)

run_single_sim <- function(i) {
  row <- grid_params[i, ]
  p_val <- row$p_rewire
  f_val <- row$F_th
  rep <- row$replica
  
  # Set F_th globally for the workers/functions
  F_th <<- f_val
  
  # Generate network
  g <- network_topology("smallworld", N, l, p_small_world = p_val)
  
  # Initialize attributes
  att <- c(rep(0, zero), rep(1, plus), rep(-1, minus)) |> sample()
  V(g)$att <- att
  
  # Setup global adjacency matrix for performance
  A <<- as_adjacency_matrix(g)
  
  # Run simulation
  g_final <- time_simulation(g, max_steps)
  
  # Compute metrics
  int_dens <- interface_density(g_final)
  seg_index <- 1 - 2 * int_dens
  conv_t <- g_final$conv_time
  
  # Return data frame
  return(data.frame(
    p_rewire = p_val,
    F_th = f_val,
    replica = rep,
    conv_time = conv_t,
    interface_density = int_dens,
    segregation_index = seg_index,
    stringsAsFactors = FALSE
  ))
}

# Run sequentially
results_list <- lapply(1:nrow(grid_params), run_single_sim)
raw_results <- do.call(rbind, results_list)

# Write raw results
write.csv(raw_results, file = file.path(out_dir, "phase_transition_raw.csv"), row.names = FALSE)

# Aggregate results (mean and variance)
summary_results <- raw_results %>%
  group_by(p_rewire, F_th) %>%
  summarise(
    mean_seg_index = mean(segregation_index),
    var_seg_index = var(segregation_index),
    mean_interface_density = mean(interface_density),
    var_interface_density = var(interface_density),
    mean_conv_time = mean(conv_time),
    .groups = "drop"
  )

# Calculate susceptibility chi = N * var_seg_index
summary_results$susceptibility <- N * summary_results$var_seg_index

# Write summary results
write.csv(summary_results, file = file.path(out_dir, "phase_transition_summary.csv"), row.names = FALSE)

# Close graphic device
if (dev.cur() > 1) dev.off()

cat(sprintf("Sweep completed! Results saved to %s/\n", out_dir))
