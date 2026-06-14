# THIS CODE IS REWRITTEN WITH AI 
# Sequential Sensitivity Analysis (Parameter Sweep) for the Schelling Model
# Import parameters and functions
source("shelling_parameters.r")
source("shelling_functions.r")
source("shelling_time_simulation.r")

#+ output directory 
out_dir <- "../output"
dir.create(out_dir, showWarnings = FALSE)

# Sequential execution setup
cat("Esecuzione sequenziale per evitare deadlock su macOS...\n")

# Fast preview toggle
PREVIEW_MODE <- FALSE # Set to FALSE for full execution

#+ set up the sweeps 
topologies <- c("lattice", "smallworld", "scalefree", "erdos")
thresholds <- seq(0.1, 0.9, by = 0.025)
replications <- if (PREVIEW_MODE) 1 else 10
max_steps <- 10000

#+ no printed pdf
pdf(NULL)

# ==============================================================================
# 1. SWEEP PRIMARIO: CONFRONTO TOPOLOGIE
# ==============================================================================
cat("Inizio dello sweep primario (Confronto Topologie)...\n")

grid_primary <- expand.grid(
  topology = topologies,
  F_th = thresholds,
  replica = 1:replications,
  stringsAsFactors = FALSE
)

run_primary_single <- function(i) {
  row <- grid_primary[i, ]
  topo <- row$topology
  f_val <- row$F_th
  rep <- row$replica
  
  # Set threshold globally for the worker process
  F_th <<- f_val
  
  # 1. Generazione rete
  g <- network_topology(topo, N, l, p_small_world = 0.5)
  
  # 2. Inizializzazione attributi
  att <- c(rep(0, zero), rep(1, plus), rep(-1, minus)) |> sample()
  V(g)$att <- att
  
  # 3. Matrice adiacenza
  A <<- as_adjacency_matrix(g)
  
  # 4. Esecuzione
  g_final <- time_simulation(g, max_steps)
  
  # 5. Estrazione metriche
  conv_t <- g_final$conv_time
  int_dens <- interface_density(g_final)
  
  # Cluster
  comp <- sub_clusters(g_final, plot_flag = FALSE)
  n_clusters <- comp$no
  m_cluster_size <- mean(comp$csize)
  if (is.nan(m_cluster_size)) m_cluster_size <- 0
  
  if (conv_t < max_steps) {
    state <- if (int_dens < 0.2) "Segregated" else "Mixed"
  } else {
    state <- "Frozen"
  }
  
  return(data.frame(
    topology = topo,
    F_th = f_val,
    replica = rep,
    conv_time = conv_t,
    interface_density = int_dens,
    num_clusters = n_clusters,
    mean_cluster_size = m_cluster_size / N,
    state = state,
    stringsAsFactors = FALSE
  ))
}

# Run sequentially
primary_results_list <- lapply(1:nrow(grid_primary), run_primary_single)
raw_results <- do.call(rbind, primary_results_list)

write.csv(raw_results, file = file.path(out_dir, "sensitivity_analysis_raw.csv"), row.names = FALSE)

# Calcola e salva statistiche riassuntive
summary_results <- data.frame(
  topology = character(),
  F_th = numeric(),
  mean_conv_time = numeric(),
  sd_conv_time = numeric(),
  mean_interface_density = numeric(),
  sd_interface_density = numeric(),
  mean_num_clusters = numeric(),
  sd_num_clusters = numeric(),
  mean_cluster_size = numeric(),
  sd_cluster_size = numeric(),
  stringsAsFactors = FALSE
)

for (topo in topologies) {
  for (f_val in thresholds) {
    subset_data <- subset(raw_results, topology == topo & F_th == f_val)
    summary_results <- rbind(summary_results, data.frame(
      topology = topo,
      F_th = f_val,
      mean_conv_time = mean(subset_data$conv_time),
      sd_conv_time = ifelse(length(subset_data$conv_time) > 1, sd(subset_data$conv_time), 0),
      mean_interface_density = mean(subset_data$interface_density),
      sd_interface_density = ifelse(length(subset_data$interface_density) > 1, sd(subset_data$interface_density), 0),
      mean_num_clusters = mean(subset_data$num_clusters),
      sd_num_clusters = ifelse(length(subset_data$num_clusters) > 1, sd(subset_data$num_clusters), 0),
      mean_cluster_size = mean(subset_data$mean_cluster_size),
      sd_cluster_size = ifelse(length(subset_data$mean_cluster_size) > 1, sd(subset_data$mean_cluster_size), 0),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(summary_results, file = file.path(out_dir, "sensitivity_analysis_summary.csv"), row.names = FALSE)


# ==============================================================================
# 2. SWEEP SECONDARIO: SMALL-WORLD REWIRING
# ==============================================================================
cat("\nInizio dello sweep secondario (Small-World Rewiring)...\n")
sw_p_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)

grid_sw <- expand.grid(
  p_rewire = sw_p_values,
  F_th = thresholds,
  replica = 1:replications,
  stringsAsFactors = FALSE
)

run_sw_single <- function(i) {
  row <- grid_sw[i, ]
  p_val <- row$p_rewire
  f_val <- row$F_th
  rep <- row$replica
  
  F_th <<- f_val
  
  g <- network_topology("smallworld", N, l, p_small_world = p_val)
  att <- c(rep(0, zero), rep(1, plus), rep(-1, minus)) |> sample()
  V(g)$att <- att
  A <<- as_adjacency_matrix(g)
  
  g_final <- time_simulation(g, max_steps)
  
  conv_t <- g_final$conv_time
  int_dens <- interface_density(g_final)
  comp <- sub_clusters(g_final, plot_flag = FALSE)
  n_clusters <- comp$no
  m_cluster_size <- mean(comp$csize)
  if (is.nan(m_cluster_size)) m_cluster_size <- 0
  
  if (conv_t < max_steps) {
    state <- if (int_dens < 0.2) "Segregated" else "Mixed"
  } else {
    state <- "Frozen"
  }
  
  return(data.frame(
    F_th = f_val,
    p_rewire = p_val,
    replica = rep,
    conv_time = conv_t,
    interface_density = int_dens,
    num_clusters = n_clusters,
    mean_cluster_size = m_cluster_size / N,
    state = state,
    stringsAsFactors = FALSE
  ))
}

# Run sequentially
sw_results_list <- lapply(1:nrow(grid_sw), run_sw_single)
sw_raw_results <- do.call(rbind, sw_results_list)

write.csv(sw_raw_results, file = file.path(out_dir, "smallworld_sweep_raw.csv"), row.names = FALSE)

# Calcola e salva statistiche riassuntive
sw_summary_results <- data.frame(
  F_th = numeric(),
  p_rewire = numeric(),
  mean_conv_time = numeric(),
  sd_conv_time = numeric(),
  mean_interface_density = numeric(),
  sd_interface_density = numeric(),
  mean_num_clusters = numeric(),
  sd_num_clusters = numeric(),
  mean_cluster_size = numeric(),
  sd_cluster_size = numeric(),
  stringsAsFactors = FALSE
)

for (p_val in sw_p_values) {
  for (f_val in thresholds) {
    subset_data <- subset(sw_raw_results, p_rewire == p_val & F_th == f_val)
    sw_summary_results <- rbind(sw_summary_results, data.frame(
      F_th = f_val,
      p_rewire = p_val,
      mean_conv_time = mean(subset_data$conv_time),
      sd_conv_time = ifelse(length(subset_data$conv_time) > 1, sd(subset_data$conv_time), 0),
      mean_interface_density = mean(subset_data$interface_density),
      sd_interface_density = ifelse(length(subset_data$interface_density) > 1, sd(subset_data$interface_density), 0),
      mean_num_clusters = mean(subset_data$num_clusters),
      sd_num_clusters = ifelse(length(subset_data$num_clusters) > 1, sd(subset_data$num_clusters), 0),
      mean_cluster_size = mean(subset_data$mean_cluster_size),
      sd_cluster_size = ifelse(length(subset_data$mean_cluster_size) > 1, sd(subset_data$mean_cluster_size), 0),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(sw_summary_results, file = file.path(out_dir, "smallworld_sweep_summary.csv"), row.names = FALSE)


# ==============================================================================
# 3. SWEEP DI VACANCY (HOLE DENSITY)
# ==============================================================================
cat("\nInizio dello sweep di vacancy (Hole Density)...\n")
hole_densities <- seq(0.05, 0.50, by = 0.05)
vacancy_topologies <- c("lattice", "smallworld", "scalefree", "erdos")
vacancy_replications <- if (PREVIEW_MODE) 1 else 10
fixed_F_th <- 0.5

grid_vac <- expand.grid(
  topology = vacancy_topologies,
  vacancy_density = hole_densities,
  replica = 1:vacancy_replications,
  stringsAsFactors = FALSE
)

run_vac_single <- function(i) {
  row <- grid_vac[i, ]
  topo <- row$topology
  rho_z <- row$vacancy_density
  rep <- row$replica
  
  F_th <<- fixed_F_th
  
  z_nodes <- round(rho_z * N)
  p_nodes <- round((N - z_nodes) / 2)
  m_nodes <- N - z_nodes - p_nodes
  
  g <- network_topology(topo, N, l, p_small_world = 0.5)
  att <- c(rep(0, z_nodes), rep(1, p_nodes), rep(-1, m_nodes)) |> sample()
  V(g)$att <- att
  A <<- as_adjacency_matrix(g)
  
  g_final <- time_simulation(g, max_steps)
  
  conv_t <- g_final$conv_time
  int_dens <- interface_density(g_final)
  comp <- sub_clusters(g_final, plot_flag = FALSE)
  n_clusters <- comp$no
  m_cluster_size <- mean(comp$csize)
  if (is.nan(m_cluster_size)) m_cluster_size <- 0
  
  return(data.frame(
    topology = topo,
    vacancy_density = rho_z,
    replica = rep,
    conv_time = conv_t,
    interface_density = int_dens,
    num_clusters = n_clusters,
    mean_cluster_size = m_cluster_size / N,
    stringsAsFactors = FALSE
  ))
}

# Run sequentially
vac_results_list <- lapply(1:nrow(grid_vac), run_vac_single)
vac_raw_results <- do.call(rbind, vac_results_list)

write.csv(vac_raw_results, file = file.path(out_dir, "vacancy_sweep_raw.csv"), row.names = FALSE)

# Calcola e salva statistiche riassuntive
vac_summary_results <- data.frame(
  topology = character(),
  vacancy_density = numeric(),
  mean_conv_time = numeric(),
  sd_conv_time = numeric(),
  mean_interface_density = numeric(),
  sd_interface_density = numeric(),
  mean_num_clusters = numeric(),
  sd_num_clusters = numeric(),
  mean_cluster_size = numeric(),
  sd_cluster_size = numeric(),
  stringsAsFactors = FALSE
)

for (topo in vacancy_topologies) {
  for (rho_z in hole_densities) {
    subset_data <- subset(vac_raw_results, topology == topo & vacancy_density == rho_z)
    vac_summary_results <- rbind(vac_summary_results, data.frame(
      topology = topo,
      vacancy_density = rho_z,
      mean_conv_time = mean(subset_data$conv_time),
      sd_conv_time = ifelse(length(subset_data$conv_time) > 1, sd(subset_data$conv_time), 0),
      mean_interface_density = mean(subset_data$interface_density),
      sd_interface_density = ifelse(length(subset_data$interface_density) > 1, sd(subset_data$interface_density), 0),
      mean_num_clusters = mean(subset_data$num_clusters),
      sd_num_clusters = ifelse(length(subset_data$num_clusters) > 1, sd(subset_data$num_clusters), 0),
      mean_cluster_size = mean(subset_data$mean_cluster_size),
      sd_cluster_size = ifelse(length(subset_data$mean_cluster_size) > 1, sd(subset_data$mean_cluster_size), 0),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(vac_summary_results, file = file.path(out_dir, "vacancy_sweep_summary.csv"), row.names = FALSE)

# Chiudiamo il device grafico nullo
dev.off()

cat("Tutte le simulazioni completate con successo! File salvati in:", out_dir, "\n")
