# compute_time_dependent_distances.R
# This script calculates the time-dependent distance matrix d_t(i,j) for each time step t
# and exports the results as a tidy CSV.

# Ensure output directory exists
dir.create("output_data", showWarnings = FALSE)

# Load the processed RData
rdata_path <- "data/processed/analysis_results.RData"
if (!file.exists(rdata_path)) {
  stop(paste("File not found:", rdata_path))
}
message("Loading processed results from: ", rdata_path)
load(rdata_path)

# Function to compute time-dependent distance data frame
compute_dt_df <- function(delta_list, time_steps, config_name) {
  species_names <- names(delta_list)
  M <- length(species_names)
  T_len <- length(time_steps)
  
  # Allocate vectors for performance
  total_rows <- M * M * T_len
  out_config <- rep(config_name, total_rows)
  out_time <- rep(0, total_rows)
  out_sp_i <- rep("", total_rows)
  out_sp_j <- rep("", total_rows)
  out_dist <- rep(0.0, total_rows)
  
  idx <- 1
  for (t_idx in 1:T_len) {
    t_val <- time_steps[t_idx]
    
    # Pre-extract the vectors for all species at this time step
    t_vectors <- lapply(species_names, function(sp) delta_list[[sp]][t_idx, ])
    names(t_vectors) <- species_names
    
    for (i in 1:M) {
      name_i <- species_names[i]
      vec_i <- t_vectors[[name_i]]
      
      for (j in 1:M) {
        name_j <- species_names[j]
        vec_j <- t_vectors[[name_j]]
        
        out_time[idx] <- t_val
        out_sp_i[idx] <- name_i
        out_sp_j[idx] <- name_j
        out_dist[idx] <- sqrt(sum((vec_i - vec_j)^2))
        idx <- idx + 1
      }
    }
  }
  
  data.frame(
    configuration = out_config,
    time = out_time,
    species_i = out_sp_i,
    species_j = out_sp_j,
    distance = out_dist,
    stringsAsFactors = FALSE
  )
}

message("Computing time-dependent distances for Configuration A (Baseline)...")
df_A <- compute_dt_df(delta_x_list_A, time_steps_A, "A_baseline")

message("Computing time-dependent distances for Configuration B (Amplitude = 0.5)...")
df_B <- compute_dt_df(delta_x_list_B, time_steps_B, "B_amplitude_0.5")

message("Computing time-dependent distances for Configuration C (Window = 100)...")
df_C <- compute_dt_df(delta_x_list_C, time_steps_C, "C_window_100")

# Combine all configurations
combined_df <- rbind(df_A, df_B, df_C)

# Write time-dependent distances to CSV
output_csv <- "output_data/time_dependent_distances.csv"
write.csv(combined_df, output_csv, row.names = FALSE)
message("Successfully exported time-dependent distances to: ", output_csv)

# ==============================================================================
# COMPUTE TIME-AVERAGED DISTANCE MATRICES
# ==============================================================================
message("Computing time-averaged distance matrices...")

# Calculate the mean over time for each pair of species in each configuration
mean_df <- aggregate(distance ~ configuration + species_i + species_j, 
                     data = combined_df, FUN = mean)

# Save the tidy long-format table
output_mean_csv <- "output_data/time_averaged_distances.csv"
write.csv(mean_df, output_mean_csv, row.names = FALSE)
message("Successfully exported time-averaged distances (long format) to: ", output_mean_csv)

# Function to pivot long format to wide matrix format (12x12)
export_wide_matrix <- function(df, config_name, filepath) {
  sub_df <- df[df$configuration == config_name, ]
  
  # Get unique species list (should be 12 core species)
  sps <- unique(sub_df$species_i)
  M <- length(sps)
  
  # Initialize matrix
  mat <- matrix(0, nrow = M, ncol = M)
  rownames(mat) <- sps
  colnames(mat) <- sps
  
  # Fill matrix
  for (row in 1:nrow(sub_df)) {
    mat[sub_df$species_i[row], sub_df$species_j[row]] <- sub_df$distance[row]
  }
  
  write.csv(mat, filepath)
}

# Export 12x12 matrices for each configuration
export_wide_matrix(mean_df, "A_baseline", "output_data/time_averaged_matrix_A_baseline.csv")
export_wide_matrix(mean_df, "B_amplitude_0.5", "output_data/time_averaged_matrix_B_amp0.5.csv")
export_wide_matrix(mean_df, "C_window_100", "output_data/time_averaged_matrix_C_win100.csv")

message("Successfully exported wide 12x12 matrices to output_data/")

