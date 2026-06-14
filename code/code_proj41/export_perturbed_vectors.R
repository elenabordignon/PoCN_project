# export_perturbed_vectors.R
# This script loads the processed analysis results and exports the perturbed vectors
# (delta_x) for all configurations into a single tidy CSV file.

# Ensure output directory exists
dir.create("output_data", showWarnings = FALSE)

# Load the processed RData
rdata_path <- "data/processed/analysis_results.RData"
if (!file.exists(rdata_path)) {
  stop(paste("File not found:", rdata_path))
}
message("Loading processed results from: ", rdata_path)
load(rdata_path)

# Helper function to convert delta_x_list to a data frame
convert_delta_list <- function(delta_list, time_steps, config_name, amplitude, window) {
  df_list <- list()
  idx <- 1
  
  for (pert_sp in names(delta_list)) {
    mat <- delta_list[[pert_sp]]
    for (resp_sp in colnames(mat)) {
      df_list[[idx]] <- data.frame(
        configuration = config_name,
        amplitude = amplitude,
        window = window,
        perturbed_species = pert_sp,
        response_species = resp_sp,
        time = time_steps,
        delta_x = mat[, resp_sp],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, df_list)
}

message("Processing Configuration A (Baseline: Amplitude = 1.0, Window = 200)...")
df_A <- convert_delta_list(delta_x_list_A, time_steps_A, "A_baseline", 1.0, 200)

message("Processing Configuration B (Amplitude = 0.5, Window = 200)...")
df_B <- convert_delta_list(delta_x_list_B, time_steps_B, "B_amplitude_0.5", 0.5, 200)

message("Processing Configuration C (Amplitude = 1.0, Window = 100)...")
df_C <- convert_delta_list(delta_x_list_C, time_steps_C, "C_window_100", 1.0, 100)

# Combine all configurations
combined_df <- rbind(df_A, df_B, df_C)

# Write to CSV
output_csv <- "output_data/perturbed_vectors.csv"
write.csv(combined_df, output_csv, row.names = FALSE)
message("Successfully exported perturbed vectors to: ", output_csv)
