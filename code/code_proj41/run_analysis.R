# run_analysis.R
# Analysis script to run the calculations and export computed matrices and objects,
# including robustness checks for different perturbation amplitudes and observation windows.

# Source the modular code
source("scripts/data_loader.R")
source("scripts/perturbation_analysis.R")
source("scripts/helpers.R")

# Ensure required directories exist
ensure_dir("data/processed")
ensure_dir("outputs/results")

log_info("Starting TNM perturbation analysis...")

# ==============================================================================
# CONFIGURATION A: Amplitude = 1.0, Window = 200 (Baseline)
# ==============================================================================
log_info("Loading Configuration A: Amplitude = 1.0, Window = 200...")
ds_A <- load_dataset_from_dir("data/raw")
X_ref_A <- ds_A$X_ref
pert_list_A <- ds_A$pert_list
core_species <- ds_A$core_species
time_steps_A <- ds_A$time_steps

log_info("Computing delta_x and D for Config A...")
delta_x_list_A <- compute_delta_x(X_ref_A, pert_list_A)
D_A <- compute_integrated_distance_matrix(delta_x_list_A)
hclust_A <- cluster_species(D_A, method = "average")

# ==============================================================================
# CONFIGURATION B: Amplitude = 0.5, Window = 200 (Robustness to amplitude)
# ==============================================================================
log_info("Loading Configuration B: Amplitude = 0.5, Window = 200...")
ds_B <- load_dataset_from_dir("data/raw_amp0.5")
X_ref_B <- ds_B$X_ref
pert_list_B <- ds_B$pert_list
time_steps_B <- ds_B$time_steps

log_info("Computing delta_x and D for Config B...")
delta_x_list_B <- compute_delta_x(X_ref_B, pert_list_B)
D_B <- compute_integrated_distance_matrix(delta_x_list_B)
hclust_B <- cluster_species(D_B, method = "average")

# ==============================================================================
# CONFIGURATION C: Amplitude = 1.0, Window = 100 (Robustness to observation window)
# ==============================================================================
log_info("Creating Configuration C: Amplitude = 1.0, Window = 100 (subsetting)...")
# Select time steps up to generation 4100 (first 100 generations of observation)
# Since the simulation starts at 4000, we want t <= 4100
t_subset <- which(time_steps_A <= 4100)

X_ref_C <- X_ref_A[t_subset, , drop = FALSE]
pert_list_C <- list()
for (name in names(pert_list_A)) {
    pert_list_C[[name]] <- pert_list_A[[name]][t_subset, , drop = FALSE]
}
time_steps_C <- time_steps_A[t_subset]

log_info("Computing delta_x and D for Config C...")
delta_x_list_C <- compute_delta_x(X_ref_C, pert_list_C)
D_C <- compute_integrated_distance_matrix(delta_x_list_C)
hclust_C <- cluster_species(D_C, method = "average")

# ==============================================================================
# ABUNDANCE CORRELATION BASELINE
# ==============================================================================
log_info("Computing abundance correlation baseline matrix...")
D_corr <- compute_correlation_baseline(X_ref_A)
hclust_baseline <- cluster_species(D_corr, method = "average")

# ==============================================================================
# ROBUSTNESS METRICS (Correlation of distance matrices)
# ==============================================================================
dist_vec_A <- D_A[lower.tri(D_A)]
dist_vec_B <- D_B[lower.tri(D_B)]
dist_vec_C <- D_C[lower.tri(D_C)]

cor_A_B <- cor(dist_vec_A, dist_vec_B, method = "pearson")
cor_A_C <- cor(dist_vec_A, dist_vec_C, method = "pearson")

log_info("Robustness - Pearson Correlation between D_A (Amp=1.0) and D_B (Amp=0.5):", cor_A_B)
log_info("Robustness - Pearson Correlation between D_A (Win=200) and D_C (Win=100):", cor_A_C)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
# Save matrices to CSV
save_matrix_csv(D_A, "outputs/results/perturbation_distance_matrix_A_1.0_200.csv")
save_matrix_csv(D_B, "outputs/results/perturbation_distance_matrix_B_0.5_200.csv")
save_matrix_csv(D_C, "outputs/results/perturbation_distance_matrix_C_1.0_100.csv")
save_matrix_csv(D_corr, "outputs/results/correlation_baseline_matrix.csv")

# Save primary distance matrix as perturbation_distance_matrix.csv for backward compatibility
save_matrix_csv(D_A, "outputs/results/perturbation_distance_matrix.csv")

# Save processed objects to RData
processed_data_path <- "data/processed/analysis_results.RData"
save(
    X_ref_A, delta_x_list_A, D_A, hclust_A, time_steps_A,
    X_ref_B, delta_x_list_B, D_B, hclust_B, time_steps_B,
    X_ref_C, delta_x_list_C, D_C, hclust_C, time_steps_C,
    D_corr, hclust_baseline,
    core_species, cor_A_B, cor_A_C,
    file = processed_data_path
)
log_info("Saved processed RData to:", processed_data_path)

# Export a simple summary text file of robustness results
summary_file <- "outputs/results/robustness_summary.txt"
cat(sprintf("Robustness analysis results:\n"), file = summary_file)
cat(sprintf("- Pearson correlation D(Amp=1.0, Win=200) vs D(Amp=0.5, Win=200): %.5f\n", cor_A_B), file = summary_file, append = TRUE)
cat(sprintf("- Pearson correlation D(Amp=1.0, Win=200) vs D(Amp=1.0, Win=100): %.5f\n", cor_A_C), file = summary_file, append = TRUE)
log_info("Saved robustness summary to:", summary_file)

log_info("Analysis and robustness checks successfully completed!")
