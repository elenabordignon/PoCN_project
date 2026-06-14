# calculate_grid_correlations_t4200.R
# Script to compute correlations and heatmaps for the 9-amplitude grid (0.2 to 1.0)
# starting at generation 4200, across two observation windows (50 and 200 generations).

# Source project modules
source("scripts/data_loader.R")
source("scripts/perturbation_analysis.R")
source("scripts/helpers.R")

ensure_dir("outputs/results_ensemble")
ensure_dir("outputs/figures_ensemble")

log_info("Starting 9-amplitude grid correlation analysis for t=4200...")

# 1. Define amplitudes and their raw data directories
amplitudes <- c("0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
dirs <- setNames(
    sprintf("data/raw_ensemble_amp%s_t4200", amplitudes),
    amplitudes
)

# Load raw datasets
log_info("Loading raw datasets for the 9-amplitude grid starting at 4200...")
datasets <- list()
for (amp_name in names(dirs)) {
    datasets[[amp_name]] <- load_dataset_from_dir(dirs[amp_name])
}

# Helper function to compute correlation matrix for a given window limit (max_t)
compute_grid_correlation_matrix <- function(max_t) {
    dist_matrices <- list()
    
    for (amp_name in names(dirs)) {
        ds <- datasets[[amp_name]]
        
        # Subset to the requested observation window (starts at 4200)
        t_subset <- which(ds$time_steps <= max_t)
        
        X_ref_sub <- ds$X_ref[t_subset, , drop = FALSE]
        pert_list_sub <- list()
        for (name in names(ds$pert_list)) {
            pert_list_sub[[name]] <- ds$pert_list[[name]][t_subset, , drop = FALSE]
        }
        
        delta_x_list <- compute_delta_x(X_ref_sub, pert_list_sub)
        D <- compute_integrated_distance_matrix(delta_x_list)
        dist_matrices[[amp_name]] <- D
    }
    
    n_amps <- length(amplitudes)
    cor_matrix <- matrix(1.0, nrow = n_amps, ncol = n_amps)
    rownames(cor_matrix) <- paste0("Amp ", names(dirs))
    colnames(cor_matrix) <- paste0("Amp ", names(dirs))
    
    for (i in 1:n_amps) {
        amp_i <- names(dirs)[i]
        D_i <- dist_matrices[[amp_i]]
        vec_i <- D_i[lower.tri(D_i)]
        
        for (j in i:n_amps) {
            amp_j <- names(dirs)[j]
            D_j <- dist_matrices[[amp_j]]
            vec_j <- D_j[lower.tri(D_j)]
            
            if (sd(vec_i) == 0 || sd(vec_j) == 0) {
                r_val <- 0.0
            } else {
                r_val <- cor(vec_i, vec_j, method = "pearson")
            }
            cor_matrix[i, j] <- r_val
            cor_matrix[j, i] <- r_val
        }
    }
    return(cor_matrix)
}

# 2. Compute correlation matrices for both windows
log_info("Computing correlation matrix for Window = 4200 to 4400 (200 generations)...")
cor_matrix_200 <- compute_grid_correlation_matrix(4400)

log_info("Computing correlation matrix for Window = 4200 to 4250 (50 generations)...")
cor_matrix_50 <- compute_grid_correlation_matrix(4250)

# 3. Save matrices to CSV
save_matrix_csv(cor_matrix_200, "outputs/results_ensemble/amplitude_grid_correlation_matrix_win200_t4200.csv")
save_matrix_csv(cor_matrix_50, "outputs/results_ensemble/amplitude_grid_correlation_matrix_win50_t4200.csv")

log_info("Grid Correlation Matrix - Delay + 200 generations:")
print(round(cor_matrix_200, 4))

log_info("Grid Correlation Matrix - Delay + 50 generations:")
print(round(cor_matrix_50, 4))

# 4. Helper function to plot heatmaps
plot_heatmap <- function(cor_matrix, file_path, title_text) {
    plot_success <- FALSE
    tryCatch({
        library(ggplot2)
        library(reshape2)
        
        melted_cor <- melt(cor_matrix)
        colnames(melted_cor) <- c("Var1", "Var2", "Correlation")
        melted_cor$Var1 <- factor(melted_cor$Var1, levels = rownames(cor_matrix))
        melted_cor$Var2 <- factor(melted_cor$Var2, levels = colnames(cor_matrix))
        
        p <- ggplot(melted_cor, aes(x = Var1, y = Var2, fill = Correlation)) +
            geom_tile(color = "white", lwd = 1.0, linetype = 1) +
            geom_text(aes(label = sprintf("%.2f", Correlation)), color = "white", size = 3.5, fontface = "bold") +
            scale_fill_gradientn(colors = c("#081d58", "#253494", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4"),
                                 limits = c(-0.3, 1), 
                                 name = "Pearson R") +
            theme_minimal(base_size = 14) +
            labs(
                title = title_text,
                x = "Perturbation Amplitude",
                y = "Perturbation Amplitude"
            ) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 11, margin = margin(b = 15)),
                axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
                axis.text.y = element_text(face = "bold"),
                panel.grid = element_blank(),
                legend.position = "right"
            )
        
        ggsave(file_path, plot = p, width = 8, height = 7, dpi = 150)
        plot_success <- TRUE
    }, error = function(e) {
        log_info(paste("ggplot2/reshape2 failed, falling back to base R...", e$message))
    })
    
    if (!plot_success) {
        png(file_path, width = 800, height = 700, res = 120)
        n_amps <- nrow(cor_matrix)
        op <- par(mar = c(5, 5, 4, 3) + 0.1)
        blue_palette <- colorRampPalette(c("#edf8b1", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"))(100)
        
        image(1:n_amps, 1:n_amps, t(cor_matrix[n_amps:1, ]), 
              col = blue_palette, 
              axes = FALSE, 
              xlab = "Perturbation Amplitude", 
              ylab = "Perturbation Amplitude",
              main = title_text)
        
        axis(1, at = 1:n_amps, labels = rownames(cor_matrix), tick = FALSE, font = 2)
        axis(2, at = 1:n_amps, labels = rev(colnames(cor_matrix)), tick = FALSE, font = 2)
        
        for (i in 1:n_amps) {
            for (j in 1:n_amps) {
                val <- cor_matrix[i, j]
                text_col <- if (val > 0.5) "white" else "black"
                text(i, n_amps - j + 1, sprintf("%.2f", val), col = text_col, font = 2, cex = 1.0)
            }
        }
        par(op)
        dev.off()
    }
}

# 5. Plot heatmaps with custom requested titles
plot_heatmap(cor_matrix_200, "outputs/figures_ensemble/amplitude_grid_heatmap_win200_t4200.png", "Correlation between perturbation of different intensities\n(Window = 200 time steps of delay + 200 generations)")
plot_heatmap(cor_matrix_50, "outputs/figures_ensemble/amplitude_grid_heatmap_win50_t4200.png", "Correlation between perturbation of different intensities\n(Window = 200 time steps of delay + 50 generations)")

log_info("All t=4200 heatmaps generated successfully!")
