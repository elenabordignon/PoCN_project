# generate_plots.R
# Script to generate visualizations for the Tangled Nature Model perturbation analysis,
# including the baseline configurations, baselines, and robustness checks.

source("scripts/helpers.R")
library(gplots)
library(pheatmap)
library(grid)

ensure_dir("outputs/figures")

log_info("Loading analysis results...")
load("data/processed/analysis_results.RData")

# Map core species to numbers 1 to 12
species_names <- colnames(X_ref_A)
species_numbers <- as.character(1:length(species_names))
names(species_numbers) <- species_names

# Simplify names in matrices and hclust objects for cleaner plotting
rownames(D_A) <- species_numbers[rownames(D_A)]
colnames(D_A) <- species_numbers[colnames(D_A)]
hclust_A$labels <- species_numbers[hclust_A$labels]

rownames(D_B) <- species_numbers[rownames(D_B)]
colnames(D_B) <- species_numbers[colnames(D_B)]
hclust_B$labels <- species_numbers[hclust_B$labels]

rownames(D_C) <- species_numbers[rownames(D_C)]
colnames(D_C) <- species_numbers[colnames(D_C)]
hclust_C$labels <- species_numbers[hclust_C$labels]

rownames(D_corr) <- species_numbers[rownames(D_corr)]
colnames(D_corr) <- species_numbers[colnames(D_corr)]
hclust_baseline$labels <- species_numbers[hclust_baseline$labels]

# Define premium color palettes
# A beautiful gradient from light ice-blue/white to deep rich navy/midnight-blue
dist_palette <- colorRampPalette(c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594", "#081d58"))(100)
# Correlation palette: red to blue for negative/positive correlations
corr_palette <- colorRampPalette(c("#08306b", "#4292c6", "#f7f7f7", "#ef3b2c", "#67000d"))(100)

# ==============================================================================
# CONFIGURATION A: Amplitude = 1.0, Window = 200 (Baseline)
# ==============================================================================
log_info("Generating Config A: perturbation distance heatmap...")
png("outputs/figures/perturbation_distance_heatmap.png", width = 800, height = 800, res = 120)
pheatmap(D_A, 
         cluster_rows = hclust_A, 
         cluster_cols = hclust_A, 
         color = dist_palette, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "Heatmap showing functional clusters", 
         angle_col = 0)
grid.text("Species ID", x = 0.04, y = 0.44, rot = 90, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()

log_info("Generating Config A: perturbation clustering dendrogram...")
png("outputs/figures/perturbation_dendrogram.png", width = 800, height = 600, res = 120)
op <- par(mar = c(5, 6, 4, 2) + 0.1)
plot(hclust_A, 
     main = "Dendrogram",
     xlab = "Species ID", 
     sub = "", 
     ylab = "", 
     col = "#084594",
     lwd = 2.5,
     hang = -1,
     axes = FALSE)
axis(2, las = 1)
title(ylab = "Integrated Manhattan Distance", line = 4.5)
rect.hclust(hclust_A, k = 3, border = c("#de2d26", "#3182bd", "#31a354"))
par(op)
dev.off()

# ==============================================================================
# CONFIGURATION B: Amplitude = 0.5, Window = 200 (Robustness: Amplitude 0.5)
# ==============================================================================
log_info("Generating Config B: perturbation distance heatmap...")
png("outputs/figures/perturbation_distance_heatmap_amp0.5.png", width = 800, height = 800, res = 120)
pheatmap(D_B, 
         cluster_rows = hclust_B, 
         cluster_cols = hclust_B, 
         color = dist_palette, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "heatmap (Config B)", 
         angle_col = 0)
grid.text("Species ID", x = 0.04, y = 0.44, rot = 90, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()

log_info("Generating Config B: perturbation clustering dendrogram...")
png("outputs/figures/perturbation_dendrogram_amp0.5.png", width = 800, height = 600, res = 120)
op <- par(mar = c(5, 6, 4, 2) + 0.1)
plot(hclust_B, 
     main = "Dendrogram (Config B)",
     xlab = "Species ID", 
     sub = "", 
     ylab = "", 
     col = "#2171b5",
     lwd = 2.5,
     hang = -1,
     axes = FALSE)
axis(2, las = 1)
title(ylab = "Integrated Manhattan Distance", line = 4.5)
rect.hclust(hclust_B, k = 3, border = c("#de2d26", "#3182bd", "#31a354"))
par(op)
dev.off()

# ==============================================================================
# CONFIGURATION C: Amplitude = 1.0, Window = 100 (Robustness: Window 100)
# ==============================================================================
log_info("Generating Config C: perturbation distance heatmap...")
png("outputs/figures/perturbation_distance_heatmap_win100.png", width = 800, height = 800, res = 120)
pheatmap(D_C, 
         cluster_rows = hclust_C, 
         cluster_cols = hclust_C, 
         color = dist_palette, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "heatmap (Config C)", 
         angle_col = 0)
grid.text("Species ID", x = 0.04, y = 0.44, rot = 90, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()

log_info("Generating Config C: perturbation clustering dendrogram...")
png("outputs/figures/perturbation_dendrogram_win100.png", width = 800, height = 600, res = 120)
op <- par(mar = c(5, 6, 4, 2) + 0.1)
plot(hclust_C, 
     main = "Dendrogram (Config C)",
     xlab = "Species ID", 
     sub = "", 
     ylab = "", 
     col = "#4292c6",
     lwd = 2.5,
     hang = -1,
     axes = FALSE)
axis(2, las = 1)
title(ylab = "Integrated Manhattan Distance", line = 4.5)
rect.hclust(hclust_C, k = 3, border = c("#de2d26", "#3182bd", "#31a354"))
par(op)
dev.off()

# ==============================================================================
# ABUNDANCE CORRELATION BASELINE
# ==============================================================================
log_info("Generating abundance correlation distance heatmap...")
C_matrix <- cor(X_ref_A, method = "pearson")
rownames(C_matrix) <- species_numbers[rownames(C_matrix)]
colnames(C_matrix) <- species_numbers[colnames(C_matrix)]

corr_palette_rev <- colorRampPalette(c("#67000d", "#ef3b2c", "#f7f7f7", "#4292c6", "#08306b"))(100)

png("outputs/figures/correlation_baseline_heatmap.png", width = 800, height = 800, res = 120)
pheatmap(C_matrix, 
         cluster_rows = hclust_A, 
         cluster_cols = hclust_A, 
         color = corr_palette_rev, 
         breaks = seq(-1, 1, length.out = 101),
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "Correlation heatmap", 
         angle_col = 0)
grid.text("Species ID", x = 0.04, y = 0.44, rot = 90, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()

log_info("Generating baseline clustering dendrogram...")
png("outputs/figures/correlation_dendrogram.png", width = 800, height = 600, res = 120)
op <- par(mar = c(5, 6, 4, 2) + 0.1)
plot(hclust_baseline, 
     main = "Dendrogram (Correlation Baseline)",
     xlab = "Species ID", 
     sub = "", 
     ylab = "", 
     col = "#4a148c",
     lwd = 2.5,
     hang = -1,
     axes = FALSE)
axis(2, las = 1)
title(ylab = "Correlation Distance (1 - R)", line = 4.5)
rect.hclust(hclust_baseline, k = 3, border = c("#de2d26", "#3182bd", "#31a354"))
par(op)
dev.off()

# ==============================================================================
# PERTURBATION TRAJECTORY PROPAGATION EXAMPLE (Config A)
# ==============================================================================
# Map core species to numbers 1 to 12
species_names <- colnames(X_ref_A)
species_numbers <- as.character(1:length(species_names))
names(species_numbers) <- species_names

log_info("Core Species Mapping:")
for (name in species_names) {
    log_info(paste("Species", species_numbers[name], ":", name))
}

# Pick the species with the highest average abundance to see the propagation
avg_abundance <- colMeans(X_ref_A)
main_species <- names(which.max(avg_abundance))
main_species_num <- species_numbers[main_species]
log_info("Selecting species for trajectory visualization:", main_species, "(Species", main_species_num, ")")

# Get delta_x matrix for this species in Config A
delta_main <- delta_x_list_A[[main_species]]

png("outputs/figures/perturbation_propagation.png", width = 800, height = 500, res = 120)
# Select top 5 species that respond most strongly
max_resp <- colSums(abs(delta_main))
top_responders <- names(sort(max_resp, decreasing = TRUE)[1:5])
top_responders_labels <- paste("Species", species_numbers[top_responders])

# Set up the plot
matplot(time_steps_A, delta_main[, top_responders], 
        type = "l", 
        lty = 1, 
        lwd = 2,
        col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),
        main = paste("Perturbation Propagation over Time (Config A)\n(Perturbing Species:", main_species_num, ")"),
        xlab = "Simulation Generations (t)",
        ylab = expression(paste("Perturbation Vector ", delta, "x(t)")))
legend("topright", legend = top_responders_labels, 
       col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"), 
       lty = 1, lwd = 2, bty = "n")
dev.off()

# Define a consistent color palette for the 12 core species
species_colors <- c(
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
    "#a65628", "#f781bf", "#999999", "#8dd3c7", "#bebada", 
    "#fb8072", "#0173b2"
)
if (length(species_names) > length(species_colors)) {
    species_colors <- rainbow(length(species_names))
}
names(species_colors) <- species_names

# Centered rolling window function to compute mean and standard deviation
compute_rolling_stats <- function(X, window_size = 256) {
    N <- nrow(X)
    M <- ncol(X)
    half_w <- floor(window_size / 2)
    
    mean_mat <- matrix(NA, nrow = N, ncol = M)
    sd_mat <- matrix(NA, nrow = N, ncol = M)
    colnames(mean_mat) <- colnames(X)
    colnames(sd_mat) <- colnames(X)
    
    for (i in 1:N) {
        start_idx <- max(1, i - half_w)
        end_idx <- min(N, i + half_w)
        
        window_data <- X[start_idx:end_idx, , drop = FALSE]
        mean_mat[i, ] <- colMeans(window_data)
        
        if (nrow(window_data) > 1) {
            sd_mat[i, ] <- apply(window_data, 2, sd)
        } else {
            sd_mat[i, ] <- 0
        }
    }
    sd_mat[is.na(sd_mat)] <- 0
    return(list(mean = mean_mat, sd = sd_mat))
}

plot_variance_band <- function(time, mean_vec, sd_vec, color = "grey80", alpha = 0.2) {
    lower <- pmax(0, mean_vec - sd_vec)
    upper <- mean_vec + sd_vec
    
    polygon_x <- c(time, rev(time))
    polygon_y <- c(lower, rev(upper))
    
    polygon(polygon_x, polygon_y, col = adjustcolor(color, alpha.f = alpha), border = NA)
}

# ==============================================================================
# UNPERTURBED TRAJECTORIES
# ==============================================================================
log_info("Generating Config A: unperturbed trajectories with variance bands (window = 256)...")
stats_ref <- compute_rolling_stats(X_ref_A, window_size = 256)

png("outputs/figures/unperturbed_trajectories.png", width = 900, height = 500, res = 120)
# Adjust margins to leave space for the legend on the right
op <- par(mar = c(5, 4, 4, 8) + 0.1)

# Compute global ymax to ensure all variance bands/trajectories fit on the plot
ymax <- max(pmax(X_ref_A, stats_ref$mean + stats_ref$sd))

# Initialize simple plot structure: Time vs Abundance
plot(time_steps_A, X_ref_A[, 1], type = "n",
     xlim = range(time_steps_A),
     ylim = c(0, ymax),
     main = "Unperturbed Abundance Trajectories",
     xlab = "Time",
     ylab = "Abundance")

# 1. Plot the transparent grey variance bands first (so they sit in the background)
for (name in species_names) {
    plot_variance_band(time_steps_A, stats_ref$mean[, name], stats_ref$sd[, name],
                       color = "grey60", alpha = 0.3)
}

# 2. Overlay the solid raw (un-averaged) lines for each species
for (name in species_names) {
    lines(time_steps_A, X_ref_A[, name],
          col = species_colors[name], lwd = 1.5)
}

# Add legend outside the plotting area on the right
legend("topright", inset = c(-0.25, 0), legend = paste("Specie", species_numbers), 
       col = species_colors, 
       lty = 1, lwd = 2, bty = "n", cex = 0.8, xpd = TRUE)
par(op)
dev.off()

# ==============================================================================
# PERPERTURBED TRAJECTORIES
# ==============================================================================
# Reconstruct the perturbed trajectory for Species 11 with robust fallback
target_pert_species <- "N_750181"
target_pert_species_num <- "11"
if (!(target_pert_species %in% names(delta_x_list_A))) {
    idx <- min(11, length(delta_x_list_A))
    target_pert_species <- names(delta_x_list_A)[idx]
    target_pert_species_num <- species_numbers[target_pert_species]
}
log_info("Generating Config A: perturbed trajectories for Species", target_pert_species_num, "(", target_pert_species, ")...")
X_pert_target <- X_ref_A + delta_x_list_A[[target_pert_species]]

png("outputs/figures/perturbed_trajectories.png", width = 900, height = 500, res = 120)
# Adjust margins to leave space for the legend on the right
op <- par(mar = c(5, 4, 4, 8) + 0.1)

# Simply plot Time vs Abundance for the perturbed trajectories
matplot(time_steps_A, X_pert_target, 
        type = "l", 
        lty = 1, 
        lwd = 1.5,
        col = species_colors,
        main = paste("Perturbed Abundance Trajectories (perturbed species :", target_pert_species_num, ")"),
        xlab = "Time",
        ylab = "Abundance")

# Add legend outside the plotting area on the right
legend("topright", inset = c(-0.25, 0), legend = paste("Specie", species_numbers), 
       col = species_colors, 
       lty = 1, lwd = 2, bty = "n", cex = 0.8, xpd = TRUE)
par(op)
dev.off()

log_info("All plots generated successfully in outputs/figures/!")
