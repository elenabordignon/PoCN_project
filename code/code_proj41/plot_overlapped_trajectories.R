# plot_overlapped_trajectories.R
# This script overlays the unperturbed reference trajectories (in shades of blue)
# and the perturbed trajectories (in shades of red) for Species 11 (N_750181)

# Load the processed results
rdata_path <- "data/processed/analysis_results.RData"
if (!file.exists(rdata_path)) {
  stop(paste("File not found:", rdata_path))
}
load(rdata_path)

# Ensure output directory exists
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

target_sp <- "N_750181"
target_sp_num <- "11"

# Reference trajectories (unperturbed)
X_ref <- X_ref_A
# Perturbed trajectories
X_pert <- X_ref + delta_x_list_A[[target_sp]]

# Define color palettes
# 12 shades of blue for unperturbed reference trajectories
blues <- colorRampPalette(c("#c6dbef", "#6baed6", "#084594"))(12)
# 12 shades of red for perturbed trajectories
reds <- colorRampPalette(c("#fcae91", "#fb6a4a", "#a50f15"))(12)

# Compute y limits to fit everything
ylim_range <- range(c(X_ref, X_pert))

png("outputs/figures/overlapped_trajectories.png", width = 900, height = 500, res = 120)
op <- par(mar = c(5, 4, 4, 10) + 0.1)

# Initialize plot
plot(time_steps_A, X_ref[, 1], type = "n", 
     xlim = range(time_steps_A), 
     ylim = ylim_range,
     main = "Reference vs Perturbed Abundances (Perturbing Species: 11)",
     xlab = "Time (Generations)",
     ylab = "Abundance")

# Plot unperturbed trajectories (blues)
for (i in 1:ncol(X_ref)) {
  lines(time_steps_A, X_ref[, i], col = blues[i], lwd = 1.3)
}

# Plot perturbed trajectories (reds)
for (i in 1:ncol(X_pert)) {
  lines(time_steps_A, X_pert[, i], col = reds[i], lwd = 1.5, lty = 2)
}

# Draw vertical line at perturbation time
abline(v = 4000, col = "grey40", lty = 3, lwd = 1.5)

# Add legend outside the plotting area
legend("topright", inset = c(-0.32, 0), 
       legend = c(paste("Unperturbed Sp.", 1:12), paste("Perturbed Sp.", 1:12)), 
       col = c(blues, reds), 
       lty = c(rep(1, 12), rep(2, 12)), 
       lwd = 2, bty = "n", cex = 0.7, xpd = TRUE, ncol = 2)

par(op)
dev.off()

message("Overlapped trajectories plot generated successfully: outputs/figures/overlapped_trajectories.png")
