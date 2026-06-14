# plot_perturbed_vectors_species11.R
# This script generates a plot of the perturbed vectors delta_x(t) for Species 11 (N_750181)

# Load the processed results
rdata_path <- "data/processed/analysis_results.RData"
if (!file.exists(rdata_path)) {
  stop(paste("File not found:", rdata_path))
}
load(rdata_path)

# Ensure output directories exist
dir.create("output_data", showWarnings = FALSE)
dir.create("tangled_flat", showWarnings = FALSE)

# Species 11
target_pert_species <- "N_750181"
target_pert_species_num <- "11"

# Get delta_x matrix for this species in Config A
delta_target <- delta_x_list_A[[target_pert_species]]
species_names <- colnames(delta_target)

# Define numbers 1 to 12
species_numbers <- as.character(1:length(species_names))
names(species_numbers) <- species_names

# Define colors
species_colors <- c(
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
    "#a65628", "#f781bf", "#999999", "#8dd3c7", "#bebada", 
    "#fb8072", "#0173b2"
)
names(species_colors) <- species_names

# Plotting function
generate_plot <- function(filepath) {
  png(filepath, width = 900, height = 500, res = 120)
  op <- par(mar = c(5, 4, 4, 8) + 0.1)
  
  # Compute y limits dynamically
  ylim <- range(delta_target)
  
  # Plot matrix
  matplot(time_steps_A, delta_target, 
          type = "l", 
          lty = 1, 
          lwd = 1.8,
          col = species_colors,
          ylim = ylim,
          main = "perturbed vector (perturbed specie: 11)",
          xlab = "Time",
          ylab = "Perturbed value")
  
  # Draw a horizontal line at 0 for reference
  abline(h = 0, col = "grey60", lty = 2, lwd = 1.2)
  
  # Add legend outside the plotting area
  legend("topright", inset = c(-0.25, 0), legend = paste("Specie", species_numbers), 
         col = species_colors, 
         lty = 1, lwd = 2, bty = "n", cex = 0.8, xpd = TRUE)
  
  par(op)
  dev.off()
}

# Generate plots in both output_data and tangled_flat
generate_plot("output_data/perturbed_vectors_species11.png")
generate_plot("tangled_flat/perturbed_vectors_species11.png")

message("Plots generated successfully:")
message("- output_data/perturbed_vectors_species11.png")
message("- tangled_flat/perturbed_vectors_species11.png")
