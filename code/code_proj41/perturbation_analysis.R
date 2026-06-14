# perturbation_analysis.R
# Core calculations for perturbation-response geometry of the Tangled Nature Model

#' Compute perturbation vectors delta_x for each species
#' @param ref_matrix Numeric matrix of reference abundances (T x M)
#' @param pert_list List of numeric matrices of perturbed abundances, indexed by species name
#' @return List of matrices (T x M) containing delta_x(t) for each perturbed species
compute_delta_x <- function(ref_matrix, pert_list) {
    delta_x_list <- list()
    species_names <- colnames(ref_matrix)
    
    for (name in species_names) {
        if (!is.null(pert_list[[name]])) {
            # delta_x(t) = x_pert(t) - x_ref(t)
            delta_x_list[[name]] <- pert_list[[name]] - ref_matrix
        } else {
            warning(paste("No perturbed trajectory found for species", name))
        }
    }
    
    return(delta_x_list)
}

#' Compute time-integrated perturbation distance matrix D using Manhattan norm
#' @param delta_x_list List of delta_x matrices (T x M) for each species
#' @return A symmetric numeric matrix of size M x M of time-integrated Manhattan distances
compute_integrated_distance_matrix <- function(delta_x_list) {
    species_names <- names(delta_x_list)
    M <- length(species_names)
    
    # Initialize distance matrix
    D <- matrix(0, nrow = M, ncol = M)
    rownames(D) <- species_names
    colnames(D) <- species_names
    
    for (i in 1:M) {
        name_i <- species_names[i]
        delta_i <- delta_x_list[[name_i]]
        
        for (j in i:M) {
            name_j <- species_names[j]
            if (i == j) {
                D[i, j] <- 0
            } else {
                delta_j <- delta_x_list[[name_j]]
                
                # Compute generation-by-generation Manhattan distance (L1 norm):
                # d_t(i,j) = sum_k |delta_i(t,k) - delta_j(t,k)|
                # We can do this efficiently in R using rowSums and abs:
                d_t <- rowSums(abs(delta_i - delta_j))
                
                # Time integral = sum of distances over the time window
                integrated_dist <- sum(d_t)
                
                D[i, j] <- integrated_dist
                D[j, i] <- integrated_dist
            }
        }
    }
    
    return(D)
}

#' Perform hierarchical clustering on a distance matrix
#' @param dist_matrix Numeric matrix of pairwise distances
#' @param method Character string specifying clustering linkage method (default "average")
#' @return An hclust object
cluster_species <- function(dist_matrix, method = "average") {
    # Convert to dist object
    r_dist <- as.dist(dist_matrix)
    hclust_res <- hclust(r_dist, method = method)
    return(hclust_res)
}

#' Compute baseline abundance correlation distance matrix
#' @param ref_matrix Numeric matrix of reference abundances (T x M)
#' @return A symmetric correlation distance matrix: D_corr(i,j) = 1 - cor(ref_i, ref_j)
compute_correlation_baseline <- function(ref_matrix) {
    # Compute Pearson correlation matrix
    C <- cor(ref_matrix, method = "pearson")
    
    # Pearson distance: D_corr(i, j) = 1 - C(i, j)
    # This maps correlation of 1 to 0, and correlation of -1 to 2.
    D_corr <- 1 - C
    
    # Avoid tiny floating-point issues making it slightly negative
    D_corr[D_corr < 0] <- 0
    
    return(D_corr)
}
