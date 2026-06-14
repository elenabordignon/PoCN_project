# data_loader.R
# Functions to read and parse Tangled Nature Model trajectory files

#' Load abundance trajectory from a .dat file
#' @param file_path Character path to the file
#' @return A data frame containing time 't' and populations of core species
load_trajectory <- function(file_path) {
    if (!file.exists(file_path)) {
        stop(paste("File not found:", file_path))
    }
    
    # Read the data, skipping commented lines if any, but keeping column names.
    # Note: R's read.table handles comment.char='#' automatically by default,
    # but the header line starts with '# t N_...' which read.table might drop or fail to parse as header.
    # Let's read the first line separately to get the headers.
    first_line <- readLines(file_path, n = 1)
    
    # Clean the header line (remove leading '#' and split by whitespace)
    header_clean <- gsub("^#\\s*", "", first_line)
    header_tokens <- strsplit(header_clean, "\\s+")[[1]]
    
    # Read the data without the first line
    data_matrix <- read.table(file_path, header = FALSE, comment.char = "#")
    colnames(data_matrix) <- header_tokens
    
    return(data_matrix)
}

#' Extract core species IDs from a trajectory file header
#' @param file_path Character path to the file
#' @return A character vector of core species column names (excluding 't')
extract_core_species <- function(file_path) {
    if (!file.exists(file_path)) {
        stop(paste("File not found:", file_path))
    }
    
    first_line <- readLines(file_path, n = 1)
    header_clean <- gsub("^#\\s*", "", first_line)
    header_tokens <- strsplit(header_clean, "\\s+")[[1]]
    
    # The first token is 't'. The remaining tokens are core species names (e.g. 'N_573792')
    core_species <- header_tokens[header_tokens != "t"]
    return(core_species)
}

#' Load a full dataset (reference and perturbations) from a raw directory
#' @param raw_dir Character path to the directory containing .dat files
#' @return A list with X_ref, pert_list, core_species, and time_steps
load_dataset_from_dir <- function(raw_dir) {
    ref_files <- list.files(path = raw_dir, pattern = "^species_ref_.*\\.dat$", full.names = TRUE)
    if (length(ref_files) == 0) {
        stop(paste("Reference file not found in", raw_dir))
    }
    ref_file <- ref_files[1]
    
    core_species <- extract_core_species(ref_file)
    ref_data <- load_trajectory(ref_file)
    time_steps <- ref_data$t
    X_ref <- as.matrix(ref_data[, core_species])
    
    pert_list <- list()
    for (species_name in core_species) {
        species_id <- gsub("N_", "", species_name)
        pattern <- paste0("^species_pert_", species_id, "_.*\\.dat$")
        pert_files <- list.files(path = raw_dir, pattern = pattern, full.names = TRUE)
        
        if (length(pert_files) > 0) {
            pert_file <- pert_files[1]
            pert_data <- load_trajectory(pert_file)
            X_pert <- as.matrix(pert_data[, core_species])
            pert_list[[species_name]] <- X_pert
        } else {
            warning(paste("No perturbed file found for species", species_id, "in", raw_dir))
        }
    }
    
    return(list(
        X_ref = X_ref,
        pert_list = pert_list,
        core_species = core_species,
        time_steps = time_steps
    ))
}

