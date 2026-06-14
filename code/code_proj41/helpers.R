# helpers.R
# General utilities and helpers for the project

#' Save a matrix to a CSV file with correct headers and row names
#' @param mat Numeric matrix to save
#' @param file_path Character path to save location
save_matrix_csv <- function(mat, file_path) {
    write.csv(mat, file = file_path, row.names = TRUE, quote = FALSE)
}

#' Log a formatted message to console
#' @param ... arguments passed to cat
log_info <- function(...) {
    cat("[INFO]", paste(...), "\n")
}

#' Create directory if it does not exist
#' @param dir_path Path to directory
ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
}
