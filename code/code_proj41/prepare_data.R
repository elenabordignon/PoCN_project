# prepare_data.R
# Copy raw output files from Tangled-Nature-master/output_results to data/raw/

src_dir <- "simulator/output_results"
dest_dir <- "data/raw"

if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
}

# Find all .dat files in the source directory
files_to_copy <- list.files(path = src_dir, pattern = "\\.dat$", full.names = TRUE)

if (length(files_to_copy) == 0) {
    stop("No .dat files found in Tangled-Nature-master/output_results! Please check the source directory.")
}

cat("Copying", length(files_to_copy), "files from", src_dir, "to", dest_dir, "...\n")

for (f in files_to_copy) {
    file.copy(from = f, to = dest_dir, overwrite = TRUE)
}

cat("Data preparation complete!\n")
