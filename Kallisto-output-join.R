library(readr)
library(dplyr)
library(progress)

cat("Scanning for kallisto output directories...\n")

# 1. Detect Kallisto folder

dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
dirs <- dirs[grepl("kallisto-output-GSM", dirs)]
if (length(dirs) == 0) {
  stop("No kallisto output directory found.")
}

cat(paste("Found", length(dirs), "samples.\n"))


# 2.  Build file paths
files <- file.path(dirs, "abundance.tsv")

# Removing non-existing files
missing <- !file.exists(files)
if (any(missing)) {
  cat("Warning: Missing files detected:\n")
  print(files[missing])
  cat("These samples will be skipped.\n")
  dirs <- dirs[!missing]
  files <- files[!missing]
}


# 3. Function to read in kallisto abundance.tsv and extract counts
read_kallisto <- function(file) {
  x <- read_tsv(file, show_col_types = FALSE)
  y <- x$est_counts
  names(y) <- x$target_id
  return(y)
}

# 4. Progress bar

if (length(files) > 0) {
    pb <- progress_bar$new(  total = length(files),
                             format="Reading [:bar] :percent :eta")
    # Read all files
    all_counts <- list()
    for (i in seq_along(files)) {
        all_counts[[i]] <- read_kallisto(files[i])
        pb$tick()
    }
    
    # 5. Build count matrix

    sample_names <- basename(dirs)

    count_matrix <- do.call(cbind, all_counts)
    colnames(count_matrix) <- sample_names

    cat("\n count matrix dimensions:",
    nrow(count_matrix), "genes x", ncol(count_matrix), "sampls\n")

    # 6. Save output

    write.csv(count_matrix, "count_matrix.csv", row.names = TRUE)

    cat("SUCCESS: count_matrix.csv created!\n")
} else {
    cat("No valid files to process.\n")
}
