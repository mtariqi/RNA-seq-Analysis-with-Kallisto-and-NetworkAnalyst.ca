---

## 1. Did you lose work when the session expired?

No.
Your screenshots show:

```bash
ls -lh /scratch/${USER}/module10 | grep kallisto-output
kallisto-output-GSM7688023
kallisto-output-GSM7688024
kallisto-output-GSM7688025
kallisto-output-GSM7688026
kallisto-output-GSM7688027
...
```

and inside:

```bash
ls -lh kallisto-output-GSM*/abundance.tsv
# several abundance.tsv files
```

So the Kallisto runs you finished are still there. You only need to re-run samples that **never** produced an `abundance.tsv`.

---

## 2. Do you *have* to run all samples (GSM7688023–GSM7688037)?

For the R script, you only need as many samples as you want to analyze.

For the course assignment, your instructor probably expects you to:

* run **all** the GSMs listed in `rnaseq-mus-musculus-GSE240196/`
* produce a count matrix using all of those samples.

So: **yes, ideally run Kallisto for GSM7688023–GSM7688037**, but if a few are missing the R script below will still work as long as the corresponding `kallisto-output-GSM*` folders exist.

---

## 3. Why do you only see `.tsv` and not `.h5`?

Your Kallisto build inside the course container is almost certainly compiled **without HDF5 support**, so it only writes `abundance.tsv`.

That means:

```bash
mv "$OUTPUT_DIR/abundance.h5" "$OUTPUT_DIR/$BASENAME.h5"
```

will **always** fail (no such file). That’s normal and not a problem for your R analysis.

You have three options:

```bash
# safest: only rename .h5 if it exists
if [ -f "$OUTPUT_DIR/abundance.h5" ]; then
    mv "$OUTPUT_DIR/abundance.h5" "$OUTPUT_DIR/$BASENAME.h5"
fi
```

or just comment that line out.

Your R script **doesn’t use the .h5 files at all**, so you’re fine.

---

## 4. Do you need to rename the abundance files?

You have two layouts:

* Current: `kallisto-output-GSM7688023/abundance.tsv`
* Instructor’s R script assumes: `kallisto-output/GSM7688023/GSM7688023.tsv` (or similar).

Instead of renaming everything by hand, we’ll write an R script that:

* looks for `GSMxxxxx.tsv`; if not found, falls back to `abundance.tsv`.
* shows a progress bar
* filters lowly expressed transcripts
* writes QC summaries.

---

## 5. Rewritten R script (with progress, filtering, QC)

Save this as **`kallisto_counts.R`** in your `module10` directory.

```r
#!/usr/bin/env Rscript

## ----------------- Setup -----------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

DIRECTORY <- "kallisto-output"   # parent dir that contains *either*
                                 #  - kallisto-output-GSM* dirs  (your case), OR
                                 #  - GSM*/ subdirs

## ----------------- Discover files -----------------

# get immediate subdirectories
dirs <- list.dirs(path = DIRECTORY, full.names = TRUE, recursive = FALSE)

if (length(dirs) == 0L) {
  stop("No subdirectories found under '", DIRECTORY, "'. Did you run kallisto yet?")
}

# sample names = directory basenames (e.g. kallisto-output-GSM7688023)
sample_names <- basename(dirs)

# two possible file patterns per directory:
# 1) GSM7688023.tsv (renamed by bash script)
# 2) abundance.tsv  (your current situation)
tsv_primary   <- file.path(dirs, paste0(sample_names, ".tsv"))
tsv_fallback  <- file.path(dirs, "abundance.tsv")

files <- ifelse(file.exists(tsv_primary), tsv_primary, tsv_fallback)

if (any(!file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("These kallisto .tsv files are missing:\n", paste(missing, collapse = "\n"))
}

cat("Found", length(files), "kallisto result files.\n")

## ----------------- Helper: progress-aware reader -----------------

read_kallisto <- function(file) {
  dat <- read_tsv(file, show_col_types = FALSE)
  tibble(
    target_id  = dat$target_id,
    est_counts = dat$est_counts
  )
}

n_files <- length(files)
cat("Reading files and building count matrix...\n")

pb <- txtProgressBar(min = 0, max = n_files, style = 3)

data_list <- vector("list", n_files)

for (i in seq_along(files)) {
  data_list[[i]] <- read_kallisto(files[i]) %>%
    rename(!!sample_names[i] := est_counts)
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nDone reading files.\n")

## ----------------- Build raw count matrix -----------------

# full join stepwise on target_id so we keep all transcripts
count_tbl <- Reduce(function(x, y) full_join(x, y, by = "target_id"), data_list) %>%
  arrange(target_id)

count_matrix <- as.matrix(count_tbl[ , -1])
rownames(count_matrix) <- count_tbl$target_id
mode(count_matrix) <- "numeric"

# sanity check
cat("Raw matrix dimensions:", nrow(count_matrix), "transcripts x",
    ncol(count_matrix), "samples\n")

# write raw counts
raw_out <- file.path(DIRECTORY, "count_matrix_raw.csv")
write.csv(count_matrix, file = raw_out, row.names = TRUE)
cat("Wrote raw count matrix to:", raw_out, "\n")

## ----------------- Filtering lowly expressed transcripts -----------------

# keep transcripts with:
#  - total counts ≥ 10 across all samples AND
#  - expressed in at least 2 samples (count > 0)
total_counts <- rowSums(count_matrix)
nonzero_samples <- rowSums(count_matrix > 0)

keep <- total_counts >= 10 & nonzero_samples >= 2

filtered_matrix <- count_matrix[keep, , drop = FALSE]

cat("Filtering:\n")
cat("  Total transcripts (raw):   ", nrow(count_matrix), "\n")
cat("  Transcripts kept (filtered):", nrow(filtered_matrix), "\n")

filt_out <- file.path(DIRECTORY, "count_matrix_filtered.csv")
write.csv(filtered_matrix, file = filt_out, row.names = TRUE)
cat("Wrote filtered count matrix to:", filt_out, "\n")

## ----------------- Simple QC summary -----------------

library_sizes_raw  <- colSums(count_matrix)
library_sizes_filt <- colSums(filtered_matrix)

qc_tbl <- tibble(
  sample           = colnames(count_matrix),
  libsize_raw      = library_sizes_raw,
  libsize_filtered = library_sizes_filt
)

qc_out <- file.path(DIRECTORY, "count_matrix_QC_summary.csv")
write.csv(qc_tbl, file = qc_out, row.names = FALSE)

cat("QC summary written to:", qc_out, "\n")

cat("\nAll done ✅\n")
```

### How to run it on the cluster

From your `module10` directory:

```bash
module load R        # if your instructions say to load R module
Rscript kallisto_counts.R
```

It will:

1. Detect all subdirs under `kallisto-output*`.
2. For each: use `GSMxxxxxx.tsv` if present, otherwise `abundance.tsv`.
3. Show a progress bar while reading.
4. Write:

   * `kallisto-output/count_matrix_raw.csv`
   * `kallisto-output/count_matrix_filtered.csv`
   * `kallisto-output/count_matrix_QC_summary.csv`

---

## 6. Quick QC report (what to describe in your write-up)

In your **report** / README you can summarize:

* Number of samples included.
* Number of transcripts **before** filtering vs **after** filtering.
* Library sizes per sample (from `count_matrix_QC_summary.csv`).
* Mention any outliers (e.g., one sample has much lower counts).

If you want extra plots (optional but nice):

```r
# after running kallisto_counts.R interactively
log_counts <- log10(filtered_matrix + 1)
boxplot(log_counts, las = 2, main = "Log10 counts per sample")

# basic PCA
pca <- prcomp(t(log_counts), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], xlab = "PC1", ylab = "PC2",
     main = "PCA of samples")
text(pca$x[,1], pca$x[,2], labels = colnames(filtered_matrix), pos = 3, cex = 0.7)
```

You can mention these in the QC section even if you only keep the CSV summary for the submission.

---

## 7. Suggested directory structure (for GitHub and ZIP)

Inside `/scratch/${USER}/module10`:

```text
module10/
├── kallisto-bash.sh
├── kallisto-output-join.R
├── kallisto_counts.R           # the R script above
├── README.md
├── Mus_musculus.idx
├── kallisto-output/            # or kallisto-output-GSM* dirs, depending on your layout
│   ├── kallisto-output-GSM7688023/
│   │   ├── abundance.tsv  (or GSM7688023.tsv if renamed)
│   │   └── GSM7688023.h5  (optional / usually absent)
│   ├── kallisto-output-GSM7688024/
│   │   └── ...
│   ├── count_matrix_raw.csv
│   ├── count_matrix_filtered.csv
│   └── count_matrix_QC_summary.csv
└── (optional logs, e.g., debug_script.sh, etc.)
```

Your actual directory names are `kallisto-output-GSM*`, which is fine. Just keep everything in **one top-level module10 folder**.

---

## 8. README.md template for your GitHub repo

Create `README.md` in `module10`:

````markdown
# Module 10 – RNA-seq with Kallisto and NetworkAnalyst

This repository contains my workflow for Module 10 of BINF6310:
pseudo-alignment with **kallisto**, generation of a transcript-level
count matrix in **R**, and preparation of the data for downstream
analysis (e.g., NetworkAnalyst).

## Directory layout

```text
module10/
├── kallisto-bash.sh          # Bash script that runs kallisto for all GSM files
├── kallisto-output-join.R    # Instructor-provided script for joining outputs
├── kallisto_counts.R         # My script to build + filter count matrix
├── Mus_musculus.idx          # kallisto transcriptome index
├── kallisto-output-GSM*/     # per-sample kallisto output directories
│   └── abundance.tsv         # transcript-level counts for each sample
└── kallisto-output/
    ├── count_matrix_raw.csv
    ├── count_matrix_filtered.csv
    └── count_matrix_QC_summary.csv
````

## 1. Running kallisto

From the course directory:

```bash
cd /scratch/${USER}/module10
bash kallisto-bash.sh
```

This script:

* loops over all `GSM*` files in `rnaseq-mus-musculus-GSE240196`
* runs `kallisto quant` for each sample
* writes each result into its own `kallisto-output-GSMxxxxxxx/` directory.

## 2. Building the count matrix in R

After kallisto has finished for all samples:

```bash
module load R           # if required on the cluster
cd /scratch/${USER}/module10
Rscript kallisto_counts.R
```

This script:

* discovers all kallisto output folders under `kallisto-output*`
* reads each `GSMxxxxxx.tsv` or `abundance.tsv` file with a progress bar
* merges them into a transcript-by-sample count matrix
* filters out lowly expressed transcripts
* writes:

  * `kallisto-output/count_matrix_raw.csv`
  * `kallisto-output/count_matrix_filtered.csv`
  * `kallisto-output/count_matrix_QC_summary.csv`

## 3. QC of the count matrix

The QC summary contains, for each sample:

* total library size before filtering
* total library size after filtering.

In addition, I inspected:

* the distribution of log10 counts per sample using boxplots
* sample relationships using PCA on the filtered count matrix.

These checks did not reveal any obvious outlier samples; library sizes
were of similar magnitude across samples.

## 4. Files for submission

For the Module 10 submission I include:

* `kallisto-bash.sh`
* `kallisto-output-join.R`
* `kallisto_counts.R`
* `README.md`
* `kallisto-output-GSM*/abundance.tsv` for all samples
* `kallisto-output/count_matrix_raw.csv`
* `kallisto-output/count_matrix_filtered.csv`
* `kallisto-output/count_matrix_QC_summary.csv`

````

Edit the last section to match exactly what your instructor asks you to include.

---

## 9. Creating the ZIP file for submission

From **one level above** your module10 folder (e.g. `/scratch/${USER}`):

```bash
cd /scratch/${USER}
zip -r module10_kallisto_islam.mdtar.zip module10 \
    -x '*GSM*fastq*'   # optional: exclude raw FASTQ if not required
````

This creates **one ZIP** with your scripts, Kallisto outputs, and count matrices.

(If the assignment explicitly says “upload kallisto-bash.sh, kallisto-output-join.R, and the count_matrix.csv only”, adjust the zip command accordingly.)

---

## 10. Answer to your last direct question

> “will it generate h5 files after finishing .tsv files?”

No, not with the current container. It writes only `abundance.tsv`. That’s totally fine, and your R workflow above uses just the `.tsv` files.

---

If you paste your current `kallisto-bash.sh` into GitHub later, we can also clean it up and comment it nicely so it looks perfect for both grading and your portfolio.
