# Metabolomics heatmap with automatically detected metadata annotations.
# Place this script, a metabolite .xlsx file, and metadata.csv in the same
# directory, then run: Rscript metabolomics_heatmap.r

required_packages <- c("readxl", "pheatmap")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(pheatmap))

script_args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(script_args)) script_args[[1]] else getwd()
input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)

find_one <- function(pattern, description) {
  paths <- list.files(input_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(paths)) stop("Could not find ", description, " in ", input_dir)
  if (length(paths) > 1L) stop("Found multiple ", description, " files: ", paste(basename(paths), collapse = ", "), ". Keep one in the input folder.")
  paths[[1]]
}
input_file <- find_one("\\.xlsx?$", "metabolite Excel")
metadata_file <- find_one("^metadata\\.csv$", "metadata.csv")

output_dir <- file.path(input_dir, "metabolomics_heatmap_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "metabolite_heatmap_clustered_metadata.png")

data <- read_excel(input_file, sheet = 1, na = c("-", "NA", ""))
if (ncol(data) < 2L) stop("The Excel sheet must contain a metabolite column and sample columns.")
metabolite_names <- make.unique(trimws(as.character(data[[1]])))
measurements <- as.data.frame(data[-1], check.names = FALSE)
rownames(measurements) <- metabolite_names
measurements[] <- lapply(measurements, function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
})

group_names <- sub("[-._](iii|ii|i)$", "", names(measurements), ignore.case = TRUE)
group_order <- unique(group_names)
if (length(group_order) < 2L) stop("Could not identify at least two sample groups from the Excel column names.")
mean_matrix <- sapply(group_order, function(group) {
  ix <- which(group_names == group)
  vals <- as.matrix(measurements[, ix, drop = FALSE])
  means <- rowMeans(vals, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  means
})
mean_matrix <- as.matrix(mean_matrix)
rownames(mean_matrix) <- metabolite_names
colnames(mean_matrix) <- group_order

log_matrix <- mean_matrix
log_matrix[!is.finite(log_matrix) | log_matrix <= 0] <- NA_real_
log_matrix <- log2(log_matrix)
if (all(is.na(log_matrix))) stop("No positive metabolite measurements are available for log2 transformation.")

z_matrix <- t(apply(log_matrix, 1, function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  s <- sd(x[ok])
  if (!is.finite(s) || s == 0) out[ok] <- 0
  else out[ok] <- (x[ok] - mean(x[ok])) / s
  out
}))
dimnames(z_matrix) <- dimnames(log_matrix)
ok <- is.finite(z_matrix)
z_matrix[ok] <- pmax(-2.5, pmin(2.5, z_matrix[ok]))

lines <- readLines(metadata_file, n = 5L, warn = FALSE, encoding = "UTF-8")
separators <- c(comma = ",", semicolon = ";", tab = "\t")
counts <- vapply(separators, function(sep) {
  vapply(lines, function(line) {
    if (!nzchar(line)) return(0)
    nchar(line) - nchar(gsub(sep, "", line, fixed = TRUE))
  }, numeric(1))
}, numeric(length(lines)))
separator <- names(separators)[which.max(colMeans(counts))]
sep_char <- unname(separators[[separator]])
metadata <- read.table(metadata_file, header = TRUE, sep = sep_char,
  dec = if (sep_char == ";") "," else ".", quote = "\"",
  na.strings = c("", "NA", "-"), stringsAsFactors = FALSE,
  check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(metadata) <- trimws(names(metadata))

id_candidates <- c("Sample_Name", "Sample", "SampleID", "Sample_ID", "Sample name", "Sample ID", "sample")
id_match <- match(tolower(id_candidates), tolower(names(metadata)), nomatch = 0L)
id_match <- id_match[id_match > 0L]
if (length(id_match)) {
  sample_col <- names(metadata)[id_match[[1]]]
} else {
  overlap <- vapply(metadata, function(x) sum(trimws(as.character(x)) %in% group_order, na.rm = TRUE), numeric(1))
  if (!length(overlap) || max(overlap) == 0L) stop("Could not identify a metadata sample ID column matching the Excel sample names.")
  sample_col <- names(metadata)[which.max(overlap)]
}
metadata[[sample_col]] <- trimws(as.character(metadata[[sample_col]]))
if (anyDuplicated(metadata[[sample_col]])) {
  warning("Metadata has duplicate sample IDs; using the first matching row.")
  metadata <- metadata[!duplicated(metadata[[sample_col]]), , drop = FALSE]
}

variable_cols <- setdiff(names(metadata), sample_col)
for (feature in variable_cols) {
  raw <- trimws(as.character(metadata[[feature]]))
  cleaned <- gsub(",", ".", raw, fixed = TRUE)
  if (tolower(feature) %in% c("time", "day", "days", "timepoint")) cleaned <- gsub("[^0-9.+-]", "", cleaned)
  converted <- suppressWarnings(as.numeric(cleaned))
  populated <- !is.na(raw) & nzchar(raw)
  if (any(populated) && all(!is.na(converted[populated]))) {
    metadata[[feature]] <- converted
  } else {
    vals <- raw
    vals[is.na(vals) | !nzchar(vals)] <- "None"
    metadata[[feature]] <- factor(vals)
  }
}
message("Metadata sample ID column: ", sample_col)
message("Metadata variables detected: ", if (length(variable_cols)) paste(variable_cols, collapse = ", ") else "(none)")

metadata_match <- match(colnames(z_matrix), metadata[[sample_col]])
if (anyNA(metadata_match)) {
  # Some exports place trailing sample letters before replicate digits (e.g.
  # FBSPa1 vs FBSP1a). For unmatched IDs, compare letters and digits separately.
  normalize_sample_id <- function(x) {
    x <- toupper(gsub("[^[:alnum:]]", "", trimws(as.character(x))))
    letters <- gsub("[^A-Z]", "", x)
    digits <- gsub("[^0-9]", "", x)
    paste0(letters, digits)
  }
  target_keys <- normalize_sample_id(colnames(z_matrix))
  metadata_keys <- normalize_sample_id(metadata[[sample_col]])
  if (anyDuplicated(metadata_keys)) stop("Sample IDs become ambiguous after ignoring punctuation and letter/digit order.")
  fallback <- match(target_keys, metadata_keys)
  use_fallback <- is.na(metadata_match) & !is.na(fallback)
  metadata_match[use_fallback] <- fallback[use_fallback]
}
if (anyNA(metadata_match)) stop("Metadata is missing sample groups: ", paste(colnames(z_matrix)[is.na(metadata_match)], collapse = ", "))
annotation_col <- as.data.frame(metadata[metadata_match, variable_cols, drop = FALSE], check.names = FALSE)
rownames(annotation_col) <- colnames(z_matrix)

annotation_colors <- list()
for (feature in names(annotation_col)) {
  if (is.factor(annotation_col[[feature]]) || is.character(annotation_col[[feature]])) {
    annotation_col[[feature]] <- factor(annotation_col[[feature]])
    lev <- levels(annotation_col[[feature]])
    palette <- if (tolower(feature) == "fermentation_batch") {
      grDevices::gray.colors(length(lev), start = 0.85, end = 0.25)
    } else {
      grDevices::hcl.colors(length(lev), palette = "Dark 3")
    }
    annotation_colors[[feature]] <- stats::setNames(palette, lev)
  }
}
message("Matched samples: ", ncol(z_matrix), "; metabolites: ", nrow(z_matrix))

cluster_matrix <- z_matrix
for (i in seq_len(nrow(cluster_matrix))) {
  x <- cluster_matrix[i, ]
  ok <- is.finite(x)
  if (any(ok)) x[!ok] <- mean(x[ok]) else x[] <- 0
  cluster_matrix[i, ] <- x
}
row_tree <- hclust(dist(cluster_matrix), method = "average")
column_tree <- hclust(dist(t(cluster_matrix)), method = "average")

png(output_file, width = 3600, height = 3200, res = 200)
tryCatch({
  p <- pheatmap(mat = z_matrix,
    color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101), na_col = "#EEEEEE",
    scale = "none", cluster_rows = row_tree, cluster_cols = column_tree,
    annotation_col = annotation_col, annotation_colors = annotation_colors,
    border_color = NA, angle_col = 45,
    main = "Metabolite Heatmap\nLog2 Data; Row-wise Z-Scores", silent = TRUE)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
}, finally = grDevices::dev.off())
size <- file.info(output_file)$size
if (is.na(size) || size == 0) stop("The output PNG is empty.")
message("Heatmap saved: ", normalizePath(output_file, winslash = "/"))
message("PNG size: ", size, " bytes")




