# Taxonomic heatmaps for 16S and optional ITS
# Run from the data directory or pass it as the first argument:
# Rscript metagenomics_heatmap.r <data_directory>

required_packages <- c("pheatmap", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
} else getwd()
data_dir <- if (length(args)) args[1] else getwd()

find_file <- function(pattern, label, required = TRUE) {
  dirs <- unique(normalizePath(c(data_dir, script_dir), winslash = "/",
                                mustWork = FALSE))
  found <- unique(unlist(lapply(dirs, function(d) {
    if (dir.exists(d)) list.files(d, pattern = pattern, full.names = TRUE,
                                  ignore.case = TRUE) else character()
  }), use.names = FALSE))
  found <- unique(normalizePath(found, winslash = "/", mustWork = FALSE))
  if (!length(found)) {
    if (required) stop("Could not find ", label, " in data folder or beside script.")
    return(NA_character_)
  }
  if (length(found) > 1) stop("Multiple files found for ", label, ": ",
                              paste(basename(found), collapse = ", "))
  found[1]
}

metadata_file <- find_file("^metadata\\.csv$", "metadata.csv")
tax_files <- c(
  "16S" = find_file("^16S.*_tax_table\\.csv$", "16S tax table"),
  "ITS" = find_file("^ITS.*_tax_table\\.csv$", "ITS tax table", required = FALSE)
)
if (is.na(tax_files[["ITS"]])) {
  tax_files <- tax_files["16S"]
  message("ITS tax table not found; processing 16S only.")
}
data_dir <- dirname(metadata_file)
output_dir <- file.path(data_dir, "metagenomics_heatmap_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_auto <- function(file) {
  readr::read_delim(
    file, delim = NULL, trim_ws = TRUE, name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, na = c("", "NA")
  )
}
parse_numbers <- function(x) {
  x <- gsub(",", ".", trimws(as.character(x)), fixed = TRUE)
  m <- regexpr("[-+]?(?:[0-9]*\\.)?[0-9]+(?:[eE][-+]?[0-9]+)?", x, perl = TRUE)
  ok <- m > 0
  out <- rep(NA_real_, length(x))
  ml <- attr(m, "match.length")
  out[ok] <- suppressWarnings(as.numeric(substr(x[ok], m[ok], m[ok] + ml[ok] - 1L)))
  out
}

metadata <- read_auto(metadata_file)
id_candidates <- c("Sample_Name", "Sample_ID", "SampleID", "Sample")
id_column <- intersect(id_candidates, names(metadata))[1]
if (is.na(id_column)) {
  stop("metadata.csv needs a sample ID column named Sample_Name, Sample_ID, SampleID, or Sample.")
}
metadata$Sample_ID <- trimws(as.character(metadata[[id_column]]))
if (anyNA(metadata$Sample_ID) || any(!nzchar(metadata$Sample_ID)) ||
    anyDuplicated(metadata$Sample_ID)) {
  stop("Metadata sample IDs are blank or duplicated.")
}
metadata_columns <- setdiff(names(metadata), c(id_candidates, "Sample_ID"))
numeric_features <- character()
categorical_features <- character()
numeric_pattern <- paste0("^[-+]?[0-9][0-9.,]*(?:[eE][-+]?[0-9]+)?",
                          "[[:space:]]*[[:alpha:]%°/]*$")
for (v in metadata_columns) {
  x <- trimws(as.character(metadata[[v]]))
  present <- !is.na(x) & nzchar(x)
  if (any(present) && all(grepl(numeric_pattern, x[present], perl = TRUE))) {
    metadata[[v]] <- parse_numbers(x)
    numeric_features <- c(numeric_features, v)
  } else {
    x[!present] <- "None"
    metadata[[v]] <- factor(x)
    categorical_features <- c(categorical_features, v)
  }
}
message("Metadata variables included: ",
        if (length(metadata_columns)) paste(metadata_columns, collapse = ", ") else "none")
message("Numeric metadata variables found: ",
        if (length(numeric_features)) paste(numeric_features, collapse = ", ") else "none")

# Treatment label is the sample ID prefix before its first number.
treatment <- sub("[0-9].*$", "", metadata$Sample_ID)
treatment <- sub("[._ -]+$", "", treatment)
treatment[!nzchar(treatment)] <- NA_character_
has_treatment <- any(!is.na(treatment))
if (has_treatment) {
  metadata$Treatment <- factor(ifelse(is.na(treatment), "Unclassified", treatment))
  message("Treatment prefixes found: ", paste(unique(na.omit(treatment)), collapse = ", "))
}

# Heatmap gradient and categorical annotation.
heatmap_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
heatmap_breaks <- seq(-2.5, 2.5, length.out = 101)

for (marker in names(tax_files)) {
  tax_file <- tax_files[[marker]]
  if (!file.exists(tax_file)) stop("Input file not found: ", tax_file)
  tax <- read_auto(tax_file)
  if (ncol(tax) < 3) stop("Tax table needs a taxon column and at least two sample columns: ", tax_file)
  taxa <- make.unique(trimws(as.character(tax[[1]])))
  sample_ids <- trimws(names(tax)[-1])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Tax-table sample IDs are blank or duplicated: ", tax_file)
  }
  unmatched <- setdiff(sample_ids, metadata$Sample_ID)
  if (length(unmatched)) {
    stop("Metadata is missing sample IDs in ", basename(tax_file), ": ",
         paste(unmatched, collapse = ", "))
  }

  abundance <- do.call(cbind, lapply(tax[-1], parse_numbers))
  rownames(abundance) <- taxa
  colnames(abundance) <- sample_ids
  storage.mode(abundance) <- "numeric"
  abundance[!is.finite(abundance) | abundance < 0] <- 0
  keep_n <- min(100L, nrow(abundance))
  if (!keep_n) stop("No taxa found in: ", tax_file)
  abundance <- abundance[order(rowMeans(abundance, na.rm = TRUE),
                               decreasing = TRUE)[seq_len(keep_n)], , drop = FALSE]
  positive <- abundance[abundance > 0]
  if (!length(positive)) stop("No positive abundances in: ", tax_file)

  log_abundance <- log2(abundance + min(positive) / 2)
  z <- t(scale(t(log_abundance)))
  z[!is.finite(z)] <- 0
  if (ncol(z) < 2) stop("Fewer than two samples found in: ", tax_file)
  message(marker, " heatmap dimensions: ", nrow(z), " taxa x ", ncol(z), " samples")

  md_rows <- match(colnames(z), metadata$Sample_ID)
  annotation <- as.data.frame(
    metadata[md_rows, metadata_columns, drop = FALSE],
    check.names = FALSE
  )
  rownames(annotation) <- colnames(z)
  if (has_treatment) annotation$Treatment <- metadata$Treatment[md_rows]

  annotation_colors <- list()
  for (v in categorical_features) {
    vals <- droplevels(factor(annotation[[v]]))
    annotation[[v]] <- vals
    lev <- levels(vals)
    cols <- if (v == "Fermentation_batch") {
      grDevices::gray.colors(length(lev), start = 0.85, end = 0.25)
    } else grDevices::hcl.colors(length(lev), palette = "Dark 3")
    annotation_colors[[v]] <- stats::setNames(cols, lev)
  }
  if (has_treatment) {
    vals <- droplevels(factor(annotation$Treatment))
    annotation$Treatment <- vals
    annotation_colors$Treatment <- stats::setNames(
      grDevices::hcl.colors(nlevels(vals), palette = "Dark 3"), levels(vals)
    )
  }

  row_tree <- if (nrow(z) > 1) hclust(dist(z), method = "average") else FALSE
  col_tree <- hclust(dist(t(z)), method = "average")
  output_file <- file.path(output_dir, paste0(
    marker, "_taxonomy_heatmap_top100_clustered_metadata.png"
  ))
  grDevices::png(output_file, width = 3600, height = 2600, res = 200)
  tryCatch({
    p <- pheatmap::pheatmap(
      mat = z, color = heatmap_palette, breaks = heatmap_breaks, scale = "none",
      cluster_rows = row_tree, cluster_cols = col_tree,
      annotation_col = if (ncol(annotation)) annotation else NULL,
      annotation_colors = annotation_colors, border_color = NA, angle_col = 45,
      fontsize = 9, fontsize_row = 8, fontsize_col = 8,
      main = paste0(marker, " Taxonomic Heatmap\nTop ", keep_n,
                    " Taxa; Log2 Relative Abundance; Row-wise Z-Scores"),
      silent = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = grDevices::dev.off())
  size <- file.info(output_file)$size
  if (is.na(size) || size == 0) stop("PNG output is empty: ", output_file)
  message("Heatmap saved: ", output_file, " (", size, " bytes)")
}
message("All available heatmaps saved in: ", output_dir)
