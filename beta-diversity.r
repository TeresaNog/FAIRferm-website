# FAIRferm beta-diversity analysis for 16S and ITS
# Multi-panel PCoA layout with observed categorical groups and detected
# numeric metadata features.
# Adds PERMANOVA using independent fermentation batches for treatments and
# within-batch restricted permutations for continuous predictors.

# Packages 
required_packages <- c("readr", "vegan", "ggplot2", "gridExtra")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Input discovery. Put metadata.csv and both tax tables in the working directory,
# or pass their directory as the first argument: Rscript Beta-diversity.r <data_dir>.
input_args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
data_dir <- if (length(input_args)) input_args[1] else getwd()

find_one_file <- function(directory, pattern, label) {
  matches <- list.files(directory, pattern = pattern, full.names = TRUE,
                        ignore.case = TRUE)
  if (!length(matches) && !identical(normalizePath(directory, mustWork = FALSE),
                                     normalizePath(script_dir, mustWork = FALSE))) {
    matches <- list.files(script_dir, pattern = pattern, full.names = TRUE,
                          ignore.case = TRUE)
  }
  if (!length(matches)) stop("Could not find ", label, " in ", directory,
                             " or beside the script.")
  if (length(matches) > 1) stop("Found multiple files for ", label, ": ",
                                paste(basename(matches), collapse = ", "))
  matches[1]
}

metadata_file <- find_one_file(data_dir, "^metadata\\.csv$", "metadata.csv")
taxon_files <- c(
  "16S" = find_one_file(data_dir, "^16S.*_tax_table\\.csv$", "16S tax table"),
  "ITS" = find_one_file(data_dir, "^ITS.*_tax_table\\.csv$", "ITS tax table")
)
data_dir <- dirname(metadata_file)
output_dir <- file.path(data_dir, "beta_diversity_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
permutations <- 9999
set.seed(123)

# Helpers
normalize_sample_id <- function(x) {
  trimws(as.character(x))
}

parse_number <- function(x) {
  x <- gsub(",", ".", trimws(as.character(x)), fixed = TRUE)
  match <- regexpr("[-+]?(?:[0-9]*\\.)?[0-9]+(?:[eE][-+]?[0-9]+)?", x)
  found <- match > 0
  result <- rep(NA_real_, length(x))
  match_length <- attr(match, "match.length")
  result[found] <- suppressWarnings(as.numeric(substr(
    x[found], match[found], match[found] + match_length[found] - 1L
  )))
  result
}

treatment_from_sample_id <- function(x) {
  # Treatment is encoded as the sample-name prefix before its first number,
  # e.g. FBSP1a -> FBSP, FBI1a -> FBI, and FBIP1a -> FBIP.
  treatment <- sub("[0-9].*$", "", trimws(as.character(x)))
  treatment <- sub("[._ -]+$", "", treatment)
  if (anyNA(treatment) || any(!nzchar(treatment))) {
    stop("Could not infer a treatment prefix from one or more sample IDs. ",
         "Sample IDs must start with the treatment label followed by a number.")
  }
  treatment
}

read_taxon_table <- function(file, metadata) {
  raw <- readr::read_delim(
    file, delim = NULL, name_repair = "minimal", trim_ws = TRUE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
  if (ncol(raw) < 3) stop("Expected taxon IDs plus sample columns in: ", file)

  taxa <- make.unique(trimws(as.character(raw[[1]])))
  sample_ids <- normalize_sample_id(names(raw)[-1])
  matched <- sample_ids %in% metadata$Sample_ID
  if (sum(matched) < 3) stop("Fewer than three samples matched between table and metadata: ", file)
  if (anyDuplicated(sample_ids[matched])) stop("Duplicate sample IDs after normalization in: ", file)

  # Convert within this function; `dat` is local and is not used outside it.
  dat <- as.data.frame(
    lapply(raw[, which(matched) + 1, drop = FALSE], parse_number),
    check.names = FALSE
  )
  abundance <- as.matrix(dat)
  storage.mode(abundance) <- "numeric"
  rownames(abundance) <- taxa
  colnames(abundance) <- sample_ids[matched]
  abundance[!is.finite(abundance)] <- 0
  abundance <- t(abundance) # samples in rows, taxa in columns

  sample_order <- metadata$Sample_ID[metadata$Sample_ID %in% rownames(abundance)]
  abundance <- abundance[sample_order, , drop = FALSE]
  md <- as.data.frame(
    metadata[match(sample_order, metadata$Sample_ID), , drop = FALSE],
    check.names = FALSE
  )
  rownames(md) <- sample_order
  stopifnot(identical(rownames(abundance), rownames(md)))
  if (any(abundance < 0)) stop("Negative abundance values found in: ", file)
  totals <- rowSums(abundance)
  if (any(totals <= 0)) stop("At least one sample has zero total abundance in: ", file)
  abundance <- sweep(abundance, 1, totals, "/")
  list(abundance = abundance, metadata = md)
}

run_permanova <- function(distance, predictor, data, strata = NULL) {
  f <- stats::as.formula(paste0("distance ~ `", gsub("`", "", predictor), "`"))
  fit <- if (is.null(strata)) {
    vegan::adonis2(f, data = data, permutations = permutations)
  } else {
    vegan::adonis2(f, data = data, permutations = permutations, strata = strata)
  }
  r <- as.data.frame(fit)[1, , drop = FALSE]
  data.frame(
    Predictor = predictor, n = attr(distance, "Size"), Df = r$Df,
    SumOfSqs = r$SumOfSqs, R2 = r$R2, F = r$F,
    p_value = r$`Pr(>F)`, stringsAsFactors = FALSE
  )
}

make_categorical_panel <- function(scores, variable, title, colors, axis_percent) {
  p <- ggplot2::ggplot(scores, ggplot2::aes(PCoA1, PCoA2, color = .data[[variable]])) +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::stat_ellipse(
      data = scores,
      ggplot2::aes(group = .data[[variable]], color = .data[[variable]]),
      type = "norm", level = 0.95, linewidth = 0.7, na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::labs(
      title = title, color = title,
      x = sprintf("PCoA1 (%.1f%%)", axis_percent[1]),
      y = sprintf("PCoA2 (%.1f%%)", axis_percent[2])
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      panel.grid.minor = ggplot2::element_blank()
    )
  p
}

make_continuous_panel <- function(scores, variable, title, gradient, axis_percent) {
  ggplot2::ggplot(
    scores,
    ggplot2::aes(PCoA1, PCoA2, color = .data[[variable]])
  ) +
    ggplot2::geom_point(size = 2.6, alpha = 0.95, na.rm = TRUE) +
    ggplot2::scale_color_gradientn(
      colours = gradient, name = title, na.value = "#B8B8B8"
    ) +
    ggplot2::labs(
      title = title,
      x = sprintf("PCoA1 (%.1f%%)", axis_percent[1]),
      y = sprintf("PCoA2 (%.1f%%)", axis_percent[2])
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Metadata; readr guesses the separator independently for each input file.
metadata <- readr::read_delim(
  metadata_file, delim = NULL, trim_ws = TRUE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)
id_candidates <- c("Sample_Name", "Sample_ID", "SampleID", "Sample")
id_column <- intersect(id_candidates, names(metadata))[1]
if (is.na(id_column)) stop("Metadata must contain a sample ID column (Sample_Name, Sample_ID, SampleID, or Sample).")
required <- c("Fermentation_batch")
missing <- setdiff(required, names(metadata))
if (length(missing)) stop("Missing required metadata columns: ", paste(missing, collapse = ", "))

metadata$Sample_ID <- normalize_sample_id(metadata[[id_column]])
metadata$Fermentation_batch <- trimws(metadata$Fermentation_batch)
if (anyNA(metadata$Sample_ID) || any(!nzchar(metadata$Sample_ID))) {
  stop("Metadata contains missing or blank sample IDs.")
}
if (anyNA(metadata$Fermentation_batch) || any(!nzchar(metadata$Fermentation_batch))) {
  stop("Metadata contains missing or blank Fermentation_batch values.")
}
culture_columns <- intersect(c("Starter_culture", "Probiotic_culture"), names(metadata))
for (v in culture_columns) {
  metadata[[v]] <- trimws(metadata[[v]])
  metadata[[v]][is.na(metadata[[v]]) | metadata[[v]] %in% c("", "-", "None")] <- NA_character_
}
metadata$Fermentation_treatment <- treatment_from_sample_id(metadata$Sample_ID)
if ("Starter_culture" %in% culture_columns) {
  metadata$Starter_culture[is.na(metadata$Starter_culture)] <- "No starter culture"
}
if ("Probiotic_culture" %in% culture_columns) {
  metadata$Probiotic_culture[is.na(metadata$Probiotic_culture)] <- "No probiotic culture"
}
if (anyDuplicated(metadata$Sample_ID)) stop("Duplicate sample names in metadata after normalization.")

# Detect numeric metadata features from the values present, including values
# such as "9 days" or decimal commas. Identifier and grouping fields are skipped.
excluded_metadata <- c(id_candidates, "Sample_ID", "Fermentation_batch",
                       culture_columns,
                       "Fermentation_treatment")
numeric_features <- character()
for (column in setdiff(names(metadata), excluded_metadata)) {
  raw_values <- trimws(as.character(metadata[[column]]))
  present <- !is.na(raw_values) & nzchar(raw_values)
  if (!any(present)) next
  parsed <- parse_number(raw_values)
  if (all(!is.na(parsed[present]))) {
    metadata[[column]] <- parsed
    numeric_features <- c(numeric_features, column)
  }
}
if (length(numeric_features)) {
  message("Numeric metadata features found: ", paste(numeric_features, collapse = ", "))
} else {
  message("No numeric metadata features found; numeric-feature tests and panels will be omitted.")
}
treatment_groups <- unique(metadata$Fermentation_treatment)
message("Fermentation treatment groups found: ", paste(treatment_groups, collapse = ", "))

# Generate colors from the treatment and culture levels present in this dataset.
make_colors <- function(x) setNames(grDevices::hcl.colors(length(unique(x)), "Dark 3"), unique(x))
treatment_colors <- make_colors(metadata$Fermentation_treatment)
starter_colors <- if ("Starter_culture" %in% culture_columns) make_colors(metadata$Starter_culture) else NULL
probiotic_colors <- if ("Probiotic_culture" %in% culture_columns) make_colors(metadata$Probiotic_culture) else NULL
if (!is.null(starter_colors)) {
  starter_colors["No starter culture"] <- "#BDBDBD"
  starter_colors["Lactiplantibacillus plantarum"] <- "#9467BD"
}
if (!is.null(probiotic_colors)) {
  probiotic_colors["No probiotic culture"] <- "#BDBDBD"
  probiotic_colors["Lacticaseibacillus rhamnosus GG"] <- "#E15759"
}
gradient <- c("#2C3E8F", "#2878B5", "#33A884", "#F0D84A", "#C93832")

# 16S and ITS analyses
all_permanova <- list()
all_permdisp <- list()
all_coordinates <- list()
all_eigenvalues <- list()
all_correlations <- list()

for (marker in names(taxon_files)) {
  file <- taxon_files[[marker]]
  if (!file.exists(file)) stop("Input file not found: ", file)
  input <- read_taxon_table(file, metadata)
  abundance <- input$abundance
  md <- input$metadata

  # Bray-Curtis distances and PCoA are calculated from all sample profiles.
  bray_samples <- vegan::vegdist(abundance, method = "bray")
  ord <- stats::cmdscale(bray_samples, k = 2, eig = TRUE)
  positive_eigenvalues <- ord$eig[ord$eig > 0]
  axis_percent <- 100 * ord$eig[1:2] / sum(positive_eigenvalues)
  scores <- data.frame(
    PCoA1 = ord$points[, 1], PCoA2 = ord$points[, 2],
    Sample_ID = rownames(ord$points), stringsAsFactors = FALSE
  )
  metadata_for_scores <- md[
    match(scores$Sample_ID, rownames(md)),
    ,
    drop = FALSE
  ]
  # Sample_ID is already present in `scores`; do not duplicate that column.
  metadata_for_scores$Sample_ID <- NULL
  scores <- cbind(scores, metadata_for_scores)
  if (anyDuplicated(names(scores))) {
    stop("Duplicate column names in PCoA plotting data for ", marker, ".")
  }
  scores$Fermentation_treatment <- factor(
    scores$Fermentation_treatment, levels = treatment_groups
  )

# Categorical panels: treatment, starter culture, and probiotic culture.
  categorical_panels <- list(
    make_categorical_panel(scores, "Fermentation_treatment", "Treatment", treatment_colors, axis_percent)
  )
  if ("Starter_culture" %in% names(scores)) {
    categorical_panels[[length(categorical_panels) + 1]] <-
      make_categorical_panel(scores, "Starter_culture", "Starter culture", starter_colors, axis_percent)
  }
  if ("Probiotic_culture" %in% names(scores)) {
    categorical_panels[[length(categorical_panels) + 1]] <-
      make_categorical_panel(scores, "Probiotic_culture", "Probiotic culture", probiotic_colors, axis_percent)
  }
  continuous_panels <- lapply(numeric_features, function(v) {
    make_continuous_panel(scores, v, v, gradient, axis_percent)
  })

  top <- gridExtra::arrangeGrob(grobs = categorical_panels, ncol = 3)
  figure_title <- grid::textGrob(
      paste0(marker, " beta diversity | Bray-Curtis PCoA | ",
             ncol(abundance), " taxa | n = ", nrow(abundance)),
      gp = grid::gpar(fontface = "bold", fontsize = 15)
  )
  figure_footer <- grid::textGrob(
      sprintf("PCoA1 %.1f%% | PCoA2 %.1f%%", axis_percent[1], axis_percent[2]),
      gp = grid::gpar(fontsize = 10)
  )
  if (length(continuous_panels)) {
    bottom <- gridExtra::arrangeGrob(grobs = continuous_panels, ncol = min(4, length(continuous_panels)))
    figure <- gridExtra::arrangeGrob(top, bottom, ncol = 1,
                                     heights = c(1, 1.05), top = figure_title,
                                     bottom = figure_footer)
  } else {
    figure <- gridExtra::arrangeGrob(top, ncol = 1, top = figure_title,
                                     bottom = figure_footer)
  }
  ggplot2::ggsave(
    file.path(output_dir, paste0(marker, "_beta_diversity_PCoA.png")),
    figure, width = 16, height = 10, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(marker, "_beta_diversity_PCoA.pdf")),
    figure, width = 16, height = 10, bg = "white"
  )

  # PERMANOVA categorical tests use independent fermentation batches,
  # averaging repeated time points within each batch.
  batch_ids <- unique(md$Fermentation_batch)
  batch_profiles <- do.call(rbind, lapply(batch_ids, function(batch_id) {
    colMeans(abundance[md$Fermentation_batch == batch_id, , drop = FALSE])
  }))
  rownames(batch_profiles) <- batch_ids
  batch_profiles <- batch_profiles / rowSums(batch_profiles)
  batch_md <- md[match(batch_ids, md$Fermentation_batch), , drop = FALSE]
  rownames(batch_md) <- batch_ids
  batch_md$Fermentation_treatment <- factor(batch_md$Fermentation_treatment,
                                             levels = treatment_groups)
  bray_batches <- vegan::vegdist(batch_profiles, method = "bray")

  cat_results <- do.call(rbind, lapply(
    c("Fermentation_treatment", culture_columns),
    function(v) run_permanova(bray_batches, v, batch_md)
  ))
  cat_results$Analysis_level <- "Fermentation-batch means"
  cat_results$p_Holm <- p.adjust(cat_results$p_value, method = "holm")

  # Continuous predictors use all samples, with permutations blocked by batch.
  bray_mat <- as.matrix(bray_samples)
  continuous_results <- if (length(numeric_features)) do.call(rbind, lapply(
    numeric_features,
    function(v) {
      keep <- complete.cases(md[[v]])
      ids <- rownames(md)[keep]
      d <- stats::as.dist(bray_mat[ids, ids, drop = FALSE])
      md_sub <- md[keep, , drop = FALSE]
      run_permanova(d, v, md_sub, strata = md_sub$Fermentation_batch)
    }
  )) else data.frame()
  if (nrow(continuous_results)) {
    continuous_results$Analysis_level <- "Samples; permutations restricted within fermentation batch"
    continuous_results$p_Holm <- p.adjust(continuous_results$p_value, method = "holm")
  }

  perm <- if (nrow(continuous_results)) {
    rbind(cat_results, continuous_results)
  } else {
    cat_results
  }
  perm$Marker <- marker
  all_permanova[[marker]] <- perm

  dispersion <- vegan::betadisper(bray_batches, batch_md$Fermentation_treatment)
  disp_test <- vegan::permutest(dispersion, permutations = permutations)
  all_permdisp[[marker]] <- data.frame(
    Marker = marker,
    F = unname(disp_test$tab[1, "F"]),
    p_value = unname(disp_test$tab[1, "Pr(>F)"]),
    n_fermentation_batches = nrow(batch_md)
  )

  scores$PCoA1_percent <- axis_percent[1]
  scores$PCoA2_percent <- axis_percent[2]
  scores$Marker <- marker
  all_coordinates[[marker]] <- scores
  if (length(numeric_features)) {
    all_correlations[[marker]] <- do.call(rbind, lapply(treatment_groups, function(group) {
      group_scores <- scores[scores$Fermentation_treatment == group, , drop = FALSE]
      do.call(rbind, lapply(numeric_features, function(feature) {
        do.call(rbind, lapply(c("PCoA1", "PCoA2"), function(axis) {
          valid <- complete.cases(group_scores[[feature]], group_scores[[axis]])
          x <- group_scores[[feature]][valid]
          y <- group_scores[[axis]][valid]
          rho <- p_value <- NA_real_
          if (length(x) >= 3 && length(unique(x)) > 1 && length(unique(y)) > 1) {
            test <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
            rho <- unname(test$estimate)
            p_value <- test$p.value
          }
          data.frame(Marker = marker, Fermentation_treatment = group,
                     Metadata_feature = feature, PCoA_axis = axis, n = length(x),
                     Spearman_rho = rho, p_value = p_value)
        }))
      }))
    }))
  }
  all_eigenvalues[[marker]] <- data.frame(
    Marker = marker,
    Axis = seq_along(ord$eig),
    Eigenvalue = ord$eig,
    Percent_positive_eigenvalue_sum = ifelse(
      ord$eig > 0, 100 * ord$eig / sum(positive_eigenvalues), 0
    )
  )

  write.csv(
    as.matrix(bray_samples),
    file.path(output_dir, paste0(marker, "_bray_curtis_distances.csv")),
    row.names = TRUE
  )
}

# Save statistical results and coordinates 
permanova_table <- do.call(rbind, all_permanova)
permdisp_table <- do.call(rbind, all_permdisp)
coordinates_table <- do.call(rbind, all_coordinates)
eigenvalue_table <- do.call(rbind, all_eigenvalues)

write.csv(permanova_table, file.path(output_dir, "beta_diversity_PERMANOVA.csv"), row.names = FALSE, na = "")
write.csv(permdisp_table, file.path(output_dir, "beta_diversity_PERMDISP.csv"), row.names = FALSE, na = "")
write.csv(coordinates_table, file.path(output_dir, "beta_diversity_PCoA_coordinates.csv"), row.names = FALSE, na = "")
write.csv(eigenvalue_table, file.path(output_dir, "beta_diversity_PCoA_eigenvalues.csv"), row.names = FALSE, na = "")
if (length(all_correlations)) {
  correlations <- do.call(rbind, all_correlations)
  correction_group <- interaction(correlations$Marker, correlations$PCoA_axis,
                                  correlations$Metadata_feature, drop = TRUE)
  correlations$p_Holm_within_marker_axis_feature <- ave(
    correlations$p_value, correction_group,
    FUN = function(p) p.adjust(p, method = "holm")
  )
  write.csv(correlations,
            file.path(output_dir, "beta_diversity_correlations_by_treatment.csv"),
            row.names = FALSE, na = "")
}

cat("Figures and results saved to:\n", output_dir, "\n\n")
print(permanova_table)
print(permdisp_table)
