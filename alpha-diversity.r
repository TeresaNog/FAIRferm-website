# Install any required packages that are not already installed
packages <- c("ggplot2", "ggpubr", "patchwork", "readr", "dplyr", "tidyr", "stringr")

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(ggpubr)
library(patchwork)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Input files
# Default names are alpha_diversity.csv and metadata.csv in the current working directory.
# You may also pass paths as arguments:
# Rscript alpha_diversity.r <alpha-diversity.csv> <metadata.csv> [output_dir]
input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) == 0) {
  alpha_file <- "alpha_diversity.csv"
  metadata_file <- "metadata.csv"
  output_dir <- file.path(dirname(alpha_file), "alpha_diversity_results")
} else if (length(input_args) == 2 || length(input_args) == 3) {
  alpha_file <- input_args[1]
  metadata_file <- input_args[2]
  output_dir <- if (length(input_args) == 3) {
    input_args[3]
  } else {
    file.path(dirname(alpha_file), "alpha_diversity_results")
  }
} else {
  stop(
    "Usage: Rscript alpha-diversity.r <alpha_diversity.csv> ",
    "<metadata.csv> [output_dir]"
  )
}

if (!file.exists(alpha_file)) {
  stop("Alpha-diversity input file not found: ", alpha_file)
}
if (!file.exists(metadata_file)) {
  stop("Metadata input file not found: ", metadata_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# Read input files; readr guesses the field separator
alpha <- read_delim(
  alpha_file,
  delim = NULL,
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE,
  trim_ws = TRUE
)

metadata <- read_delim(
  metadata_file,
  delim = NULL,
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE,
  trim_ws = TRUE
)

# Keep sample IDs as supplied in alpha_diversity.csv
alpha <- alpha %>%
  mutate(Sample = as.character(SampleID))

# Join diversity data to metadata
data <- alpha %>%
  left_join(
    metadata,
    by = c("Sample" = "Sample_ID")
  )

if (any(is.na(data$Fermentation_batch))) {
  unmatched <- data$Sample[is.na(data$Fermentation_batch)]
  stop(
    "Some samples could not be matched to the metadata: ",
    paste(unmatched, collapse = ", ")
  )
}

# Derive treatment groups from the culture names present in the metadata.
# Empty cells and dash markers indicate that the corresponding culture is absent.
clean_culture_name <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "-")] <- NA_character_
  x
}

data <- data %>%
  mutate(
    Starter_culture = clean_culture_name(Starter_culture),
    Probiotic_culture = clean_culture_name(Probiotic_culture),
    Fermentation_treatment = case_when(
      !is.na(Starter_culture) & !is.na(Probiotic_culture) ~
        paste(Starter_culture, Probiotic_culture, sep = " + "),
      !is.na(Starter_culture) ~ Starter_culture,
      !is.na(Probiotic_culture) ~ Probiotic_culture,
      TRUE ~ "No culture"
    ),
    Fermentation_treatment = factor(
      Fermentation_treatment,
      levels = unique(Fermentation_treatment)
    ),
    Starter_culture = if_else(
      is.na(Starter_culture), "No starter culture", Starter_culture
    ),
    Probiotic_culture = if_else(
      is.na(Probiotic_culture), "No probiotic culture", Probiotic_culture
    )
  )

treatment_groups <- levels(droplevels(data$Fermentation_treatment))
message(
  "Fermentation treatment groups found: ",
  paste(treatment_groups, collapse = ", ")
)
# Require the two alpha-diversity metrics used by the analysis.
diversity_columns <- c("16S", "ITS")
missing_diversity_columns <- setdiff(diversity_columns, names(data))
if (length(missing_diversity_columns) > 0) {
  stop(
    "Required alpha-diversity columns are missing: ",
    paste(missing_diversity_columns, collapse = ", ")
  )
}

for (column in diversity_columns) {
  if (!is.numeric(data[[column]])) {
    data[[column]] <- readr::parse_number(
      as.character(data[[column]]),
      locale = readr::locale(decimal_mark = ",")
    )
  }
}

# Find numeric metadata features automatically. Identifier and culture columns
# are excluded because their values describe groups rather than measurements.
excluded_metadata_columns <- c(
  "Sample_ID", "Fermentation_batch", "Starter_culture", "Probiotic_culture"
)
continuous_candidates <- setdiff(names(metadata), excluded_metadata_columns)
continuous_features <- character()

for (column in continuous_candidates) {
  values <- data[[column]]
  if (inherits(values, c("Date", "POSIXct", "POSIXlt"))) next

  raw_values <- trimws(as.character(values))
  present <- !is.na(values) & nzchar(raw_values)
  if (!any(present)) next

  looks_numeric <- is.numeric(values) || all(grepl(
    "^[-+]?[0-9][0-9.,]*[[:space:]]*[[:alpha:]%°/]*$",
    raw_values[present]
  ))
  if (!looks_numeric) next

  if (is.numeric(values)) {
    parsed_values <- values
  } else {
    has_comma_decimal <- any(grepl(",", raw_values[present], fixed = TRUE))
    has_dot_decimal <- any(grepl(".", raw_values[present], fixed = TRUE))
    decimal_mark <- if (has_comma_decimal && !has_dot_decimal) "," else "."
    parsed_values <- readr::parse_number(
      raw_values,
      locale = readr::locale(decimal_mark = decimal_mark)
    )
  }

  if (all(!is.na(parsed_values[present]))) {
    data[[column]] <- parsed_values
    continuous_features <- c(continuous_features, column)
  }
}

message(
  "Numeric metadata features found: ",
  if (length(continuous_features) > 0) {
    paste(continuous_features, collapse = ", ")
  } else {
    "none"
  }
)

# Long format for plotting and statistical analysis
data_long <- data %>%
  pivot_longer(
    cols = all_of(diversity_columns),
    names_to = "Diversity_metric",
    values_to = "Shannon"
  )

categorical_features <- c(
  "Fermentation_treatment",
  "Starter_culture",
  "Probiotic_culture"
)

feature_labels <- c(
  Fermentation_treatment = "Fermentation treatment",
  Starter_culture = "Starter culture",
  Probiotic_culture = "Probiotic culture",
  Time = "Time (days)",
  pH = "pH",
  Acidity = "Acidity (g lactic acid/100 mL brine)",
  Temperature = "Temperature (°C)"
)
additional_feature_labels <- setdiff(continuous_features, names(feature_labels))
if (length(additional_feature_labels) > 0) {
  feature_labels <- c(
    feature_labels,
    stats::setNames(additional_feature_labels, additional_feature_labels)
  )
}
# Generate colors for the treatment groups found in this dataset.
batch_colors <- setNames(
  grDevices::hcl.colors(length(treatment_groups), palette = "Dark 3"),
  treatment_groups
)

# Generate culture palettes from the observed culture names.
starter_levels <- unique(as.character(data$Starter_culture))
probiotic_levels <- unique(as.character(data$Probiotic_culture))

starter_colors <- setNames(
  grDevices::hcl.colors(length(starter_levels), palette = "Dark 3"),
  starter_levels
)
starter_colors[starter_levels != "No starter culture"] <- "#7B2CBF"
starter_colors[starter_levels == "No starter culture"] <- "#D9D9D9"

probiotic_colors <- setNames(
  grDevices::hcl.colors(length(probiotic_levels), palette = "Dark 3"),
  probiotic_levels
)
probiotic_colors[probiotic_levels != "No probiotic culture"] <- "#F28E2B"
probiotic_colors[probiotic_levels == "No probiotic culture"] <- "#D9D9D9"
# Overall statistical tests
overall_results <- list()

for (metric in c("16S", "ITS")) {
  metric_data <- filter(data_long, Diversity_metric == metric)
  metric_results <- list()

  for (feature in c(categorical_features, continuous_features)) {
    valid <- metric_data %>%
      filter(!is.na(.data[[feature]]), !is.na(Shannon))

    statistic <- NA_real_
    p_value <- NA_real_

    if (feature %in% categorical_features) {
      test_name <- "Kruskal-Wallis"
      effect <- "H statistic"
      groups_for_test <- droplevels(as.factor(valid[[feature]]))
      if (nrow(valid) >= 2 && nlevels(groups_for_test) >= 2) {
        test <- kruskal.test(x = valid$Shannon, g = groups_for_test)
        statistic <- unname(test$statistic)
        p_value <- test$p.value
      }
    } else {
      test_name <- "Spearman correlation"
      effect <- "Spearman rho"
      if (
        nrow(valid) >= 3 &&
        length(unique(valid[[feature]])) > 1 &&
        length(unique(valid$Shannon)) > 1
      ) {
        test <- cor.test(
          valid[[feature]], valid$Shannon,
          method = "spearman", exact = FALSE
        )
        statistic <- unname(test$estimate)
        p_value <- test$p.value
      }
    }
    metric_results[[feature]] <- data.frame(
      Diversity_metric = metric,
      Metadata_feature = feature_labels[[feature]],
      Test = test_name,
      n = nrow(valid),
      Statistic = statistic,
      Effect_size = effect,
      p_value = p_value
    )
  }

  metric_results <- bind_rows(metric_results) %>%
    mutate(p_Holm_across_features = p.adjust(p_value, method = "holm"))

  overall_results[[metric]] <- metric_results
}

overall_statistics <- bind_rows(overall_results)

write_csv(
  overall_statistics,
  file.path(output_dir, "alpha_diversity_overall_statistics.csv")
)

# Pairwise Wilcoxon tests for categorical features
pairwise_results <- list()

for (metric in c("16S", "ITS")) {
  metric_data <- filter(data_long, Diversity_metric == metric)

  for (feature in categorical_features) {
    valid <- metric_data %>%
      filter(!is.na(.data[[feature]]), !is.na(Shannon))

    groups <- split(valid$Shannon, as.character(valid[[feature]]))
    groups <- groups[lengths(groups) >= 2]


    if (length(groups) >= 2) {
      pairs <- combn(names(groups), 2, simplify = FALSE)

      pair_results <- lapply(pairs, function(pair) {
        test <- wilcox.test(
          groups[[pair[1]]],
          groups[[pair[2]]],
          exact = FALSE
        )

        data.frame(
          Diversity_metric = metric,
          Metadata_feature = feature_labels[[feature]],
          Group_1 = pair[1],
          Group_2 = pair[2],
          Test = "Wilcoxon rank-sum",
          n_1 = length(groups[[pair[1]]]),
          n_2 = length(groups[[pair[2]]]),
          p_value = test$p.value
        )
      })

      pairwise_results[[paste(metric, feature, sep = "_")]] <-
        bind_rows(pair_results) %>%
        mutate(p_Holm_within_feature = p.adjust(p_value, method = "holm"))
    }
  }
}

pairwise_statistics <- if (length(pairwise_results) > 0) {
  bind_rows(pairwise_results)
} else {
  data.frame(
    Diversity_metric = character(), Metadata_feature = character(),
    Group_1 = character(), Group_2 = character(), Test = character(),
    n_1 = integer(), n_2 = integer(), p_value = double(),
    p_Holm_within_feature = double()
  )
}

write_csv(
  pairwise_statistics,
  file.path(output_dir, "alpha_diversity_pairwise_tests.csv")
)

# Spearman correlations for each continuous variable within each fermentation group
batch_correlation_results <- list()

for (metric in c("16S", "ITS")) {
  metric_data <- filter(data_long, Diversity_metric == metric)

  for (feature in continuous_features) {
    group_results <- lapply(treatment_groups, function(group) {
      valid <- metric_data %>%
        filter(
          Fermentation_treatment == group,
          !is.na(.data[[feature]]),
          !is.na(Shannon)
        )

      rho <- NA_real_
      p_value <- NA_real_
      if (
        nrow(valid) >= 3 &&
        length(unique(valid[[feature]])) > 1 &&
        length(unique(valid$Shannon)) > 1
      ) {
        test <- cor.test(
          valid[[feature]], valid$Shannon,
          method = "spearman", exact = FALSE
        )
        rho <- unname(test$estimate)
        p_value <- test$p.value
      }

      data.frame(
        Diversity_metric = metric,
        Metadata_feature = feature_labels[[feature]],
        Fermentation_treatment = group,
        n = nrow(valid),
        Spearman_rho = rho,
        p_value = p_value
      )
    })

    batch_correlation_results[[paste(metric, feature, sep = "_")]] <-
      bind_rows(group_results) %>%
      mutate(p_Holm_within_feature = p.adjust(p_value, method = "holm"))
  }
}

batch_correlations <- if (length(batch_correlation_results) > 0) {
  bind_rows(batch_correlation_results)
} else {
  data.frame(
    Diversity_metric = character(), Metadata_feature = character(),
    Fermentation_treatment = character(), n = integer(),
    Spearman_rho = double(), p_value = double(),
    p_Holm_within_feature = double()
  )
}

write_csv(
  batch_correlations,
  file.path(
    output_dir,
    "alpha_diversity_correlations_by_fermentation_batch.csv"
  )
)

# Function to create plots for one diversity metric
make_figure <- function(metric) {
  metric_data <- filter(data_long, Diversity_metric == metric)

  # Box plots: fermentation batch, starter culture, and probiotic culture
  boxplot_list <- lapply(categorical_features, function(feature) {
    valid <- metric_data %>%
      filter(!is.na(.data[[feature]]), !is.na(Shannon))

    group_names <- sort(unique(as.character(valid[[feature]])))
    comparisons <- if (length(group_names) >= 2) {
      combn(group_names, 2, simplify = FALSE)
    } else {
      list()
    }

    p <- ggplot(
      valid,
      aes(x = .data[[feature]], y = Shannon, fill = .data[[feature]])
    ) +
      geom_boxplot(
        width = 0.65,
        outlier.shape = NA,
        color = "#263238"
      ) +
      geom_jitter(
        width = 0.08,
        size = 1.8,
        alpha = 0.85,
        color = "#202733"
      ) +
      labs(
        title = feature_labels[[feature]],
        x = NULL,
        y = "Shannon",
        fill = feature_labels[[feature]]
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    if (feature == "Fermentation_treatment") {
      p <- p + scale_fill_manual(values = batch_colors, drop = FALSE)
    } else if (feature == "Starter_culture") {
      p <- p + scale_fill_manual(
        values = starter_colors,
        na.value = "#D9D9D9",
        drop = FALSE
      )
    } else {
      p <- p + scale_fill_manual(
        values = probiotic_colors,
        na.value = "#D9D9D9",
        drop = FALSE
      )
    }

    # Pairwise Wilcoxon significance brackets
    if (length(comparisons) > 0) {
      y_max <- max(valid$Shannon, na.rm = TRUE)
      y_range <- diff(range(valid$Shannon, na.rm = TRUE))
      if (y_range == 0) y_range <- 1

      p <- p +
        stat_compare_means(
          comparisons = comparisons,
          method = "wilcox.test",
          p.adjust.method = "holm",
          label = "p.signif",
          hide.ns = FALSE,
          label.y = y_max + y_range * (0.12 + 0.12 * seq_along(comparisons))
        ) +
        expand_limits(
          y = y_max + y_range * (0.25 + 0.12 * length(comparisons))
        )
    }

    p
  })

  # Scatter plots with treatment-specific points and regression lines
  scatter_list <- lapply(continuous_features, function(feature) {
    valid <- metric_data %>%
      filter(
        !is.na(.data[[feature]]),
        !is.na(Shannon),
        !is.na(Fermentation_treatment)
      )

    ggplot(
      valid,
      aes(
        x = .data[[feature]],
        y = Shannon,
        color = Fermentation_treatment
      )
    ) +
      geom_point(size = 2, alpha = 0.85) +
      geom_smooth(
        aes(group = Fermentation_treatment),
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linewidth = 0.8
      ) +
      scale_color_manual(values = batch_colors, drop = FALSE) +
      labs(
        title = feature_labels[[feature]],
        x = NULL,
        y = "Shannon",
        color = "Fermentation treatment"
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )
  })

  # Center the three box plots in the upper row; place four scatter plots below
  top_row <- wrap_plots(
    plot_spacer(),
    boxplot_list[[1]],
    boxplot_list[[2]],
    boxplot_list[[3]],
    plot_spacer(),
    ncol = 5
  )

  figure_title <- plot_annotation(
    title = paste(metric, "alpha diversity")
  )

  if (length(scatter_list) > 0) {
    bottom_row <- wrap_plots(scatter_list, ncol = 4)
    figure <- top_row / bottom_row +
      plot_layout(heights = c(1, 1)) +
      figure_title
  } else {
    figure <- top_row + figure_title
  }

  figure
}

# Create, save, and display the 16S and ITS figures
figure_16S <- make_figure("16S")
figure_ITS <- make_figure("ITS")

ggsave(
  file.path(output_dir, "alpha_diversity_16S_metadata_plots.png"),
  figure_16S,
  width = 13,
  height = 8.7,
  dpi = 200
)

ggsave(
  file.path(output_dir, "alpha_diversity_ITS_metadata_plots.png"),
  figure_ITS,
  width = 13,
  height = 8.7,
  dpi = 200
)

# Display figures in the RStudio Plots pane
print(figure_16S)
print(figure_ITS)

message("Figures and statistical tables saved to: ", output_dir)
