# FAIRferm-website

FAIR-compliant workflows and bioinformatics pipelines for time-series analysis of microbiome and metabolomics data in plant-based food fermentation

### FAIRferm analysis scripts

This repository contains four R scripts for alpha diversity, beta diversity, metabolomics, and metagenomics. The analyses read input files from a data directory, detect the available metadata variables, and save figures and tables in result folders.

## 1. Requirements

These instructions are for Ubuntu or another Debian-based Linux distribution. Install R and the system libraries commonly needed to compile R packages from CRAN:

```bash
sudo apt update
sudo apt install -y r-base r-base-dev   build-essential gfortran   libcurl4-openssl-dev libssl-dev libxml2-dev   libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

Check that R is available:

```bash
R --version
Rscript --version
```

The scripts install their own R packages when needed. You can also install all packages in advance:

```bash
Rscript -e 'install.packages(c(
  "ggplot2", "ggpubr", "patchwork", "readr", "dplyr", "tidyr",
  "stringr", "vegan", "gridExtra", "readxl", "pheatmap"
), repos = "<https://cloud.r-project.org>")'
```

The first run needs an internet connection to download any missing R packages.

## 2. Repository and input-file layout

Keep the scripts in a `scripts/` directory and the data files together in a `data/` directory. For example:

```
your-project/
├── data/
│   ├── 16S_tax_table.csv
│   ├── ITS*_tax_table.csv
│   ├── alpha_diversity.csv
│   ├── metabolites.xlsx
│   └── metadata.csv
└── scripts/
    ├── alpha-diversity.r
    ├── beta-diversity.r
    ├── metabolomics_heatmap.r 
    └── metagenomics_heatmap.r 
```

Only include the files needed for the analyses you intend to run. The beta-diversity analysis requires both tax tables. The metagenomics heatmap can run with the 16S table alone; the ITS table is optional.

Linux filenames are case-sensitive. The commands below use the script filenames you specified; make sure the committed filenames match exactly, including the lower-case `.r` extension.

## 3. Input data formats

### Metadata:
`metadata.csv`

The metadata file should have one row per sample. The separator is detected automatically by the scripts that read metadata. Sample IDs must correspond to the IDs in the analysis data.

For alpha diversity, the metadata must contain `Sample_ID`, `Fermentation_batch`, `Starter_culture`, and `Probiotic_culture`. The beta-diversity and metagenomics scripts accept a sample-ID column named `Sample_Name`, `Sample_ID`, `SampleID`, or `Sample`; beta diversity also requires `Fermentation_batch`. The metabolomics script detects common sample-ID column names and includes all other metadata columns as annotations. Numeric features are inferred from their values; other fields are treated as categories.

### Alpha diversity: 
`alpha_diversity.csv`

Provide the columns:

- `SampleID`: sample IDs matching the `Sample_ID` column in metadata.
- `16S`: 16S alpha-diversity values.
- `ITS`: ITS alpha-diversity values.

The script also uses `Fermentation_batch`, `Starter_culture`, and `Probiotic_culture` from the metadata to form group comparisons. Numeric metadata features are detected automatically and included in statistical analyses and scatterplots.

### Tax tables: 
`16S_tax_table.csv`, `ITS_tax_table.csv`

Each tax table should contain taxa in rows: the first column contains taxon names, and every remaining column is a sample. Sample-column names must match metadata sample IDs.

The beta-diversity script expects both the 16S and ITS tax tables. The metagenomics script requires a file beginning with `16S` and ending in `_tax_table.csv`; it also processes an ITS file matching `ITS*_tax_table.csv` when one is present. For example, `ITS1_tax_table.csv` is accepted.

### Metabolomics: 
`metabolites.xlsx`

The first worksheet should have metabolite names in the first column and sample measurements in the remaining columns. The script averages replicate columns whose names end with `-i`, `-ii`, or `-iii`, then matches the resulting sample names to metadata. Positive measurements are log2-transformed and standardized by metabolite for the heatmap.

## 4. Run the analyses

From the repository root, run the scripts with the data directory as shown.

### Alpha diversity

```bash
Rscript scripts/alpha-diversity.r \
  data/alpha_diversity.csv \
  data/metadata.csv \
  data/alpha_diversity_results
```

The output folder contains:

- `alpha_diversity_16S_metadata_plots.png`
- `alpha_diversity_ITS_metadata_plots.png`
- Overall statistical results, pairwise tests, and correlations by fermentation treatment in CSV files.

### Beta diversity

```bash
Rscript scripts/beta-diversity.r data
```

The script infers treatment labels from the sample-name prefix before the first number (for example, `FBSP1a` becomes `FBSP`). Use consistent sample names that begin with the treatment prefix followed by a number.

The `data/beta_diversity_results/` folder contains 16S and ITS Bray–Curtis PCoA figures (PNG and PDF), distance matrices, coordinates, eigenvalues, PERMANOVA and PERMDISP results, and correlations by treatment when numeric metadata features are found.

### Metabolomics heatmap

```bash
Rscript scripts/metabolomics_heatmap.r data
```

The script reads the Excel workbook and metadata from `data/`, detects the metadata sample-ID column and variables, and saves:

```
data/metabolomics_heatmap_results/
└── metabolite_heatmap_clustered_metadata.png
```

Keep only one metabolite `.xlsx` workbook in the data directory so the script can identify it unambiguously.

### Metagenomics heatmaps

```bash
Rscript scripts/metagenomics.heatmap.r data
```

The script detects the 16S tax table and, when present, the ITS tax table. Treatment labels are inferred from the sample-name prefix before the first number. It detects numeric and categorical metadata features and includes them as heatmap annotations.

The `data/metagenomics_heatmap_results/` folder contains:

- `16S_taxonomy_heatmap_top100_clustered_metadata.png`
- `ITS_taxonomy_heatmap_top100_clustered_metadata.png` (when an ITS table is present)

## 5. Run all four from the Linux terminal

After placing all required inputs in `data/`, run:

```bash
Rscript scripts/alpha-diversity.r data/alpha_diversity.csv data/metadata.csv data/alpha_diversity_results
Rscript scripts/beta-diversity.r data
Rscript scripts/metabolomics_heatmap.r data
Rscript scripts/metagenomics.heatmap.r data
```

Each script reports the location of its output files when it finishes. Check that sample IDs match across the relevant input files if an analysis reports missing or unmatched samples.

## 6. Run the analyses in RStudio

1. Open the repository as an RStudio Project. If the repository does not have an `.Rproj` file, open RStudio and choose **File → New Project → Existing Directory**, then select the repository folder.
2. Set the working directory to the `data/` folder. In RStudio, choose **Session → Set Working Directory → Choose Directory...** and select `data/`. Alternatively, run this in the Console from the repository root:
    
    ```r
    setwd("data")
    getwd()  # confirm that this is the data folder
    ```
    
3. Run one script at a time from the RStudio Console. Since the working directory is `data/`, use the relative paths below:
    
    ```r
    source("../scripts/alpha-diversity.r")
    source("../scripts/beta-diversity.r")
    source("../scripts/metabolomics_heatmap.r")
    source("../scripts/metagenomics.heatmap.r")
    ```
    
    You can also open a script in the editor and click **Source**. Run the alpha-diversity script only after `alpha_diversity.csv` and `metadata.csv` are in the working directory. The other three scripts also read their input files from the working directory.
    
4. The scripts install any missing R packages from CRAN on first use. Accept the installation prompts if RStudio asks for permission. When installation finishes, run the script again if it did not continue automatically.
5. Find the generated files in the result folders inside `data/`, such as `data/beta_diversity_results/` and `data/metagenomics_heatmap_results/`. Use the **Files** pane in RStudio to open the folders or view the PNG files.

RStudio's **Source** button runs a script in the current R session, so its working directory determines where the script looks for input data. If the script reports that an input file is missing, check the working directory with `getwd()` and confirm that the filenames and capitalization match exactly.
