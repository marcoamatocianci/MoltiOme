
# MoltiOme: Multi-Omics Data Integration and Visualization Toolkit




`MoltiOme` is an R package designed for the efficient integration, normalization, annotation, and visualization of multi-omics datasets, with a primary focus on gene expression (RNA-seq) and chromatin accessibility (ATAC-seq) data.

This package provides a suite of tools to handle the entire workflow, from loading raw data outputs (e.g., from Cell Ranger) to generating insightful plots that correlate different omics layers. It is built on top of popular Bioconductor and tidyverse packages like `GenomicRanges`, `data.table`, and `ggplot2` to provide a robust and familiar analysis environment.

## Key Features

-   **Load Raw Data**: Easily load sparse matrix data from standard bioinformatics outputs (`matrix.mtx`, `features.tsv`, `barcodes.tsv`).
-   **Data Preparation**: Convert sparse data into a user-friendly `data.table` format for analysis.
-   **Normalization**: Normalize count data using the standard Counts Per Million (CPM) method.
-   **Genomic Annotation**: Annotate genomic features (like ATAC-seq peaks) with overlapping protein-coding genes from a GTF file.
-   **Feature Augmentation**: Calculate summary statistics like sums and means for features across samples.
-   **Multi-Omics Integration**: Identify and combine data for genes that are common across both expression and accessibility datasets.
-   **Statistics**: Calculate summary statistics on the success of the integration, showing how many genes and peaks were successfully merged.
-   **Visualization**: Generate plots to visualize feature intensity across chromosomes or to directly compare normalized values between omics types.

## Installation

You can install the development version of `MoltiOme` from GitHub using the `devtools` package.

```r
# First, install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install MoltiOme from GitHub
devtools::install_github("marcoamatocianci/MoltiOme")
```

## Example Workflow

This example demonstrates a complete analysis pipeline using the core functions of the `MoltiOme` package.

```r
# Load the package and its dependencies
library(MoltiOme)
library(Matrix)
library(data.table)
library(GenomicRanges)
```

### 1. Data Loading and Preparation

Start by loading your raw data and converting it into an analysis-ready `data.table`.

```r
# Assume you have your raw files: raw_mtx, features, barcodes
# For this example, we'll create simulated data
raw_mtx <- Matrix(c(150, 2000, 50, 5, 200, 2500, 65, 0), nrow = 4, sparse = TRUE)
features <- data.table(V1 = c("GENE_1", "GENE_2", "PEAK_101", "PEAK_102"), V2 = c("TP53", "EGFR", "Peak101", "Peak102"), V3 = c("RNA", "RNA", "ATAC", "ATAC"), V4 = c("chr17", "chr7", "chr17", "chr7"), V5 = c(7661779, 55019017, 7660000, 55018000), V6 = c(7687550, 55211628, 7660500, 55018500))
barcodes <- data.table(V1 = c("SampleA", "SampleB"))

# Add names to the matrix
named_mtx <- CollapseSparseCell(raw_mtx, features, barcodes)

# Convert to a data.table
master_dt <- SparsetoDataTable(named_mtx)
```

### 2. Data Extraction, Normalization, and Annotation

Next, separate the data by omics type, normalize it, and annotate your genomic features.

```r
# Extract RNA and ATAC data
rna_data <- ExtractOmicsData(master_dt, features, "RNA")
atac_data <- ExtractOmicsData(master_dt, features, "ATAC")

# Normalize the RNA data
normalized_rna <- NormalizeCPM(rna_data)
# (You would typically normalize ATAC data as well, e.g., using TF-IDF or CPM)

# Augment feature tables with summary stats for plotting
rna_features <- features[V3 == "RNA"]
atac_features <- features[V3 == "ATAC"]
rna_features <- AddRowMeansToFeature(rna_data, rna_features)
atac_features <- AddRowMeansToFeature(atac_data, atac_features)

# Create GRanges objects
gr_expression <- CreateGRangeObj(rna_features, "Gene_Expression")
gr_atac <- CreateGRangeObj(atac_features, "Peaks")

# Annotate ATAC peaks with gene information (requires a GTF object)
# gtf <- rtracklayer::import("path/to/your/genes.gtf")
# gr_atac_annotated <- MapProteinCodingOverlaps(gtf, gr_atac)
```

### 3. Integration and Final Visualization

Finally, combine the data and create a summary plot.

```r
# For this example, let's assume gr_atac_annotated exists
# gr_atac_annotated <- ... (from previous step)
# normalized_atac_data <- ... (from previous step)

# Generate a final plot comparing the two omics types
# final_plot <- FinalPlot(
#   gr_atac_annotated = gr_atac_annotated,
#   normalized_rna_data = normalized_rna,
#   normalized_atac_data = normalized_atac_data,
#   color = atac.Sample # Color points by sample name
# )
# print(final_plot)
```

## Detailed Tutorial

For a more detailed, step-by-step guide on how to use the package, please see the package vignette.

```r
# This will open the detailed tutorial in your browser
vignette("MoltiOme-workflow", package = "MoltiOme")
```

## Function Overview

| Function                 | Description                                                                          |
| ------------------------ | ------------------------------------------------------------------------------------ |
| `CollapseSparseCell`     | Assigns feature and barcode names to a raw sparse matrix.                            |
| `SparsetoDataTable`      | Converts a named sparse matrix into a dense `data.table`.                            |
| `ExtractOmicsData`       | Subsets the main data table to get data for a specific omics type (e.g., "RNA").      |
| `NormalizeCPM`           | Performs log2(CPM+1) normalization on a count table.                                 |
| `AddRowSumToFeature`     | Calculates row sums (total counts) and adds them to the features table.              |
| `AddRowMeansToFeature`   | Calculates row means (average counts) and adds them to the features table.           |
| `CreateGRangeObj`        | Creates a `GRanges` object from a features table for genomic analysis.               |
| `MapProteinCodingOverlaps`| Annotates a `GRanges` object (e.g., ATAC peaks) with overlapping gene IDs.           |
| `findProteinCoding`      | Extracts the IDs of all protein-coding genes from a GTF annotation object.           |
| `CommonGenesData`        | Combines expression and ATAC data for genes present in both datasets.                |
| `PeakStats` / `GeneStats`| Calculate statistics on how many features were successfully integrated.                |
| `FinalPlot`              | Generates a scatter plot to directly compare normalized RNA and ATAC values.         |
| `LogCPMCombinedPlot`     | Creates scatter plots of mean ATAC vs. RNA values, separated by chromosome.          |
| `PlotATACDistribution`   | Plots the intensity of ATAC peaks across genomic positions for each chromosome.      |
| `PlotGeneExpression`     | Plots the intensity of gene expression across genomic positions for each chromosome. |

## License

This package is licensed under the MIT License.

[1] https://github.com/TransBioInfoLab/pathwayMultiomics
[2] https://github.com/PacktPublishing/R-Bioinformatics-Cookbook/blob/master/README.md
[3] https://cran.r-project.org/web/packages/BGData/readme/README.html
[4] https://cran.r-project.org/web/packages/rOCEAN/readme/README.html
[5] https://ucdavis-bioinformatics-training.github.io/2022_February_Introduction_to_R_for_Bioinformatics/installing-packages.html
[6] https://github.com/mikelove/awesome-multi-omics/blob/master/README.md
[7] https://cran.r-project.org/web/packages/epiomics/readme/README.html
[8] https://rdrr.io/bioc/IMAS/man/IMAS-package.html
[9] https://du-bii.github.io/module-3-Stat-R/stat-R_2021/tutorials/data-exploration_pavkovic_2019/tuto_data-exploration_pavkovic.html
[10] https://cran.r-project.org/web/packages/qgg/readme/README.html
