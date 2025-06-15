
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
## Docker Analysis Workflow

To run the `MoltiOme` analysis pipeline within a reproducible Docker environment, you can use the pre-built image available on Docker Hub. This method ensures that all dependencies are correctly managed and allows you to run the analysis on any system with Docker installed, which is ideal for complex genomic analysis workflows[1][2][3].

### Prerequisites

1.  **Docker**: Ensure Docker is installed and running on your system.
2.  **Data Directory**: Create a main project directory on your local machine. Inside this directory, create two subfolders:
    *   `data`: Place your Cell Ranger output folder (`filtered_feature_bc_matrix`) and your GTF annotation file (`Homo_sapiens.GRCh38.114.gtf.gz`) inside this `data` folder.
    *   `outputs`: This empty folder will be used to store the results generated by the analysis.

Your directory structure should look like this:
```
my_project/
├── data/
│   ├── filtered_feature_bc_matrix/
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── Homo_sapiens.GRCh38.114.gtf.gz
└── outputs/
```

### Running the Analysis

Navigate to your main project directory (`my_project`) in your terminal and execute the following command. This will download the `mcianci/examproject` image, run the analysis using the R scripts inside, and save the output files to your local `outputs` folder.

```bash
docker run --rm \
  -v "$(pwd)/data":/data \
  -v "$(pwd)/outputs":/results \
  mcianci/examproject
```

### Command Explanation

*   `docker run --rm`: This command runs the container and includes the `--rm` flag to automatically remove it once the process is complete, keeping your system clean.
*   `-v "$(pwd)/data":/data`: This is a **bind mount** that shares data from the host to the container[8][9]. It maps your local `data` directory (referenced by the absolute path `$(pwd)/data`) to the `/data` directory inside the container. This allows the R script inside the container to read your input files[8].
*   `-v "$(pwd)/outputs":/results`: This second bind mount maps your local `outputs` directory to the `/results` directory inside the container. Any files saved to `/results` by the script will appear in your local `outputs` folder, making the results accessible on your host machine after the container finishes[8][9].
*   `mcianci/examproject`: This is the name of the Docker image that contains the `MoltiOme` package, all its R dependencies, and the script to execute the analysis pipeline.

After the container finishes running, your `outputs` directory will contain the generated plots and data files from the analysis.

## Detailed Tutorial

For a more detailed, step-by-step guide on how to use the package, please see the package vignette.

```r
# This will open the detailed tutorial in your browser
vignette("helpMoltiOme", package = "MoltiOme")
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
| `LogCPMCombinedPlot`     | Creates scatter plots of mean ATAC vs. RNA values, separated by chromosome.          |
| `PlotATACDistribution`   | Plots the intensity of ATAC peaks across genomic positions for each chromosome.      |
| `PlotGeneExpression`     | Plots the intensity of gene expression across genomic positions for each chromosome. |


