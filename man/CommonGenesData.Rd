% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonGenesData.R
\name{CommonGenesData}
\alias{CommonGenesData}
\title{Extract Data for Genes with Both Expression and ATAC Overlap}
\usage{
CommonGenesData(data_table, expression_data, atac_data, gr_expression, gr_atac)
}
\arguments{
\item{data_table}{A data.frame or data.table containing all features (not used in current implementation).}

\item{expression_data}{A data.frame or data.table with gene expression data. Must have a `row.names` column with gene IDs.}

\item{atac_data}{A data.frame or data.table with ATAC-seq data. Must have a `row.names` column with peak IDs.}

\item{gr_expression}{A `GRanges` object for expression features. Must have a `Gene_Expression_id` metadata column.}

\item{gr_atac}{A `GRanges` object for ATAC features. Must have a `gene_id` (overlapping gene) and `Peaks_id` metadata columns.}
}
\value{
A data.frame or data.table containing rows from both `expression_data` (for common genes) and
        `atac_data` (for peaks linked to those genes), combined using `rbind`.

Combined expression and ATAC data for genes with both data types.
}
\description{
Identifies genes present in both expression and ATAC data (based on overlaps in GRanges objects),
and returns the corresponding rows from the expression and ATAC data tables.
}
