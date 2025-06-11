#' Extract Data for Genes with Both Expression and ATAC Overlap
#'
#' Identifies genes present in both expression and ATAC data (based on overlaps in GRanges objects),
#' and returns the corresponding rows from the expression and ATAC data tables.
#'
#' @param data_table A data.frame or data.table containing all features (not used in current implementation).
#' @param expression_data A data.frame or data.table with gene expression data. Must have a `row.names` column with gene IDs.
#' @param atac_data A data.frame or data.table with ATAC-seq data. Must have a `row.names` column with peak IDs.
#' @param gr_expression A `GRanges` object for expression features. Must have a `Gene_Expression_id` metadata column.
#' @param gr_atac A `GRanges` object for ATAC features. Must have a `gene_id` (overlapping gene) and `Peaks_id` metadata columns.
#'
#' @return A data.frame or data.table containing rows from both `expression_data` (for common genes) and
#'         `atac_data` (for peaks linked to those genes), combined using `rbind`.
#'
#' @import data.table
#' @return Combined expression and ATAC data for genes with both data types.
#' @export
CommonGenesData <- function(data_table,expression_data,atac_data,gr_expression,gr_atac){

  common_genes <- intersect(gr_expression$Gene_Expression_id, gr_atac$gene_id)
  final_peaks <- gr_atac[gr_atac$gene_id %in% common_genes]$Peaks_id
  common1 <- rbind(expression_data[row.names %in% common_genes,],atac_data[row.names %in% final_peaks,])
  # common2 <- data_table[row.names %in% union(final_peaks,common_genes),]
  return(common1)
}
