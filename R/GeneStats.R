#' Calculate Gene Annotation Statistics
#'
#' Computes statistics on the overlap between gene expression features and ATAC-seq peaks.
#'
#' @param expression A  GRanges object with a `Gene_Expression_id` column for gene expression features.
#' @param atac A GRanges object with a `gene_id` column for ATAC-seq features.
#'
#' @return A data.frame with one row and the following columns:
#'   - `tot_genes`: Total number of genes in the expression data
#'   - `merged_num`: Number of genes with at least one overlapping ATAC peak
#'   - `non_merged`: Number of genes without any overlapping ATAC peak
#'   - `prc_merged`: Percentage of merged genes out of total
#'   - `prc_unmerged`: Percentage of non-merged genes out of total
#'
#' @import data.table
#' @export
GeneStats <- function(expression,atac){
  common_genes <- intersect(expression$Gene_Expression_id, atac$gene_id)
  GenesNotMerged <- expression[!expression$Gene_Expression_id %in% common_genes]$Gene_Expression_id

  genes_stats <- data.frame(row.names = "gene_stats",tot_genes=length(expression$Gene_Expression_id), merged_num=length(common_genes),non_merged=length(GenesNotMerged),prc_merged=length(common_genes)/length(expression$Gene_Expression_id)*100, prc_unmerged=length(GenesNotMerged)/length(expression$Gene_Expression_id)*100)

  return(genes_stats)
}
