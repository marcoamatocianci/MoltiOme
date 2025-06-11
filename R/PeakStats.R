#' Calculate Peak Annotation Statistics
#'
#' Computes statistics on the overlap between ATAC-seq peaks and gene expression features.
#'
#' @param expression A data.frame or data.table (or GRanges) with a `Gene_Expression_id` column for gene expression features.
#' @param atac A data.frame or data.table (or GRanges) with `gene_id` (overlapping gene) and `Peaks_id` columns for ATAC-seq features.
#'
#' @return A data.frame with one row and the following columns:
#'   - `tot_peaks`: Total number of ATAC peaks
#'   - `merged_num`: Number of peaks overlapping genes with expression data
#'   - `non_merged`: Number of peaks not overlapping any expressed gene
#'   - `prc_tot`: Percentage of merged peaks out of total
#'   - `prc_unmerged`: Percentage of non-merged peaks out of total
#'
#' @import data.table
#' @export
PeakStats <- function(expression,atac){
  common_genes <- intersect(expression$Gene_Expression_id, atac$gene_id)
  PeakNotMerged <- atac[!atac$gene_id %in% common_genes]$Peaks_id
  final_peaks <- atac[atac$gene_id %in% common_genes]$Peaks_id
  total_peaks <- atac$gene_id
  final_genes <- common_genes
  peak_stats <- data.frame(row.names = "peak_stats",tot_peaks=length(total_peaks),merged_num=length(final_peaks),non_merged=length(PeakNotMerged),prc_merged=length(final_peaks)/length(total_peaks)*100,prc_unmerged=length(PeakNotMerged)/length(total_peaks)*100)

  return(peak_stats)

}
