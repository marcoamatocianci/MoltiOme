#' Plot Mean ATAC vs. RNA Values by Chromosome
#'
#' For each chromosome, this function generates a scatter plot comparing the mean
#' normalized ATAC-seq values against the mean normalized RNA-seq values for
#' features linked by a common gene.
#'
#' @param atac_data A data frame or data.table of normalized ATAC values, with
#'   a `row.names` column containing peak IDs.
#' @param expression_data A data frame or data.table of normalized expression values,
#'   with a `row.names` column containing gene IDs.
#' @param atac_features A data frame or data.table of ATAC feature metadata. Must
#'   have a `V1` column containing peak IDs.
#' @param expression_features A data frame or data.table of expression feature metadata.
#'   Must have a `V1` column containing gene IDs.
#' @param gr_atac A `GRanges` object of annotated ATAC-seq peaks. Must have
#'   metadata columns `gene_id` and `Peaks_id`.
#'
#' @return A named list of `ggplot` objects, where each element corresponds to a
#'   chromosome and contains the scatter plot.
#'
#' @import dplyr
#' @import ggplot2
#' @import GenomicRanges
#' @export
#'
LogCPMCombinedPlot <- function(atac_data, expression_data, atac_features, expression_features, gr_atac){

  atac_features <- AddRowMeansToFeature(atac_data,atac_features)
  expression_features <- AddRowMeansToFeature(expression_data,expression_features)

  p <- list()
  for (chr in gr_atac@seqnames@values) {
    data <- left_join(as.data.frame(gr_atac[gr_atac@seqnames==chr])[!is.na(gr_atac$gene_id),],expression_features, by=join_by("gene_id"=="V1")) %>%
      left_join(atac_features,by = join_by("Peaks_id"=="V1"), suffix = c(".atac",".rna"))

    p[[chr]] <- ggplot(data)+
      geom_point(aes(x=mean.rna, y=mean.atac ))+
      theme_minimal()
  }

   return(p)
}
