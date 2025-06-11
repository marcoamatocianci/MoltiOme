#' Plot ATAC-seq Peak Intensity Distribution for Top Chromosomes
#'
#' Generates scatter plots of ATAC-seq peak intensity distribution across genomic positions (chromosomes in a GRanges object.)
#'
#' @param gr_atac A `GRanges` object for ATAC-seq peaks, with metadata column `Peaks_sum` (peak intensity).
#'
#' @return A named list of `ggplot2` objects, one per chromosome, visualizing peak intensity by position and size.
#'
#' @import GenomicRanges
#' @import ggplot2
#' @import scales
#' @export
PlotATACDistribution <- function(gr_atac){
  plot_list <- list()
  for (chr in gr_atac@seqnames@values){
    grsub <- gr_atac[gr_atac@seqnames==chr]
    df <- data.frame(start=grsub@ranges@start, end=(grsub@ranges@width+grsub@ranges@start-1) ,size=grsub@ranges@width, intensity=mcols(grsub)$Peaks_sum)

    p <- ggplot(df, aes(x = start, y = intensity, color = size)) +
      geom_point(alpha = 0.7, size=0.2) +
      scale_color_gradient(low = "#72A0C1", high = "#C51E3A", name = "Intensity") +
      scale_x_continuous(labels = scales::comma) +
      labs(
        title = paste("Peak Intensity Distribution:", chr),
        x = "Genomic Position",
        y = "Peak Intensity"
      ) +
      theme_minimal()

    plot_list[[chr]] <- p
  }
  return(plot_list)
}
