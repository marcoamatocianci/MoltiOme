#' Plot Gene Expression Intensity Distribution by Chromosome
#'
#' Generates scatter plots of gene expression intensity across genomic positions
#' for each chromosome in a GRanges object.
#'
#' @param gr_expression A `GRanges` object for gene expression features, with metadata column `Gene_Expression_sum` (expression intensity).
#'
#' @return A named list of `ggplot2` objects, one per chromosome, visualizing expression intensity by position.
#'
#' @import GenomicRanges
#' @import ggplot2
#' @import scales
#' @export
PlotGeneExpression <- function(gr_expression){

  plot_list <- list()
  for (chr in gr_expression_coding@seqnames@values){
    grsub <- gr_expression_coding[gr_expression_coding@seqnames==chr]
    df <- data.frame(start=grsub@ranges@start, end=(grsub@ranges@width+grsub@ranges@start-1) ,size=grsub@ranges@width, intensity=mcols(grsub)$Gene_Expression_sum)

    p <- ggplot(df, aes(x = start, y = intensity, color = intensity)) +
      geom_point(alpha = 0.7, size=1) +
      scale_color_gradient(low = "#E9D66B", high = "#856088", name = "Intensity") +
      scale_x_continuous(labels = scales::comma) +
      labs(
        title = paste("Expression Intensity Distribution:", chr),
        x = "Genomic Position",
        y = "Expression Intensity"
      ) +
      theme_minimal()

    plot_list[[chr]] <- p
  }

}
