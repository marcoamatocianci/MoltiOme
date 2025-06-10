#' Create a Scatter Plot of Paired RNA vs. ATAC Values
#'
#' This function generates a scatter plot to visualize the relationship between
#' normalized ATAC-seq and RNA-seq values for features linked by a common gene.
#' It integrates annotated peak data with expression data, reshapes it into a
#' long format, and plots the paired values.
#'
#' @param gr_atac_annotated A `GRanges` object containing annotated ATAC-seq peaks.
#'   Must include metadata columns `gene_id` (the overlapping gene) and `Peaks_id`.
#'   Non-annotated peaks (where `gene_id` is NA) will be ignored.
#' @param normalized_rna_data A data frame or data.table of normalized RNA-seq data.
#'   Must have a `row.names` column containing gene IDs that match `gene_id`
#'   in `gr_atac_annotated`. Columns should correspond to samples.
#' @param normalized_atac_data A data frame or data.table of normalized ATAC-seq data.
#'   Must have a `row.names` column containing peak IDs that match `Peaks_id`
#'   in `gr_atac_annotated`. Columns should correspond to samples.
#' @param color A variable to map to the color aesthetic in `ggplot2`. This can be a
#'   static value (e.g., `"blue"`) or a column name from the reshaped data
#'   (e.g., `atac.Name` to color points by sample).
#'
#' @return A `ggplot` object showing the scatter plot of ATAC vs. RNA values.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import GenomicRanges
#' @export
FinalPlot <- function(gr_atac_annotated, normalized_rna_data, normalized_atac_data, color){
   data_plot <- left_join(as.data.frame(gr_atac_annotated)[!is.na(gr_atac_annotated$gene_id),],normalized_rna_data, by=join_by("gene_id"=="row.names")) %>% left_join(normalized_atac_data,by = join_by("Peaks_id"=="row.names"), suffix = c(".atac",".rna")) %>%
     select(contains(c("Peaks_id","gene_id",".rna",".atac")))

   data_long <- tidyr::pivot_longer(data_plot,
                                    cols = -c("Peaks_id","gene_id"))
   data_long$Name <- sapply(strsplit(data_long$name,"\\."),"[",1)
   data_long$Type <- sapply(strsplit(data_long$name,"\\."),"[",2)
   data_long$name <- NULL
   data_long <- do.call(cbind,split(data_long, data_long$Type))

   ggplot(data_long)+
     geom_point(aes(x=rna.value, y=atac.value, colour = color))+
     theme_minimal()
}
