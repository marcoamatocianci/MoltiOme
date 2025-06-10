#' Normalize Count Data to CPM (Counts Per Million) and Log2 Transform
#'
#' Converts a count matrix to Counts Per Million (CPM), applies log2 transformation, and preserves row names.
#'
#' @param data_table A data.frame or data.table with a column named `row.names` containing feature identifiers and count columns for samples.
#'
#' @return A data.table with normalized (log2 CPM) values and a `row.names` column for feature IDs.
#'
#' @import data.table
#' @export
NormalizeCPM <- function(data_table){
  rownames <- data_table$row.names
  data_table$row.names <- NULL
  if(!is.null(data_table$sum)){
    data_table$sum <- NULL
    warning("sum column is excluded!")
  }
  if(!is.null(  data_table$mean)){
    data_table$mean<- NULL
    warning("mean column is excluded!")
  }
  data_table_frac <- sweep(data_table,2,colSums(data_table),"/")
  data_table_norm <- log2((data_table_frac*10^6)+1)
  rownames(data_table_norm) <- rownames
  return(setDT(data_table_norm,keep.rownames = "row.names"))
}
