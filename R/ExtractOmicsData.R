#' Extract Omics Data by Feature Type
#'
#' Extracts a subset of rows from a data table based on feature IDs.
#'
#' @param data_table A data.frame or data.table containing omics data. Must have a column named `row.names` with feature IDs.
#' @param features A data.frame with at least column `V3`.
#' @param feature_id A character string specifying the type of feature to extract (matching `V3` in `features`).
#'
#' @return A subset of `data_table` containing only the rows corresponding to the selected features.
#'
#' @import data.table
#' @export
ExtractOmicsData <- function(data_table,features,feature_id){
  data <- data_table[data_table$row.names %in% features[features$V3 == feature_id,]$V1,]
  return(data)
}
