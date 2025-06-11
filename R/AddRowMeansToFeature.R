#' Add Row Mean from Data Table to Feature Table
#'
#' Calculates the row Mean of numeric columns in a data table and merges these mean into the features table.
#'
#' @param data_table A `data.table` with a `row.names` column and numeric columns for which row mean should be calculated.
#' @param features A `data.table` with a `V1` column corresponding to feature IDs (matching `row.names` in `data_table`).
#'
#' @return A `data.table` containing all columns from `features` and a new `sum` column with the row means from `data_table`.
#'
#' @import data.table
#' @export
AddRowMeansToFeature <- function (data_table,features){
  if (!is.null(features$mean)) {
    features$mean <- NULL
  }
  matrix_mean<- data_table[, mean := rowMeans(.SD), .SDcols = is.numeric]
  features_mean <- features[matrix_mean[, .(row.names, mean)], on = .(V1=row.names)]
  return(features_mean)
}
