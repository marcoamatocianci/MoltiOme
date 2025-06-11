#' Add Row Sums from Data Table to Feature Table
#'
#' Calculates the row sums of numeric columns in a data table and merges these sums into the features table.
#'
#' @param data_table A `data.table` with a `row.names` column and numeric columns for which row sums should be calculated.
#' @param features A `data.table` with a `V1` column corresponding to feature IDs (matching `row.names` in `data_table`).
#'
#' @return A `data.table` containing all columns from `features` and a new `sum` column with the row sums from `data_table`.
#'
#' @import data.table
#' @export
AddRowSumToFeature <- function (data_table,features){
  if (!is.null(features$sum)) {
    features$sum <- NULL
  }
  matrix_sum <- data_table[, sum := rowSums(.SD), .SDcols = is.numeric]
  features_sum <- features[matrix_sum[, .(row.names, sum)], on = .(V1=row.names)]
  return(features_sum)
}
