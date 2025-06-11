#' Convert a Sparse Matrix to a data.table
#'
#' This function converts a sparse matrix (e.g., from the Matrix package) to a \code{data.table},
#' preserving row names as a separate column named \code{row.names}.
#'
#' @param mtx A sparse matrix to be converted.
#'
#' @return A \code{data.table} with the same data as the input matrix. The row names are stored in a column called \code{row.names}.
#'
#' @examples
#' \dontrun{
#' # mtx <- Matrix::Matrix(data, sparse=TRUE)
#' # dt <- SparsetoDataTable(mtx)
#'}
#' @import data.table
#' @export
SparsetoDataTable <- function(mtx){
  full_mat <- as.matrix(mtx)
  data_table <- as.data.table(full_mat,keep.rownames = "row.names")
  return(data_table)
}
