#' Assign Feature and Barcode Names to a Sparse Matrix
#'
#' This function assigns row and column names to a sparse matrix, typically from single-cell data,
#' using provided feature and barcode vectors.
#' @param mtx A sparse matrix (e.g., from Matrix package) to be annotated.
#' @param features A data.frame or data.table containing feature (gene/peak) identifiers in the first column (`V1`).
#' @param barcodes A vector or list containing barcode (cell) identifiers. If a list, the first element is used.
#'
#' @return The input matrix with row and column names assigned.
#' @examples
#' \dontrun{
#' # Example usage:
#' # mtx <- Matrix::Matrix(data, sparse=TRUE)
#' # features <- data.frame(V1 = c("Gene1", "Gene2"))
#' # barcodes <- list(c("Cell1", "Cell2"))
#' # mtx_named <- CollapseSparseCell(mtx, features, barcodes)
#' }
#'
#' @export
CollapseSparseCell <- function (mtx, features, barcodes){
  # Map rownames and colnames
  rownames(mtx) <- features$V1
  colnames(mtx) <- barcodes[[1]]
  return(mtx)
}
