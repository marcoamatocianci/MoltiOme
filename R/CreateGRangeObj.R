#' Create Genomic Ranges Object from Feature Data
#'
#' Constructs a `GRanges` object from genomic feature data while adding metadata columns.
#' Removes features with mitocondrial locations and cleans chromosome names.
#'
#' @param specific_feature A data.frame containing genomic features. Must contain columns:
#'   - `V4`: Chromosome (e.g., "chr1")
#'   - `V5`: Start position
#'   - `V6`: End position
#'   - `V2`: Gene symbol
#'   - `V1`: Feature ID
#'   - Optional `sum`: Aggregate value column
#' @param data_type Character string specifying data type (used for metadata column names)
#'
#' @return A `GRanges` object with genomic ranges and metadata columns containing:
#'   - `<data_type>_Symbol`: Gene symbols from `V2`
#'   - `<data_type>_id`: Feature IDs from `V1`
#'   - (Optional) `<data_type>_sum`: Values from `sum` column if present
#'
#' @import GenomicRanges
#' @import IRanges
#' @import dplyr
#' @export
CreateGRangeObj <- function (specific_feature,data_type){
  feature_sub <- specific_feature %>%
    filter(V4 != "") # mitocondrial genes have no chromosome location

  gr <- GRanges(
    seqnames = gsub("chr","",feature_sub$V4),
    ranges = IRanges(
      start = feature_sub$V5,
      end = feature_sub$V6
    ),
  )
  #adding the gene symbol
  mcols(gr)[[gsub(" ","_",paste0(data_type, "_Symbol"))]]<- feature_sub$V2
  mcols(gr)[[gsub(" ","_",paste0(data_type, "_id"))]] <- feature_sub$V1
  if(!is.null(feature_sub$sum)){
    mcols(gr)[[gsub(" ","_",paste0(data_type, "_sum"))]] <- feature_sub$sum }
  return(gr)
}
