#' Find Protein-Coding Gene IDs from GTF Annotation
#'
#' Extracts the gene IDs of all protein-coding genes from a GTF annotation object.
#'
#' @param gtf A `GRanges` object containing GTF annotation data. Must include metadata columns:
#'   - `type`: Feature type (should include "gene")
#'   - `gene_biotype`: Gene classification (should include "protein_coding")
#'   - `gene_id`: Unique gene identifiers
#'
#' @return A character vector of `gene_id`s for protein-coding genes.
#'
#' @import GenomicRanges
#' @export
findProteinCoding <- function(gtf){
  #Selecting only protein coding gene and gene features
  coding_genes <- gtf[gtf$type == "gene" &
                        gtf$gene_biotype == "protein_coding"]$gene_id
  return(coding_genes)
}
