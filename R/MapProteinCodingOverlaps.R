#' Annotate Genomic Features with Overlapping Protein-Coding Genes
#'
#' Identifies overlaps between genomic features and protein-coding genes from GTF annotation,
#' then adds gene IDs to feature metadata.
#'
#' @param gtf A `GRanges` object containing GTF annotation data. Must includ e metadata columns:
#'   - `type`: Feature type (e.g., "gene")
#'   - `gene_biotype`: Gene classification (e.g., "protein_coding")
#'   - `gene_id`: Unique gene identifiers
#' @param gr_feature A `GRanges` object containing genomic features to annotate (e.g., peaks, regions)
#'
#' @return The input `gr_feature` with added `gene_id` metadata column indicating overlapping
#'   protein-coding genes. Non-overlapping features receive `NA`.
#'
#' @import GenomicRanges
#' @import IRanges
#' @export
MapProteinCodingOverlaps <- function(gtf,gr_feature){
  #Selecting only protein coding gene and gene features
  gtf_coding <- gtf[gtf$type == "gene" &
                      gtf$gene_biotype == "protein_coding"]
  #Finding overlaps between feature and protein coding genes
  overlaps <- findOverlaps(gr_feature, gtf_coding)

  #Annotating the Peaks with the Genes for which they overlap
  mcols(gr_feature)$gene_id <- NA
  mcols(gr_feature)$gene_id[queryHits(overlaps)] <- gtf_coding$gene_id[subjectHits(overlaps)]

  return(gr_feature)
}
