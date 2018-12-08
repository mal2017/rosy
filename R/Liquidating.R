#' Find signal from a bam within regions.
#'
#' @param gr A GRList object.
#' @param reads A bam file path, GRanges list, or other object holding reads.
#' @return RangedSummarizedExperiment
liquidate_internal <- function(features, reads, paired = F) {
  res <- GenomicAlignments::summarizeOverlaps(features = features,
                                       reads = reads,
                                       inter.feature = T,
                                       singleEnd = !paired)
}

#' Scale signal in RangedSummarizedExperiment.
#'
#' @param rse A RangedSummarizedExperiment
#' @return RangedSummarizedExperiment
normalize_internal <- function(rse) {
  mat <- BiocGenerics::counts(rse)
  norm_factor <- 1/(colSums(mat)/1e6)
  BiocGenerics::counts(rse) <- t(t(mat)*norm_factor)
  rse
}
