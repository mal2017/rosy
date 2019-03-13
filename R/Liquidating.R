#' Find signal from a bam within regions.
#'
#' Generic behavior is handled by GenomicAlignments. This is a thin wrapper.
#'
#' @param gr A GRList object.
#' @param reads A bam file path, GRanges list, or other object holding reads.
#' @return RangedSummarizedExperiment
liquidate_internal <- function(features, reads, paired = F) {
  GenomicAlignments::summarizeOverlaps(features = features,
                                       reads = reads,
                                       inter.feature = T,
                                       singleEnd = !paired)
}

#' Scale signal in RangedSummarizedExperiment by MRIP.
#'
#' @param rse A RangedSummarizedExperiment
#' @return RangedSummarizedExperiment
normalize_mrip_internal <- function(rse,assay="counts") {
  if ("readInPeakNorm" %in% names(SummarizedExperiment::assays(rse))) {
    stop("This normalization has alreadby been performed.")
  }
  mat <- SummarizedExperiment::assay(rse,i = assay)
  norm_factor <- 1/(colSums(mat)/1e6)
  SummarizedExperiment::assay(rse,i = "readInPeakNorm") <- t(t(mat)*norm_factor)
  rse
}

#' Scale signal in RangedSummarizedExperiment per BP.
#'
#' @param rse A RangedSummarizedExperiment
#' @return RangedSummarizedExperiment
normalize_perbp_internal <- function(rse, assay="readInPeakNorm") {
  if ("readInPeakNormPerBP" %in% names(SummarizedExperiment::assays(rse))) {
    stop("This normalization has alreadby been performed.")
  }
  mat <- SummarizedExperiment::assay(rse,i = assay)
  norm_factor <- 1/unlist(lapply(GenomicRanges::width(SummarizedExperiment::rowRanges(rse)),sum))
  SummarizedExperiment::assay(rse,i = "readInPeakNormPerBP") <- mat*norm_factor
  rse
}
