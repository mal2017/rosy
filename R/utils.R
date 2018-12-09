#' Combine RangedSummarizedExperiments
#' @rdname bindRseCols
#' @param rses A list of RangedSummarizedExperiments
#' @export
bindRseCols <- function(rses) {

  if (!is.list(rses)) {
    rses <- list(rses)
  }

  stopifnot(all(lapply(rses,
                       FUN = function(x) equivalentRanges(rses[[1]], x))))

  mat <- do.call(cbind, lapply(rses, FUN = SummarizedExperiment::assay))

  colnames(mat) <- make.unique(colnames(mat))

  SummarizedExperiment(list(counts = mat),
                       rowRanges = SummarizedExperiment::rowRanges(rses[[1]]))
}

#' Check two RangedSummarizedExperiments use the same ranges
#' @rdname equivalentRanges
#' @param rse1 A RangedSummarizedExperiment
#' @param rse2 A RangedSummarizedExperiment
#' @export
equivalentRanges <- function(rse1, rse2) {
  identical(SummarizedExperiment::rowRanges(rse1),
            SummarizedExperiment::rowRanges(rse2))
}
