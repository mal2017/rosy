#' Combine RangedSummarizedExperiments
#' @rdname bindRseCols
#' @param rse1 RangedSummarizedExperiment
#' @param rse2 RangedSummarizedExperiment
#' @export
bindRseCols <- function(rse1, rse2) {

  stopifnot(!anyDuplicated(c(colnames(rse1),colnames(rse2))))

  stopifnot(equivalentAssays(rse1, rse2))

  mats <- lapply(names(SummarizedExperiment::assays(rse1)),
                 FUN = function(x) {
                   cbind(SummarizedExperiment::assay(rse1, i=x),
                         SummarizedExperiment::assay(rse2, i=x))
          })

  names(mats) <- names(SummarizedExperiment::assays(rse1))

  SummarizedExperiment::SummarizedExperiment(mats,
                 rowRanges = SummarizedExperiment::rowRanges(rse1))
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

#' Check two RangedSummarizedExperiments have the same
#' assays.
#' @rdname equivalentAssays
#' @param rse1 A RangedSummarizedExperiment
#' @param rse2 A RangedSummarizedExperiment
#' @export
equivalentAssays <- function(rse1, rse2) {
  equivalentRanges(rse1,rse2) &
    identical(names(SummarizedExperiment::assays(rse1)),
              names(SummarizedExperiment::assays(rse2)))
}

#' Rename GRangesList in human readable coords.
seGrlNames <-  function(grl) {
  lapply(grl, FUN=function(x) {
    paste0(getChrNameFromGr(x),
           ":",
           min(GenomicRanges::start(x)),
           "-",
           max(GenomicRanges::end(x)))
  })
}

#' Get 1 chr name from a GR and confirm all are from same level.
getChrNameFromGr <- function(gr) {
  sn <- unique(as.character(GenomicRanges::seqnames(gr)))

  stopifnot(length(sn) == 1)

  sn
}
