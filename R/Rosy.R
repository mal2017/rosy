#' Run rosy on a single sample.
#' @param regions GRanges.
#' @param reads Object or path to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
rosy_internal <- function(regions, reads, paired = F,
                 colData = NULL, txdb = NULL, stitchDist = "optimize") {

  if (any(lapply(list(regions, reads, paired, colData), length) > 1)) {
    stop("Supply only a single sample to this function.")
  }

  # TODO: 0. get optimum stitching distance if called for

  st <- stitch(regions, stitchDist = 12500)

  # TODO parallelize this
  rse <- liquidate_internal(features = st,
                            reads = reads,
                            paired = paired)

  # TODO: 1. call supers
  # TODO: 2. get distances and constituent distances



  rse
}
