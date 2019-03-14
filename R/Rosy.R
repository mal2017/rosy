#' Run rose stitching + quant for a single sample.
#'
#' Will take only a single set of regions and reads.
#'
#' @param regions GRanges.
#' @param reads Object or path to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
rose_single_internal <- function(regions, reads, paired = F,
                 colData = NULL, stitchDist = 12500) {

  if (length(list(regions)) > 1 | length(list(reads)) > 1) {
    stop("Provide only a single set of reads and regions to this function.")
  }

  if (is.null(colData)) {
    stop("Please provide a dataframe to the coldata argument.")
  }

  stopifnot(nrow(colData) == 1)

  st <- stitch(unlist(regions), stitchDist = stitchDist)

  rse <- liquidate_internal(features = st,
                            reads = unlist(reads),
                            paired = paired)

  rse <- normalize_mrip_internal(rse)

  rse <- normalize_perbp_internal(rse)

  SummarizedExperiment::colData(rse) <- S4Vectors::DataFrame(colData)

  rownames(rse) <- as.character(SummarizedExperiment::rowRanges(rse))

  rse
}

# TODO add a multiple rose wrapper
