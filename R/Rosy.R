#' Run rosy on a single sample.
#' @param regions GRanges.
#' @param reads Object or path to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
rosy_internal <- function(regions, reads, paired = F,
                 colData = NULL, txdb = NULL, stitchDist = 12500) {

  if (any(lapply(list(regions, paired, colData), length) > 1)) {
    stop(paste("Only the `reads` arg allows for vectors and lists.",
         "Make sure single objects are supplied for other args."))
  }

  if (is.vector(reads) & length(reads) > 1 & is.null(names(reads))) {
    stop("You must name the items in your `reads` vector.")
  } else if (is.vector(reads) & length(reads) > 1 & !is.null(colData)) {
    if (!all(names(reads) %in% rownames(colData))) {
      stop("Not all of your named samples are present in `colData`!")
    }
    colData <- colData[names(reads),]
  }

  # TODO: 0. get optimum stitching distance if called for
  st <- stitch(regions, stitchDist = stitchDist)

  # TODO parallelize this
  rse <- liquidate_internal(features = st,
                            reads = reads,
                            paired = paired)

  # MRIP normalization
  rse <- normalize_internal(rse)

  # add coldata
  SummarizedExperiment::colData(rse) <- S4Vectors::DataFrame(colData)

  # TODO: 1. call supers
  # TODO: 2. get distances and constituent distances

  rse
}
