#' Run rose stitching + quant for a single sample.
#'
#' Will take only a single set of regions and reads.
#'
#' @param regions GRanges or path to bed.
#' @param reads Object or path to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
rose_single_internal <- function(regions, reads, paired = F,
                 colData = NULL, stitchDist = 12500, se_cutoff = "optimize") {

  if (is.list(regions) | is.list(reads) | is.vector(regions)) {
    stop("Provide only a single set of reads and regions to this function.")
  }

  if (is.null(colData)) {
    stop("Please provide a dataframe to the coldata argument.")
  }

  stopifnot(nrow(colData) == 1)

  st <- stitch(regions, stitchDist = stitchDist)

  rse <- liquidate_internal(features = st,
                            reads = unlist(reads),
                            paired = paired)

  rse <- normalize_mrip_internal(rse)

  rse <- normalize_perbp_internal(rse)

  SummarizedExperiment::colData(rse) <- S4Vectors::DataFrame(colData)

  rownames(rse) <- seGrlNames(SummarizedExperiment::rowRanges(rse))

  rse$stitchDist <- stitchDist

  call_supers_internal(rse, set_rnk = ifelse(se_cutoff == "optimize",F, se_cutoff))
}

