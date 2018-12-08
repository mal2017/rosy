#' Run rosy.
#' @rdname stitch
#' @param regions GRanges, GRangesList, or paths to files.
#' @param reads Objects or paths to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
#' @export
rosy <- function(regions, reads, paired = F, colData = None) {
  st <- stitch(regions, stitchDist = 12500)
  args_for_cts <- mapply(list, reads, paired, SIMPLIFY = F)
  # TODO parallelize this
  rses <- lapply(args_for_cts,
                 FUN = function(x) {
                   liquidate_internal(features = st,
                                      reads = x[[1]],
                                      paired = x[[2]],
                                      normalize = T)
                 })
  rses
}
