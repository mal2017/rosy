#' Stitch ranges within `stitchDist` of each other.
#' @rdname stitch
#' @param object GRanges, GRangesList, or paths to files.
#' @export
setGeneric("stitch", function(object, stitchDist=12500) standardGeneric("stitch"))

#' Run rosy.
#' @rdname rosy
#' @param regions GRanges, GRangesList, or path to file.
#' @param reads Object or path to file containing reads.
#' @param paired Logical indicating paired-end reads.
#' @param colData dataframe with sample information.
#' @export
setGeneric("rosy", function(regions, reads, paired,
                            colData, txdb, stitchDist=12500) standardGeneric("rosy"))



#' @rdname stitch
setMethod("stitch", signature(object = "GRanges", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            stitch_internal(gr = object, sd = stitchDist)
          })

#' @rdname stitch
setMethod("stitch", signature(object = "GRangesList", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            object <- GenomicRanges::reduce(unlist(object))
            stitch_internal(gr = object, sd = stitchDist)
          })

#' @rdname stitch
setMethod("stitch", signature(object = "character", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            grs <- lapply(object,rtracklayer::import)
            grs <- GenomicRanges::GRangesList(grs)
            stitch(object = grs, stitchDist = stitchDist)
          })

#' @rdname rosy
setMethod("rosy", signature(regions = "character", reads = "character",
                            paired = "logical", colData = "data.frame",
                            txdb = "TxDb",
                            stitchDist = "numeric"),
          function(regions, reads, paired,
                   colData, txdb, stitchDist=12500) {
            rosy_internal(regions = regions,
                          reads = reads,
                          paired = paired,
                          txdb = txdb,
                          colData = colData,
                          stitchDist = stitchDist)
          })
