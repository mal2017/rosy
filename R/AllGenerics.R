#' Stitch ranges within `stitchDist` of each other.
#' @rdname stitch
#' @param object GRanges, GRangesList, or paths to files.
#' @export
setGeneric("stitch", function(object, stitchDist=12500) standardGeneric("stitch"))

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
