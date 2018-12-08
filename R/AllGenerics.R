#' @rdname stitch
#' @export
setGeneric("stitch", function(object, stitchDist=12500) standardGeneric("stitch"))

#' Stitch ranges within stitchDist.
#' @rdname stitch
#' @param object GRanges object
#' @export
setMethod("stitch", signature(object = "GRanges", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            stitch_internal(gr = object, sd = stitchDist)
          })

