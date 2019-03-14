#' Stitch ranges within `stitchDist` of each other.
#' @rdname stitch
#' @param object GRanges, GRangesList, or paths to files.
#' @export
setGeneric("stitch", function(object, stitchDist=12500) standardGeneric("stitch"))

#' Run rosy.
#' @rdname rosy
#' @param regions GRanges, list of GRanges, GRangesList, path to file, or list of paths.
#' @param reads Object or path to file containing reads, or list.
#' @param paired Logical indicating paired-end reads or list.
#' @param colData dataframe with sample information.
#' @param stitchDist Enhancer stitching distance.
#' @param se_cutoff Either 'optimize' or an integer indicating the top n enhancers to consider as SE.
#' @param exclude GRanges of regions to exclude
#' @export
setGeneric("rosy", function(regions, reads, paired=F,
                            colData, txdb=NULL, stitchDist=12500, se_cutoff="optimize", exclude=NULL) standardGeneric("rosy"))



#' @rdname stitch
setMethod("stitch", signature(object = "GRanges", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            GenomicRanges::strand(object) <- "*"
            stitch_internal(gr = object, sd = stitchDist)
          })

#' @rdname stitch
setMethod("stitch", signature(object = "GRangesList", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            GenomicRanges::strand(object) <- "*"
            object <- GenomicRanges::reduce(unlist(object))
            stitch_internal(gr = object, sd = stitchDist)
          })

#' @rdname stitch
setMethod("stitch", signature(object = "character", stitchDist = "numeric"),
          function(object, stitchDist = 12500) {
            grs <- lapply(object,rtracklayer::import)
            grs <- GenomicRanges::GRangesList(grs)
            GenomicRanges::strand(grs) <- "*"
            stitch(object = grs, stitchDist = stitchDist)
          })

#' @rdname rosy
setMethod("rosy", signature(regions = "GRanges", reads = "character",
                            colData = "data.frame"),
          function(regions, reads, paired=F,
                   colData, txdb=NULL, stitchDist=12500, se_cutoff="optimize", exclude=NULL) {
            stopifnot(length(list(regions)) == 1 & length(list(reads)) & ncol(colData == 1))
            if (!is.null(exclude)) {
              regions <- IRanges::subsetByOverlaps(regions,exclude, invert = T)
            }

            rse <- rose_single_internal(regions = regions,
                          reads = reads[1],
                          paired = paired[1],
                          colData = colData[1,],
                          stitchDist = stitchDist[1])
            if (!is.null(txdb)) {
              SummarizedExperiment::rowData(rse)[,"ClosestTx"] <- closest_gene_internal(GenomicRanges::GRanges(rownames(rse)), txdb)
            }
            return(rse)
          })
