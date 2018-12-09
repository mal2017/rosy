#' Stitch GRanges within specified distance
#'
#' Returns a GRanges object stitched together such that
#' ranges within the specified stitch distance of each other
#' are now reduced. In the process, metadata  and
#' strand information is lost.
#'
#' @param gr A GRanges object.
#' @param sd Stitching distance.
stitch_internal <- function(gr, sd = 0) {
  gr <- GenomicRanges::trim(gr)
  gr_expanded <- extend_merge_internal(gr = gr, width = sd)
  GenomicRanges::strand(gr_expanded) <- "*"
  GenomicRanges::strand(gr) <- "*"

  hits <- GenomicRanges::findOverlaps(gr, gr_expanded)

  gr$stitched_region <- S4Vectors::subjectHits(hits)

  GenomicRanges::split(gr, as.factor(gr$stitched_region))
}


#' Extend GRanges by the specified distance
#'
#' A wrapper around `GenomicRanges::flank`.
#'
#' Returns a GRanges object with flanking regions added to every distinct
#' region. Metadata and strand info is lost in this process.
#'
#' @param gr A GRanges object.
#' @param sd Stitching distance.
extend_merge_internal <- function(gr, width = 0) {
  GenomicRanges::strand(gr) <- "*"
  gr_start <- c(gr, GenomicRanges::flank(gr, width = width, start = T, both = T))
  gr_start <- GenomicRanges::reduce(gr_start)
  gr_end <- c(gr,GenomicRanges::flank(gr, width = width, start=F, both = T))
  gr_end <- GenomicRanges::reduce(gr_end)
  gr2 <- GenomicRanges::reduce(c(gr_start, gr_end))
  GenomicRanges::trim(gr2)
}

get_total_constituent_width_internal <- function(gr) {
  sum(width(grl))
}

get_stitch_stats_internal <- function(gr) {
  tbl <- tibble::tibble(STEP = seq(0, 1000, by = 500))
                 #GRL = NULL,
                 #NUM_REGIONS = NULL,
                 #TOTAL_CONSTIT_WIDTH  = get_total_constituent_width_internal(gr),
                 #TOTAL_STITCH_REGION_WIDTH = NULL,
                 #MEAN_CONSTIT_WIDTH_PER_STITCH_REGION = NULL,
                 #MED_CONSTIT_WIDTH_PER_STITCH_REGION = NULL,
                 #MEAN_STITCH_REGION_WIDTH = NULL,
                 #MED_STITCH_REGION_WIDTH = NULL,
                 #STITCH_FRACTIONS = NULL,
                 #MEAN_STITCH_FRACTION =  NULL,
                 #MED_STITCH_FRACTION = NULL)

  tbl <- dplyr::mutate(tbl, GRL = purrr::map(STEP, .f = function(x) stitch_internal(gr, sd = x)))
  tbl <- dplyr::mutate(tbl, NUM_REGIONS = purrr::map_int(GRL,.f=length))
  tbl <- dplyr::mutate(tbl, TOTAL_CONSTIT_WIDTH = get_total_constituent_width_internal(gr))
  tbl
}
