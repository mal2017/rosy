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
  # don't want to think abt strands
  GenomicRanges::strand(gr) <- "*"

  # not beyound actual chromosome bounds
  gr <- GenomicRanges::reduce(gr)
  gr <- GenomicRanges::trim(gr)

  # expand, remove overlaps, trim again
  gr_expanded <- extend_merge_internal(gr = gr, width = sd)
  gr_expanded <- GenomicRanges::reduce(gr_expanded)
  gr_expanded <- GenomicRanges::trim(gr_expanded)

  # get hits
  hits <- GenomicRanges::findOverlaps(gr, gr_expanded)

  # init new col with empty
  mcols(gr)$stitched_region <- NA

  # to those cols with any overlaps, add subject hits
  mcols(gr)[queryHits(hits),"stitched_region"] <- S4Vectors::subjectHits(hits)

  # split into grl by the subject hits
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

get_total_constituent_width_internal <- function(grl) {
  sum(width(grl))
}

#' Get the total with of a GRL element.
get_grl_entry_width_internal <- function(x) {
  max(c(start(x), end(x))) - min(c(start(x), end(x)))
}

#' Get all GRL entry widths.
get_stitch_region_widths_internal <- function(grl) {
  vapply(grl, get_grl_entry_width_internal,FUN.VALUE = c(1))
}

get_total_stitch_region_width_internal <- function(grl) {

}

named_container <- function(vals,nms) {
  stopifnot(length(vals) == length(nms))
  names(vals) <- nms
  vals
}

#' Get the stitch stats for a range of distances so rosy can
#' internally choose an optimal distance.
get_stitch_stats_internal <- function(gr) {
  tbl <- tibble::tibble(STEP = seq(0, 2000, by = 500))
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


  tot_constit_width <- get_total_constituent_width_internal(gr)


  # for each possible stitching configuration:
  tbl <- dplyr::mutate(tbl, GRL = purrr::map(STEP, .f = function(x) stitch_internal(gr, sd = x)))

  # a list of vectors of widths
  stitch_region_widths <- named_container(tbl$GRL, tbl$STEP)

  tbl <- dplyr::mutate(tbl, NUM_REGIONS = purrr::map_int(GRL,.f=length))
  tbl <- dplyr::mutate(tbl, TOTAL_CONSTIT_WIDTH = tot_constit_width)
  #tbl <- dplyr::mutate(tbl,TOTAL_STITCH_REGION_WIDTH = purrr::map(STEP, .f = function(x) sum(stitch_region_widths[[as.character(x)]])))


  tbl
}
