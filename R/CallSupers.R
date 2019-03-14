call_supers_internal <- function(rse, set_rnk =F, assay="readInPeakNormPerBP") {

  mat <- matrix(nrow = nrow(rse), ncol = ncol(rse))
  colnames(mat) <- colnames(rse)
  rownames(mat) <- rownames(rse)

  for (s in colnames(rse)) {
    signals <- SummarizedExperiment::assay(rse[,s], i = assay)[,1]
    if (set_rnk) {
      cutoff <- signals[order(signals, decreasing = T)[set_rnk]]
    } else {
      cutoff <- call_super_cutoff_internal(signals)
    }
    mat[,s] <- signals >= cutoff
  }
  SummarizedExperiment::assay(rse, i = "isSuper") <- mat
  rse
}


#' Call 'supers' on a vector of values.
#' @importFrom magrittr %>%
call_super_cutoff_internal <- function(vec, target_slope = (max(vec) - min(vec)) / length(vec)) {
  z <- tibble::tibble(signal=vec) %>%
    dplyr::mutate(rank = dplyr::dense_rank(signal)) %>%
    dplyr::arrange(rank) %>%
    dplyr::mutate(slope = NA)

  for (i in 1:(nrow(z) - 1)) {
    z[[i, "slope"]] <- z[[i + 1, "signal"]] - z[[i, "signal"]]
  }

  # all points ranked lower than this are not super
  cutoff_rnk <- floor(optimize(pts_below,
                 lower = 1,
                 upper =  nrow(z),
                 v = z$signal,
                 a = target_slope)[["minimum"]])

  z[[cutoff_rnk,"signal"]]
}

#' optimize for the value of this function when calling supers
pts_below <- function(x,v,a) {
  y <- v[x] # y pt
  b <- y - a*x  # intercept
  xs <- 1:length(v) # x pts
  theoretical_y <- (a*xs) + b
  sum(v < theoretical_y)
}
