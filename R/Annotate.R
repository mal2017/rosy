#' Call the closest gene for a gr and a txdb
closest_gene_internal <- function(gr, txdb) {
  proms <- GenomicFeatures::promoters(txdb,
                                      upstream = 1000,
                                      downstream=1000,
                                      use.names = T)
  genes <- proms[GenomicRanges::nearest(gr,proms)]

  return(names(genes))
}

#' rdname sds2tbl
#'
#' Convert an SE dataset to a tibble for plotting and downstream
#' querying.
#'
#' @importFrom magrittr %>%
#' @export
sds2tbl <- function(sds) {
  # TODO check validity as SDS

  tl <- lapply(colnames(sds), FUN=function(x) sds[,x]) %>%
        lapply(sds2tbl_single_internal)

  dplyr::bind_rows(tl)
}


#' Convert an SE dataset to a tibble for plotting and downstream
#' querying.
#'
#' @importFrom magrittr %>%
sds2tbl_single_internal <- function(sds) {

  stopifnot(ncol(sds)==1)

  rngs <- SummarizedExperiment::rowRanges(sds)

  result <- tibble::tibble(locus = rownames(sds),
                           sample = colnames(sds),
                           isSuper = SummarizedExperiment::assay(sds,i = "isSuper")[,1],
                           rawCount = SummarizedExperiment::assay(sds,i = "counts")[,1],
                           readInPeakNormPerBP = SummarizedExperiment::assay(sds,i = "readInPeakNormPerBP")[,1],
                           nConstitPeaks = vapply(rngs, length, FUN.VALUE = c(88)),
                           constitPeaks = vapply(rngs,FUN=function(x) paste(as.character(x),collapse = ","),FUN.VALUE=c("")))

  result <- SummarizedExperiment::colData(sds) %>%
    dplyr::as_tibble(rownames = "sample") %>%
    dplyr::left_join(result, ., by="sample")


  result <- SummarizedExperiment::rowData(sds) %>%
    dplyr::as_tibble() %>%
    dplyr::bind_cols(result,.)

  # this must be last because I rearrange
  result <- result %>%
    dplyr::arrange(desc(readInPeakNormPerBP)) %>%
    dplyr::mutate(rankInSample = dplyr::row_number())

  return(result)
}
