% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Liquidating.R
\name{liquidate_internal}
\alias{liquidate_internal}
\title{Find signal from a bam within regions.}
\usage{
liquidate_internal(features, reads, paired = F)
}
\arguments{
\item{reads}{A bam file path, GRanges list, or other object holding reads.}

\item{gr}{A GRList object.}
}
\value{
RangedSummarizedExperiment
}
\description{
Generic behavior is handled by GenomicAlignments. This is a thin wrapper.
}
