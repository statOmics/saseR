% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertASpli.R
\name{convertASpli}
\alias{convertASpli}
\title{Converting SummarizedExperiment with different ASpli-counts to
gene, bin or junction level counts.}
\usage{
convertASpli(ASpliSE, type = "none", filter = TRUE, ...)
}
\arguments{
\item{ASpliSE}{An SummarizedExperiment object obtained by the
`BamtoAspliCounts` function, which contains gene, bin and junction level
counts in the metadata slot.}

\item{type}{character vector specifying the counts to be extracted. Can be
`gene`, `bin` and `junction`.}

\item{filter}{logical value specifying to filter bin counts based on the
type of bin. Regions that contain both intronic and exonic regions are
filtered as the intron and exon bins are already defined in other features.}

\item{...}{extra parameters for filtering.}
}
\value{
A `SummarizedExperiment` instance, representing the `gene`, `bin` or
`junction` counts as specified.
}
\description{
convertASpli is used to obtain gene, bin or junction level counts from
a SummarizedExperiment object containing these three in the metadata slot.
}
\examples{
data(saseRExample, package = "saseR")

SEgenes <- convertASpli(ASpliSE, type = "gene")
SEbins <- convertASpli(ASpliSE, type = "bin")
SEjunctions <- convertASpli(ASpliSE, type = "junction")

}
