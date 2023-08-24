#' @keywords internal
"_PACKAGE"

#' saseR package for aberrant expression and splicing analysis from bulk
#' RNA-seq data
#'
#' saseR is a highly performant and fast algorithm to detect aberrant expression and splicing.
#' The main functions are:
#'
#' \itemize{
#' \item \code{\link{BamtoAspliCounts}} - Process BAM files to ASpli counts
#' \item \code{\link{convertAspli}} - Get gene, bin or junction counts from ASpli SummarizedExperiment
#' \item \code{\link{calculateOffsets}} - Create an offsets assays for aberrant expression or splicing analysis
#' \item \code{\link{saseRfindEncodingDim}} - Estimate the optimal number of latent factors to include when estimating the mean expression
#' \item \code{\link{saseRfit}} - Parameter estimation of the negative binomial distribution and compute p-values for aberrant expression and splicing
#' }
#'
#' For information upon how to use these functions, check out our vignette at \url{https://github.com/statOmics/saseR/blob/main/vignettes/Vignette.Rmd}.
#'
#' @references
#'
#' Segers, A. et al. (2023).
#' Juggling offsets unlocks RNA-seq tools for fast scalable differential usage, aberrant splicing and expression analyses. bioRxiv.
#' \url{https://doi.org/10.1101/2023.06.29.547014}
#'
#' @author Alexandre Segers, Jeroen Gilis, Mattias Van Heetvelde, Elfride De Baere, Lieven Clement
#'
#' @docType package
#' @name saseR
## usethis namespace: start
## usethis namespace: end
NULL
