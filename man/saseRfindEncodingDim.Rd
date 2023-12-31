% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saseRfindEncodingDim.R
\name{saseRfindEncodingDim}
\alias{saseRfindEncodingDim}
\title{Determine the optimal number of latent factors to detect outlier gene
expression or splicing}
\usage{
saseRfindEncodingDim(
  se,
  method = "GD",
  analysis,
  dimensions = seq(2, min(100, ncol(se) - 2, nrow(se) - 1), 2),
  freq = 0.01,
  zScore = 3,
  sdlog = log(1.6),
  lnorm = TRUE,
  inj = "both",
  BPPARAM = bpparam(),
  aggregation,
  scale = TRUE,
  ...
)
}
\arguments{
\item{se}{A `SummarizedExperiment` instance generated with the
SummarizedExperiment function of the SummarizedExperiment package.
In the assay slot, provide the expression counts as an
ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
`colData` is a `DataFrame` describing the samples in the experiment. Also,
include `offsets` in the assays slot, which can be calculated with the
`calculateOffsets` function.
Finally, specify the experimental `design` as a formula in the metadata
slot. This formula must be based on the colData, and should be `~1` if no
known covariates are included.}

\item{method}{The method used to estimate the optimal number of latent
factors included in the regression framework. Default is `GD`, which uses
the Gavish and Donoho threshold for singular values. `DAE` will use a
Denoising autoencoder, but will require longer computation time without
increased performance.}

\item{analysis}{`AE` for aberrant expression analysis and `AS` for aberrant
splicing analysis. Used to insert corrupted counts when using the `DAE`
method.}

\item{dimensions}{A vector containing the number of latent factors that are
compared in the grid search for the optimal number of latent factors. Only
used when using the `DAE` method.}

\item{freq}{the frequency of corrupted counts injected in the dateset. Only
used when using the `DAE` method.}

\item{zScore}{The mean magnitude of the corrupted counts in the dataset.
Only used when using the `DAE` method.}

\item{sdlog}{Standard deviation of the distribution of the corrupted counts
on the log-scale. Only used when using the `DAE` method and when lnorm is
TRUE.}

\item{lnorm}{If TRUE, the corrupted counts are simulated from a log-normal
distribution with a mean of log(zScore) and standard deviation of sdlog.
Only used when using the `DAE` method.}

\item{inj}{Injection of overexpression ('high'), underexpression ('low') or
both ('both') types of corrupted counts. Only used when using the `DAE`
method.}

\item{BPPARAM}{object of class \code{bpparamClass} that specifies the
back-end to be used for computations. See
\code{bpparam} in \code{BiocParallel} package for details.}

\item{aggregation}{character vector representing the column in the rowData
to be used to calculate offsets when injecting corrupted counts according to
aberrant splicing. Only used when method is `DAE` and when analysis is `AE`.}

\item{scale}{boolean. If TRUE, the deviance residuals upon will be scaled to
mean 0 and sd = 1 before estimating the latent factors.}

\item{...}{Extra arguments for .fitRUV.}
}
\value{
An updated 'SummarizedExperiment' instance, now including:
`optimalEncDim` in the metadata slot, representing the estimated optimal
number of latent factors. `LatentFactorControl` in the metadata slot,
which represents the method used to estimate the optimal number of latent
factors (`GD` for Gavish and Donoho threshold, `DAE` for denoising
autoencoder). `deviances` in the assays slot when using the `GD` method,
representing the deviance residuals calculated using edgeR with an intercept
and known covariates. `svd_u`, `svd_d` and `svd_v` matrices in the metadata
slot when using the `GD` method, which represent the singular value
decomposition of the deviance residuals. `encDimTable` a data.table in the
metadata slot which represents the area under the curve to search for
corrupted counts at the different dimensions when using the `DAE` method.
}
\description{
saseRfindEncodingDim is used to search for the optimal number of latent
factors included in the regression to search for aberrant expression or
splicing. This optimal number of latent factors can be estimated by using
the Gavish and Donoho threshold for singular values, or by using a
denoising autoencoder, as described in our corresponding paper
}
\examples{

data(saseRExample, package = "saseR")

SEgenes <- saseRfindEncodingDim(SEgenes, method = "GD")
SEbins <- saseRfindEncodingDim(SEbins, method = "GD")
SEjunctions <- saseRfindEncodingDim(SEjunctions, method = "GD")


}
