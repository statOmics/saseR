#' @title Searching for aberrant expression or splicing.

#'

#' @description

#'

#' saseRfit is used to perform a negative binomial regression and to search
#' for aberrant expression or splicing in the context of rare Mendelian
#' disorders. To model aberrant splicing, it uses adapted offsets to model
#' proportions. It uses a predefined or an estimated optimal number of latent
#' factors in the regression. It returns the p-values for the prioritisation
#' of aberrant expression or splicing.

#' @param se A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' `colData` is a `DataFrame` describing the samples in the experiment.
#' Finally, specify the experimental design as a formula in the metadata slot.
#' This formula must be based on the colData, or equal to ~1 in case of
#' modelling only an intercept. If 'optimalEncDim' is available in the metadata
#' slot, this number will be used as the number of latent factors.
#'
#'
#' @param analysis `AE` for aberrant expression analysis, `AS` for aberrant
#'  splicing analysis, and `proportion` for any other kind of analysis using
#'  proportions.


#' @param dimensions A single integer representing the number of latent factors
#' that are used in the regression framework. If not specified, the
#' 'optimalEncDim' in the metadata slot will be used. If nor 'optimalEncDim' is
#' specified, the Gavish and Donoho threshold for singular values will be used.

#' @param padjust the method used to compute FDR-adjusted p-values. Look at the
#' documentation of the `p.adjust` package for different options. Both
#' Benjamini-Hochberg ('BH') or Benjamini-Yekutieli ('BY') are recommended.

#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.
#'
#' @param fit character value. Both `fast` or `edgeR` can be used.
#' `fast` uses a quadratic mean-variance relationship for fast estimation of
#' the mean when many covariates or latent factors are included, without
#' compromising upon power. `edgeR` can be used to perform parameter estimation,
#' although not scaling well with many covariates or latent factors.
#'
#' @param scale logical value to scale the deviance residuals before
#' performing the singular value decomposition.
#'
#' @param robustPCA logical value to perform robust PCA using the `PcaHubert`
#' function, which decreases the influence of outlier values on the estimated
#' latent factors.

#' @return An updated `SummarizedExperiment` instance, now including a matrix
#' of p-values (`pValue`) and a matrix of FDR-adjusted p-values
#' (`pValueAdjust`) in the assay slot. Also `pValueAdjustMethod` is available
#' in the metadata slot, stating which method was used to obtain FDR-adjusted
#' pvalues. If analysis is `AS` or `proportion`, also `pValuesLocus`, an
#' aggregated p-value is provided based on the `locus` column in the rowData.
#' This is calculated by taking the minimal p-values.
#'
#' @import ASpli
#' @import edgeR
#' @import MASS
#' @import pracma
#' @import SummarizedExperiment
#' @import precrec
#' @import PRROC
#' @import BiocGenerics
#' @import methods
#' @import havok
#' @import GenomicRanges
#' @import DESeq2
#' @importFrom rrcov PcaHubert
#' @importFrom limma lmFit strsplit2
#' @importFrom data.table data.table
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats aggregate median model.matrix p.adjust pnbinom  pnorm  qnbinom rlnorm rmultinom runif
#'
#'
#' @examples
#'
#' \dontrun{
#' gtfFileName <- aspliExampleGTF()
#' BAMFiles <- aspliExampleBamList()
#' targets <- data.frame(
#'     row.names = paste0('Sample',c(1:12)),
#'     bam = BAMFiles,
#'     f1 = rep("A",12),
#'     stringsAsFactors = FALSE)
#' genomeTxDb <- makeTxDbFromGFF(gtfFileName)
#' features <- binGenome(genomeTxDb)
#'
#' ASpliSE <- BamtoAspliCounts(
#'     features = features,
#'     targets = targets,
#'     minReadLength = 100,
#'     libType = "SE",
#'     BPPARAM = MulticoreParam(1L)
#' )
#'
#' SEgenes <- convertASpli(ASpliSE, type = "gene")
#' SEbins <- convertASpli(ASpliSE, type = "bin")
#' SEjunctions <- convertASpli(ASpliSE, type = "junction")
#'
#' metadata(SEgenes)$design <- ~1
#' metadata(SEbins)$design <- ~1
#' metadata(SEjunctions)$design <- ~1
#'
#' SEgenes <- calculateOffsets(SEgenes, method = "TMM")
#' SEbins <- calculateOffsets(SEbins,
#'                            method = "AS",
#'                            aggregation = "locus")
#' SEjunctions <- calculateOffsets(SEjunctions,
#'                                 method = "AS",
#'                                 aggregation = "symbol",
#'                                 mergeGeneASpli = TRUE)
#'
#' SEgenes <- saseRfindEncodingDim(SEgenes, method = "GD")
#' SEbins <- saseRfindEncodingDim(SEbins, method = "GD")
#' SEjunctions <- saseRfindEncodingDim(SEjunctions, method = "GD")
#'
#' SEgenes <- saseRfit(SEgenes,
#'                     analysis = "AE",
#'                     padjust = "BH",
#'                     fit = "fast")
#' SEbins <- saseRfit(SEbins,
#'                    analysis = "AS",
#'                    padjust = "BH",
#'                    fit = "fast")
#' SEjunctions <- saseRfit(SEjunctions,
#'                         analysis = "AS",
#'                         padjust = "BH",
#'                         fit = "fast")
#'}
#' @export


saseRfit <- function(se,
                     analysis,
                     dimensions = NULL,
                     padjust = "BH",
                     BPPARAM = bpparam(),
                     fit = "fast",
                     scale = TRUE,
                     robustPCA = FALSE){

    if(is.null(dimensions)){
        if(is.null(metadata(se)[['optimalEncDim']])){
            warning("No number of latent factors given, nor calculated with
                    'saseRfindEncodingDim'. Therefore, the Gavish and Donoho
                    threshold will be used to estimate the optimal number of
                    latent factors.")

            se <- .GavishDonoho(se)
        }

    } else{
        if(!is.null(metadata(se)[['optimalEncDim']])){
            warning("Manually imputed number of dimensions will overwrite the
                    calculated number. Set dimensions to 'NULL' if not
                    desired.")
        }
        metadata(se)[['oldoptimalEncDim']] <- metadata(se)[['optimalEncDim']]
        metadata(se)[['oldLatentFactorControl']] <-
            metadata(se)[['LatentFactorControl']]
        metadata(se)[['optimalEncDim']] <- dimensions
        metadata(se)[['LatentFactorControl']] <- "Manual"

    }

    if(!(metadata(se)[['LatentFactorControl']] %in% c("GavishDonoho"))){
        DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
        fit_DGE <- glmFit(y = DGE)

        # Principal component analysis on the deviance residuals
        deviances <- .CalculateDeviances(se = se,
                                         fit_DGE = fit_DGE)
        if (scale == TRUE) {
            deviances <- (deviances - rowMeans(deviances))/rowSds(deviances)

        }
        if (robustPCA == TRUE){
            SingularVectors <- PcaHubert(t(deviances),
                                         k = dimensions,
                                         kmax = ncol(se),
                                         scale = FALSE)

            metadata(se)[['svd_u']] <- SingularVectors$scores
            metadata(se)[['svd_d']] <- SingularVectors$eigenvalues
            metadata(se)[['svd_v']] <- SingularVectors$loadings

        } else {

            SingularVectors <- svd(t(deviances))

            assay(se,'deviances', withDimnames = FALSE) <- deviances

            metadata(se)[['svd_u']] <- SingularVectors$u
            metadata(se)[['svd_d']] <- SingularVectors$d
            metadata(se)[['svd_v']] <- SingularVectors$v
        }
    }
    dimensions <- metadata(se)[['optimalEncDim']]


    #Orthonormalisation
    se <- .orthonormalisation(se=se)

    # Calculate for each number of latent factors (= dimensions) that is
    # given as input the area under the curve of all merged p-values when
    # doing a RUV analysis with as true positives the corrupted counts.
    # This grid search is done to obtain the optimal number of latent factors
    # to include in the model.

    beta_initial <- .initial_beta(se=se)

    se <- .fitRUV(se = se,
                  analysis = analysis,
                  dimensions = dimensions,
                  beta_initial = beta_initial,
                  padjust = padjust,
                  fit = fit)

    se <- .adjustPValues(se, method = padjust)

    return(se)



}


.orthonormalisation <- function(se){
    if (getDesign(se) != ~1){
        env <- environment()
        formula <- update(getDesign(se), ~ . + metadata(se)[['svd_u']][,
                                1:metadata(se)[['optimalEncDim']]], env = env)
        environment(formula) <- env
        design <- model.matrix(formula, data = colData(se))
    } else {
        design <- model.matrix(~1+metadata(se)[['svd_u']][,
                        1:metadata(se)[['optimalEncDim']]], data = colData(se))
    }
    qr <- gramSchmidt(A = design)

    metadata(se)[["Q"]] <- qr$Q
    metadata(se)[["R"]] <- qr$R
    return(se)

}


.fitRUV <- function(se,
                    analysis,
                    dimensions,
                    beta_initial,
                    padjust,
                    fit){


    # Extract the number of latent factors that is estimated to be optimal
    # by prior knowledge or hyperparameter optimisation.
    if(fit == "fast"){
    number_of_known_confounders <- dim(model.matrix(getDesign(se), data = colData(se)))[2]
    total_dimensions <- dimensions + number_of_known_confounders

    reduced_orthonormal_vector_space <- metadata(se)[["Q"]][,1:total_dimensions]
    beta_initial_reduced_dimensions <- beta_initial[,1:total_dimensions]

    # Form a design matrix with the known covariates and latent factors.
    # Then, an edgeR regression is performed with these known covariates
    # and latent factors.


    mu <- .adapted_IRLS(se = se,
                        analysis = analysis,
                        initial_beta=beta_initial_reduced_dimensions,
                        reduced_orthonormal_vector_space =
                            reduced_orthonormal_vector_space)

    assay(se, 'mu', withDimnames=FALSE) <- mu


    DGE <- .fitEdgeRDisp(se=se, offsets=mu, design=~1)
    theta <- 1/DGE$tagwise.dispersion
    rowData(se)$theta <- theta

    } else if (fit == "edgeR"){
        if (dimensions == 0){
            design <- model.matrix(getDesign(se), data = colData(se))
        }
        else if (getDesign(se) != ~1){
            env <- environment()
            formula <- update(getDesign(se), ~ . + metadata(se)[['svd_u']][,1:dimensions], env = env)
            environment(formula) <- env
            design <- model.matrix(formula, data = colData(se))
        } else {
            design <- model.matrix(~1+metadata(se)[['svd_u']][,1:dimensions], data = colData(se))
        }
        DGE <- .fitEdgeRDisp(se=se, design = design)
        fit_DGE <- glmFit(y = DGE)
        assay(se, 'mu', withDimnames=FALSE) <- fit_DGE$fitted.values
        rowData(se)$theta <- 1/fit_DGE$tagwise.dispersion
    }
    #Calculate 2-sided p-values for every count based on previous estimates.
    se <- .calcPValues(se, mu = mu, theta = theta)


    if(padjust != "none"){
        se <- .adjustPValues(se, method = padjust)
    }

    if(analysis %in% c("AS","proportion")){
        se <- .locusPValues(se)
    }

    # Store the p-values, estimated means and dispersion to the
    # OutriderDataSet.

    return(se)


}


.CalculateDeviances <- function(se, fit_DGE){

    # Returns the deviance residuals

    deviance <- nbinomUnitDeviance(counts(se),
                                   mean = fit_DGE$fitted.values,
                                   dispersion = fit_DGE$dispersion)
    deviance_scaled <- sqrt(abs(deviance))*sign(counts(se)-fit_DGE$fitted.values)

    return(deviance_scaled)
}


.fitEdgeRDisp <- function(se, offsets = NULL, design=~1, ...){
    # Fit edgeR model on intercept and known covariates. This is
    # further used to calculate the deviances.

    DGE <- DGEList(counts=counts(se))
    if(is.null(offsets)){
        DGE$offset <- log(assays(se)$offsets)
    } else {
        DGE$offset <- log(offsets)
    }

    if(!(design == ~1)){
        design <- model.matrix(getDesign(se), data = colData(se))
    } else{
        design <- model.matrix(~1, data = colData(se))
    }

    DGE <- estimateDisp(DGE, design = design)
    return(DGE)
}

.initial_beta <- function(se){

    log_data <- log((counts(se)+1)/assays(se)$offsets)
    beta_initial<- lmFit(log_data, design = metadata(se)[["Q"]])$coefficients

    return(beta_initial)
}



.adapted_IRLS <- function(se, analysis, initial_beta, reduced_orthonormal_vector_space, maxiter=50){

    if(analysis == "AE"){
        max <- matrix(.Machine$double.xmax,ncol = ncol(se),nrow = nrow(se))

    } else if(analysis %in% c("AS","proportion")){
        max <- as.matrix(assays(se)$offsets)
    }


    mu <- exp(t(reduced_orthonormal_vector_space%*%t(initial_beta)))*assays(se)$offsets

    for(i in c(1:maxiter)){
        z <- log(mu/assays(se)$offsets) + counts(se)/mu - 1

        beta <- t(reduced_orthonormal_vector_space)%*%t(z)
        mu <- pmin(exp(t(reduced_orthonormal_vector_space%*%beta))*assays(se)$offsets*0.2 + 0.8* mu,
                   max)
    }

    return(mu)
}


.calcPValues <- function(se, mu, theta){
    left_quantile <- pnbinom(q = counts(se), mu = mu, size = theta)
    right_quantile <- pnbinom(q = counts(se)-1, mu = mu, size = theta)
    PValues <- pmin(left_quantile*2, (1-right_quantile)*2,1)
    assay(se, 'pValue', withDimnames=FALSE) <- PValues
    return(se)
}


.locusPValues <- function(se){
    pvalues_before_group <- data.frame(assays(se)$pValue,"locus" = rowData(se)$locus)

    pvalues_grouped <- pvalues_before_group %>% group_by(locus) %>%
        summarise(across(colnames(se),function(i) min(i)))

    rownames <- pvalues_grouped$locus
    pvalues_grouped[,1] <- NULL
    rownames(pvalues_grouped) <- rownames
    metadata(se)$pValuesLocus <- as.matrix(pvalues_grouped)
    return(se)

}

.adjustPValues <- function(se, method = "BH"){
    assay(se, 'pValueAdjust', withDimnames=FALSE) <-
        apply(X = assay(se, 'pValue', withDimnames=FALSE), MARGIN = 2, FUN = p.adjust, method = method)
    metadata(se)[['pValueAdjustMethod']] <- method
    return(se)
}





