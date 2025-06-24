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
#' compromising upon power. `edgeR` can be used to perform parameter estimation
#' , although not scaling well with many covariates or latent factors.
#'
#' @param scale logical value to scale the deviance residuals before
#' performing the singular value decomposition.
#'
#' @param robustPCA This argument is in beta phase and should be used carefully.
#' Logical value to perform robust PCA using the `PcaHubert`
#' function, which decreases the influence of outlier values on the estimated
#' latent factors.
#' @param ignore_samples This argument is in beta phase and should be used carefully.
#' Character vector with names of samples that are not
#' used when performing the latent factor estimation and estimation of the
#' regression coefficients.

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
#' @import S4Vectors
#' @import MASS
#' @import SummarizedExperiment
#' @import PRROC
#' @import BiocGenerics
#' @import methods
#' @import GenomicRanges
#' @import DESeq2
#' @import GenomicFeatures
#' @import IRanges
#' @import MatrixGenerics
#' @importFrom dplyr `%>%` group_by summarise across mutate
#' @importFrom limma lmFit strsplit2
#' @importFrom data.table data.table .N
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats model.matrix p.adjust pnbinom pnorm qnbinom rlnorm
#' rmultinom runif
#'
#'
#' @examples
#' data(saseRExample, package = "saseR")
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
#' @export


saseRfit <- function(se,
                     analysis,
                     dimensions = NULL,
                     padjust = "BH",
                     BPPARAM = bpparam(),
                     fit = "fast",
                     scale = TRUE,
                     robustPCA = FALSE,
                     ignore_samples = NULL){


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
        if(!is.null(ignore_samples)){
            DGE <- .fitEdgeRDisp(
                se = se[,-which(rownames(colData(se)) %in% ignore_samples)],
                design = getDesign(se))

        } else {
        DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
        }
        fit_DGE <- glmFit(y = DGE)

        # Principal component analysis on the deviance residuals
        deviances <- .CalculateDeviances(se = se,
                                         fit_DGE = fit_DGE)
        if (scale == TRUE) {
            deviances <- (deviances - rowMeans(deviances))/rowSds(deviances)
            deviances[is.nan(deviances)] <- 0
        }

        assay(se,'deviances', withDimnames = FALSE) <- deviances


        if(!is.null(ignore_samples)){
            index_samples <- colnames(se)
            se_removed <- se[,ignore_samples]
            rowData(se_removed)$theta <- NULL
            se <- se[,-which(rownames(colData(se)) %in% ignore_samples)]
        }

        if (robustPCA == TRUE){
            if (!requireNamespace("rrcov", quietly = TRUE)) {
                stop("The 'rrcov' package is required for the robustPCA 
                    function. Please install it with 
                     install.packages('rrcov').")
            }
            SingularVectors <- rrcov::PcaHubert(t(deviances),
                                         k = dimensions,
                                         kmax = ncol(se),
                                         scale = FALSE)

            metadata(se)[['svd_u']] <- SingularVectors$scores
            rownames(metadata(se)[['svd_u']]) <- colnames(se)
            metadata(se)[['svd_d']] <- SingularVectors$eigenvalues
            metadata(se)[['svd_v']] <- SingularVectors$loadings

        } else {
            SingularVectors <- svd(t(assays(se)$deviances))

            metadata(se)[['svd_u']] <- SingularVectors$u
            rownames(metadata(se)[['svd_u']]) <- colnames(se)
            metadata(se)[['svd_d']] <- SingularVectors$d
            metadata(se)[['svd_v']] <- SingularVectors$v
        }
    } else {
        if(!is.null(ignore_samples)){
            index_samples <- colnames(se)
            se_removed <- se[,ignore_samples]
            se <- se[,-which(rownames(colData(se)) %in% ignore_samples)]
        }
    }


    dimensions <- metadata(se)[['optimalEncDim']]

    se <- .full_design(se=se)

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



    se <- .calcPValues(se, mu = assays(se)$mu, theta = rowData(se)$theta)
    if(!is.null(ignore_samples)){
        se_removed <- .calcPValuesIgnoredSamples(se_removed, se, dimensions)
        se <- cbind(se,se_removed)[,index_samples]
        metadata(se)$svd_u <- rbind(metadata(se)$svd_u,
                                    metadata(se)$svd_u_removed)[index_samples,]
        metadata(se)$svd_u_removed <- NULL
    }

    if(padjust != "none"){
        se <- .adjustPValues(se, method = padjust)
    }

    if(analysis %in% c("AS","proportion")){
        se <- .locusPValues(se)
    }

    return(se)



}


.full_design <- function(se){
    if(metadata(se)[['optimalEncDim']] == 0){
    design <- model.matrix(getDesign(se),
                               data = colData(se))
    } else if (getDesign(se) != ~1){
        env <- environment()
        formula <- update(getDesign(se), ~ . + metadata(se)[['svd_u']][,
                                seq_len(metadata(se)[['optimalEncDim']])],
                                env = env)
        environment(formula) <- env
        design <- model.matrix(formula, data = colData(se))
    } else {
        design <- model.matrix(~1+metadata(se)[['svd_u']][,
                        seq_len(metadata(se)[['optimalEncDim']])],
                        data = colData(se))
    }

    metadata(se)[["full_design"]] <- design
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
    number_of_known_confounders <- dim(model.matrix(getDesign(se),
                                                    data = colData(se)))[2]
    total_dimensions <- dimensions + number_of_known_confounders

    final_design <- metadata(se)[["full_design"]][,seq_len(total_dimensions)]
    beta_initial_reduced_dimensions <- beta_initial[,seq_len(total_dimensions)]

    # Form a design matrix with the known covariates and latent factors.
    # Then, an edgeR regression is performed with these known covariates
    # and latent factors.


    se <- .adapted_IRLS(se = se,
                        analysis = analysis,
                        initial_beta=beta_initial_reduced_dimensions,
                        final_design =
                            final_design)



    DGE <- .fitEdgeRDisp(se=se, offsets=assays(se)$mu, design=~1)
    theta <- 1/DGE$tagwise.dispersion
    rowData(se)$theta <- theta


    } else if (fit == "edgeR"){
        if (dimensions == 0){
            design <- model.matrix(getDesign(se), data = colData(se))
        }
        else if (getDesign(se) != ~1){
            env <- environment()
            formula <- update(getDesign(se), ~ . +
                                 metadata(se)[['svd_u']][,seq_len(dimensions)],
                              env = env)
            environment(formula) <- env
            design <- model.matrix(formula, data = colData(se))
        } else {
            design <- model.matrix(~ 1 +
                                metadata(se)[['svd_u']][,seq_len(dimensions)],
                                data = colData(se))
        }
        DGE <- .fitEdgeRDisp(se=se, design = design)
        fit_DGE <- glmFit(y = DGE)
        assay(se, 'mu', withDimnames=FALSE) <- fit_DGE$fitted.values
        rowData(se)$theta <- 1/fit_DGE$dispersion
        metadata(se)$coefficients <- fit_DGE$coefficients

    }
    #Calculate 2-sided p-values for every count based on previous estimates.


    # Store the p-values, estimated means and dispersion to the
    # OutriderDataSet.

    return(se)


}


.CalculateDeviances <- function(se, fit_DGE){

    # Returns the deviance residuals
    mu <- assays(se)$offsets * exp(fit_DGE$coefficients %*%
                                        t(model.matrix(getDesign(se),
                                        data = colData(se))))
    deviance <- nbinomUnitDeviance(counts(se),
                                   mean = mu,
                                   dispersion = fit_DGE$dispersion)
    deviance_scaled <- sqrt(abs(deviance))*sign(counts(se)-mu)

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
        if(!is.matrix(design)){
          design <- model.matrix(design, data = colData(se))
        }

    } else{
        design <- model.matrix(~1, data = colData(se))
    }

    DGE <- estimateDisp(DGE, design = design)
    return(DGE)
}

.initial_beta <- function(se){

    log_data <- log((counts(se)+1)/assays(se)$offsets)
    beta_initial<- lmFit(log_data,
                         design = metadata(se)$full_design)$coefficients

    return(beta_initial)
}



.adapted_IRLS <- function(se,
                          analysis,
                          final_design,
                          initial_beta,
                          maxiter=50){

    if(analysis == "AE"){
        max <- matrix(.Machine$double.xmax,ncol = ncol(se),nrow = nrow(se))

    } else if(analysis %in% c("AS","proportion")){
        max <- as.matrix(assays(se)$offsets)
    }
    design <- final_design
    inverse <- solve(t(design)%*%design)
    mu <- exp(t(design%*%t(initial_beta)))*assays(se)$offsets

    for(i in seq_len(maxiter)){
        z <- log(mu/assays(se)$offsets) + counts(se)/mu - 1

        beta <- inverse %*% t(design)%*%t(z)
        mu <- pmin(exp(t(design%*%beta))*assays(se)$offsets*0.2 + 0.8* mu,
                   max)
    }
    assay(se, 'mu', withDimnames=FALSE) <- mu
    metadata(se)$coefficients <- beta

    return(se)
}


.calcPValues <- function(se, mu, theta){
    left_quantile <- pnbinom(q = counts(se), mu = mu, size = theta)
    right_quantile <- pnbinom(q = counts(se)-1, mu = mu, size = theta)
    PValues <- pmin(left_quantile*2, (1-right_quantile)*2,1)
    assay(se, 'pValue', withDimnames=FALSE) <- PValues
    return(se)
}

.calcPValuesIgnoredSamples <- function(se_removed, se, dimensions){
    svd_u_removed <- t(assays(se_removed)$deviances) %*% metadata(se)$svd_v %*%
        diag(1/metadata(se)$svd_d)
    rownames(svd_u_removed) <- colnames(se_removed)

    metadata(se_removed)$svd_u_removed <- svd_u_removed

    mu_removed <- assays(se_removed)$offsets *
        exp(t(cbind(1,svd_u_removed[, seq_len(dimensions)]) %*%
                  metadata(se)$coefficients))
    assays(se_removed)$mu <- mu_removed

    se_removed <- .calcPValues(se_removed, mu_removed, rowData(se)$theta)

    return(se_removed)

}

.locusPValues <- function(se){
    pvalues_before_group <- data.frame(assays(se)$pValue,
                                       "locus" = rowData(se)$locus)

    pvalues_grouped <- pvalues_before_group %>% group_by(locus) %>%
        summarise(across(colnames(se),function(i) min(i)))

    rownames <- pvalues_grouped$locus
    pvalues_grouped[,1] <- NULL
    pvalues_grouped <- pvalues_grouped %>% as.matrix()
    rownames(pvalues_grouped) <- rownames
    metadata(se)$pValuesLocus <- pvalues_grouped
    return(se)

}

.adjustPValues <- function(se, method = "BH"){
    assay(se, 'pValueAdjust', withDimnames=FALSE) <-
        apply(X = assay(se, 'pValue', withDimnames=FALSE),
              MARGIN = 2,
              FUN = p.adjust,
              method = method)
    metadata(se)[['pValueAdjustMethod']] <- method
    return(se)
}





