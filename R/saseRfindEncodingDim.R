#' @title Determine the optimal number of latent factors to detect outlier gene
#' expression or splicing

#'

#' @description

#'

#' saseRfindEncodingDim is used to search for the optimal number of latent

#' factors included in the regression to search for aberrant expression or

#' splicing. This optimal number of latent factors can be estimated by using
#' the Gavish and Donoho threshold for singular values, or by using a
#' denoising autoencoder, as described in our corresponding paper

#' @param se A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' `colData` is a `DataFrame` describing the samples in the experiment. Also,
#' include `offsets` in the assays slot, which can be calculated with the
#' `calculateOffsets` function.
#' Finally, specify the experimental `design` as a formula in the metadata slot.
#'  This formula must be based on the colData, and should be `~1` if no known
#'  covariates are included.
#'
#' @param method The method used to estimate the optimal number of latent
#'  factors included in the regression framework. Default is `GD`, which uses
#'  the Gavish and Donoho threshold for singular values. `DAE` will use a
#'  Denoising autoencoder, but will require longer computation time without
#'  increased performance.
#'
#' @param analysis `AE` for aberrant expression analysis and `AS` for aberrant
#'  splicing analysis. Used to insert corrupted counts when using the `DAE`
#'  method.

#' @param dimensions A vector containing the number of latent factors that are
#' compared in the grid search for the optimal number of latent factors. Only
#' used when using the `DAE` method.

#' @param freq the frequency of corrupted counts injected in the dateset. Only
#' used when using the `DAE` method.

#' @param zScore The mean magnitude of the corrupted counts in the dataset. Only
#' used when using the `DAE` method.

#' @param sdlog Standard deviation of the distribution of the corrupted counts
#' on the log-scale. Only used when using the `DAE` method and when lnorm is
#' TRUE.

#' @param lnorm If TRUE, the corrupted counts are simulated from a log-normal
#' distribution with a mean of log(zScore) and standard deviation of sdlog.
#' Only used when using the `DAE` method.


#' @param inj Injection of overexpression ('high'), underexpression ('low') or
#' both ('both') types of corrupted counts. Only used when using the `DAE`
#' method.


#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.

#' @param aggregation character vector representing the column in the rowData
#' to be used to calculate offsets when injecting corrupted counts according to
#' aberrant splicing. Only used when method is `DAE` and when analysis is `AE`.
#' @param ... Extra arguments for .fitRUV.


#' @return An updated 'SummarizedExperiment' instance, now including:
#' `optimalEncDim` in the metadata slot, representing the estimated optimal
#' number of latent factors. `LatentFactorControl` in the metadata slot,
#' which represents the method used to estimate the optimal number of latent
#' factors (`GD` for Gavish and Donoho threshold, `DAE` for denoising
#' autoencoder). `deviances` in the assays slot when using the `GD` method,
#' representing the deviance residuals calculated using edgeR with an intercept
#' and known covariates. `svd_u`, `svd_d` and `svd_v` matrices in the metadata
#' slot when using the `GD` method, which represent the singular value
#' decomposition of the deviance residuals. `encDimTable` a data.table in the
#' metadata slot which represents the area under the curve to search for
#' corrupted counts at the different dimensions when using the `DAE` method.
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
#'
#'}
#' @export

saseRfindEncodingDim <- function(se,
                                 method = "GD",
                                 analysis,
                               dimensions=seq(
                                   2, min(100, ncol(se) - 2, nrow(se) - 1), 2),
                               freq=1E-2,
                               zScore=3,
                               sdlog=log(1.6),
                               lnorm=TRUE,
                               inj='both',
                               BPPARAM=bpparam(),
                               aggregation,
                               ...) {

    if(!(method %in% c("GD","DAE"))){
        stop("Method to estimate the optimal number of latent factors has to
             be GD (GavishDonoho threshold) or DAE (denoising autoencoder). To
             manually impute the number of latent factors, include these when
             fitting the model.")
    }

    if(method == "GD"){
        .GavishDonoho(se)

    } else if(method == "DAE"){

        if(!(exists("analysis"))){
            stop("Please specify an analysis to make. Options are AE for aberrant
             expression and AS for aberrant splicing.")

        } else if(!(analysis %in% c("AE", "AS"))){
            stop("Only AE and AS in the analysis argument are possible to estimate
             the optimal number of latent factors with a denoising autoencoder")

        } else if(analysis == "AE"){
            .DAE_AE(se,
                    analysis,
                    dimensions,
                    freq,
                    zScore,
                    sdlog,
                    lnorm,
                    inj,
                    BPPARAM)
        } else if(analysis == "AS"){
            .DAE_AS(se,
                    analysis,
                    aggregation,
                    dimensions,
                    freq,
                    zScore,
                    sdlog,
                    lnorm,
                    inj,
                    BPPARAM)

        }

    }



}

.GavishDonoho <- function(se){
    DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
    fit_DGE <- glmFit(y = DGE)

    # Principal component analysis on the deviance residuals
    deviances <- .CalculateDeviances(se = se,
                                     fit_DGE = fit_DGE)
    SingularVectors <- svd(t(deviances))
    dimensions <- .OptimalHardThreshold(SingularVectors)

    metadata(se)[['optimalEncDim']] <- dimensions
    metadata(se)[['LatentFactorControl']] <- "GavishDonoho"

    assay(se,'deviances', withDimnames = FALSE) <- deviances

    metadata(se)[['svd_u']] <- SingularVectors$u
    metadata(se)[['svd_d']] <- SingularVectors$d
    metadata(se)[['svd_v']] <- SingularVectors$v

    return(se)

}



.DAE_AE <- function(se,
                    analysis,
                               dimensions=seq(
                                   2, min(100, ncol(se) - 2, nrow(se) - 1), 2),
                               freq=1E-2,
                               zScore=3,
                               sdlog=log(1.6),
                               lnorm=TRUE,
                               inj='both',
                               BPPARAM=bpparam(),
                               ...){

    assay(se,'trueOffsets', withDimnames = FALSE) <- assays(se)$offsets

    # Find the optimal number of hyperparameters when using RUV with
    # deviance residuals.
    # Estimate size factors

    se <- calculateOffsets(se,
                           method = "geommean")

    # Inject corrupted counts
    se <- .injectOutliersHyperparam(se,
                                    freq=freq,
                                    zScore=zScore,
                                    inj=inj,
                                    lnorm=lnorm,
                                    sdlog=sdlog)

    se <- calculateOffsets(se,
                           method = "TMM")

    # Fit edgeR model on intercept and known covariates. This is
    # further used to calculate the deviances.

    DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
    fit_DGE <- glmFit(y = DGE)

    # Principal component analysis on the deviance residuals
    deviances <- .CalculateDeviances(se = se,
                                     fit_DGE = fit_DGE)
    SingularVectors <- svd(t(deviances))

    assay(se,'deviances', withDimnames = FALSE) <- deviances

    metadata(se)[['svd_u']] <- SingularVectors$u
    metadata(se)[['svd_d']] <- SingularVectors$d
    metadata(se)[['svd_v']] <- SingularVectors$v

    metadata(se)[['optimalEncDim']] <- max(dimensions)

    se <- .orthonormalisation(se=se)



    beta_initial <- .initial_beta(se=se)
    eval <- bplapply(X=dimensions, ..., BPPARAM=BPPARAM,
                     FUN=function(i, ..., .evalAucPRLoss=NA){
                         .evalAutoCorrection(se = se,
                                             analysis = analysis,
                                             beta_initial = beta_initial,
                                             dimensions=i,
                                             BPPARAM=BPPARAM)})


    # Input the grid search results in the OutriderDataSet
    metadata(se)[['encDimTable']] <- data.table(
        encodingDimension= dimensions,
        evaluationLoss= unlist(eval),
        evalMethod='aucPR')

    #Obtain the optimal number of latent factors
    metadata(se)[['optimalEncDim']] <- NULL
    metadata(se)[['optimalEncDim']] <- .getBestQ(se)

    metadata(se)[['LatentFactorControl']] <- "DAE_AE"


    # Input the original counts again as count matrix, as some values were
    # changed with corrupted counts.
    se <- .clearHPoutput(se)

    return(se)
}


.DAE_AS <- function(se,
                    analysis,
                    aggregation = "locus",
                    dimensions=seq(
                        2, min(100, ncol(se) - 2, nrow(se) - 1), 2),
                    freq=1E-2,
                    zScore=3,
                    sdlog=log(1.6),
                    lnorm=TRUE,
                    inj='both',
                    BPPARAM=bpparam(),
                    ...){

    assay(se,'trueOffsets', withDimnames = FALSE) <- assays(se)$offsets

    # Find the optimal number of hyperparameters when using RUV with
    # deviance residuals.
    # Estimate size factors
    # Inject corrupted counts
    se <- .injectOutliersHyperparamSplicing(se,
                                            freq=freq,
                                            zScore=zScore,
                                            inj=inj,
                                            lnorm=lnorm,
                                            sdlog=sdlog)


    se <- calculateOffsets(se,
                           method = "proportion",
                           aggregation = "locus",
                           zeroCountOffsets = 1,
                           zeroOffsets = 1)


    # Fit edgeR model on intercept and known covariates. This is
    # further used to calculate the deviances.

    DGE <- .fitEdgeRDisp(se = se, design = getDesign(se))
    fit_DGE <- glmFit(y = DGE)

    # Principal component analysis on the deviance residuals
    deviances <- .CalculateDeviances(se = se,
                                     fit_DGE = fit_DGE)
    SingularVectors <- svd(t(deviances))

    assay(se,'deviances', withDimnames = FALSE) <- deviances

    metadata(se)[['svd_u']] <- SingularVectors$u
    metadata(se)[['svd_d']] <- SingularVectors$d
    metadata(se)[['svd_v']] <- SingularVectors$v

    metadata(se)[['optimalEncDim']] <- max(dimensions)

    se <- .orthonormalisation(se=se)


    # Calculate for each number of latent factors (= dimensions) that is
    # given as input the area under the curve of all merged p-values when
    # doing a RUV analysis with as true positives the corrupted counts.
    # This grid search is done to obtain the optimal number of latent factors
    # to include in the model.

    beta_initial <- .initial_beta(se=se)

    eval <- bplapply(X=dimensions, ..., BPPARAM=BPPARAM,
                     FUN=function(i, ..., .evalAucPRLoss=NA){
                         .evalAutoCorrection(se = se,
                                             analysis = analysis,
                                             beta_initial = beta_initial,
                                             dimensions=i,
                                             BPPARAM=BPPARAM)})

    # Input the grid search results in the OutriderDataSet
    metadata(se)[['encDimTable']] <- data.table(
        encodingDimension= dimensions,
        evaluationLoss= unlist(eval),
        evalMethod='aucPR')

    #Obtain the optimal number of latent factors
    metadata(se)[['optimalEncDim']] <- NULL
    metadata(se)[['optimalEncDim']] <- .getBestQ(se)

    metadata(se)[['LatentFactorControl']] <- "DAE_AS"

    # Input the original counts again as count matrix, as some values were
    # changed with corrupted counts; clear some outputs which were calculated
    # on these inputed counts.
    se <- .clearHPoutput(se)


    return(se)
}


.evalAutoCorrection <- function(se,
                                analysis,
                                dimensions,
                                beta_initial,
                                BPPARAM,
                                ...){
    # Search for corrupted counts by RUV with simultaneous parameter estimation
    se <- .fitRUV(se = se,
                  analysis = analysis,
                  dimensions = dimensions,
                  beta_initial = beta_initial,
                  padjust = "none",
                  ...)

    if(analysis == "AE"){
        eloss <- .evalAucPRLoss(se)

    } else if (analysis == "AS"){
        eloss <- .evalAucPRLossSplicing(se)
    }
    print(paste0('Evaluation loss: ', eloss,' for q=',dimensions))
    return(eloss)
}

.clearHPoutput <- function(se){

    assay(se,'deviances') <- NULL

    metadata(se)[['svd_u']] <- NULL
    metadata(se)[['svd_d']] <- NULL
    metadata(se)[['svd_v']] <- NULL

    metadata(se)[["Q"]] <- NULL
    metadata(se)[["R"]] <- NULL

    assay(se,'counts', withDimnames = FALSE) <- assay(se, 'trueCounts')
    assay(se,'offsets', withDimnames = FALSE) <- assay(se, 'trueOffsets')

    assay(se, 'trueCounts') <- NULL
    assay(se, 'trueOffsets') <- NULL

    assay(se, 'mu') <- NULL
    rowData(se)[["theta"]] <- NULL

    assay(se, 'pValue') <- NULL
    metadata(se)[['pValuesLocus']] <- NULL

    return(se)
}



.OptimalHardThreshold <- function(svd){
    d <- svd$d
    m <- dim(svd$u)[1]
    n <- dim(svd$v)[1]
    if(m > n){
        mtemp <- n
        n <- m
        m <- mtemp
    }
    dimensions <- sum(d >= (optimal_SVHT_coef(m/n,sigma_known = F) * median(d)))
    return(dimensions)
}



# Injection of corrupted counts based on a negative binomial distribution
# for the search of the optimal number of latent dimensions.
# Adapted from Outrider - Brechtmann et al.
.injectOutliersHyperparam <- function(se, freq, zScore, inj, lnorm, sdlog){

    # copy true counts to be able to acces them later
    assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)

    # generate index of injected corrupted counts
    size <- prod(dim(se))
    index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2),
                    replace = TRUE)
    index <- matrix(index, nrow = nrow(se))
    switch(inj,
           low = { index <- -abs(index) },
           high = { index <- abs(index) }
    )

    # Generate z-values for corrupted counts according to a lognormal
    # distribution if these are not specified
    tmpzScore <- matrix(0, ncol=ncol(se), nrow=nrow(se))
    if(isTRUE(lnorm)){
        tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
    } else {
        tmpzScore[index!=0] <- zScore
    }
    zScore <- tmpzScore

    # Define maximal count that should be injected
    max_out <- min(10*max(counts(se), na.rm=TRUE), .Machine$integer.max)

    # compute size factor normalized counts.
    # don't use it on the se to not influence the later calculation.
    se <- calculateOffsets(se, method = "geommean")
    sf_r <- getSizeFactors(se)
    # extract counts
    counts <- counts(se)
    # list of locations where corrupted counts have to be injected
    list_index <- which(index != 0, arr.ind = TRUE)
    # Loop over the list_index to inject corrupted counts

    DGE <- DGEList(counts=counts(se))
    DGE <- calcNormFactors(DGE)
    DGE <- estimateDisp(DGE) # estimate dispersion estimates
    fit_DGE <- glmFit(y = DGE)

    for(i in seq_len(nrow(list_index))){
        idxCol <- list_index[i,'col']
        idxRow <- list_index[i,'row']

        # Extract dispersion of negative binomial regression and define
        # the count that will be injected as corrupted count
        theta <- 1/fit_DGE$dispersion[idxRow]

        prob_z <- pnorm(zScore[idxRow, idxCol])

        if (index[idxRow, idxCol]==1){
            art_out <- qnbinom(p=prob_z, mu=counts(se)[idxRow,idxCol], size=theta)
        }
        else {
            art_out <- qnbinom(p=1-prob_z, mu=counts(se)[idxRow,idxCol], size=theta)
        }

        # only insert outliers if they are different from before
        # and not too large
        if(art_out < max_out & counts[idxRow, idxCol] != art_out){
            counts[idxRow, idxCol] <- art_out
        }else{
            index[idxRow, idxCol] <- 0
        }
    }

    # save the new count matrix, the index of corrupted counts and the
    # magnitude of the corrupted count via the z-value
    assay(se, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts),
                                                      nrow=nrow(se))
    assay(se, 'trueCorruptions', withDimnames=FALSE) <- index
    assay(se, 'injectedZscore', withDimnames=FALSE) <- zScore
    return(se)
}


# Adapted from FRASER - Mertes et al.
.injectOutliersHyperparamSplicing <- function(se, freq, zScore, inj, lnorm, sdlog, deltaMin = 0.2){
    # copy true counts to be able to acces them later
    assay(se, 'trueCounts', withDimnames=FALSE) <- counts(se)

    size <- length(unique(rowData(se)$locus))*ncol(se)
    index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2),
                    replace = TRUE)
    index <- matrix(index, ncol = ncol(se))
    rownames(index) <- unique(rowData(se)$locus)
    colnames(index) <- colnames(se)

    switch(inj,
           low = { index <- -abs(index) },
           high = { index <- abs(index) }
    )

    list_index <- which(index != 0, arr.ind = TRUE)
    binOutliers <- matrix(0,ncol = ncol(se),nrow = nrow(se))
    rownames(binOutliers) <- rownames(se)
    colnames(binOutliers) <- colnames(se)

    locusOutliers <- matrix(0, ncol = ncol(se), nrow = nrow(se))
    rownames(locusOutliers) <- rowData(se)$locus
    colnames(locusOutliers) <- colnames(se)

    for(i in seq_len(nrow(list_index))){
        idxCol <- list_index[i,'col']
        idxRow <- list_index[i,'row']

        locus_outlier <- rownames(index)[idxRow]

        binnames <- rownames(se)[rowData(se)$locus == locus_outlier]

        outlier_bin <- sample(binnames, 1)
        other_bins <- binnames[!(binnames %in% outlier_bin)]

        currentProportionOutlier <- assays(se)$counts[outlier_bin,idxCol] / assays(se)$offsets[outlier_bin,idxCol]
        ProportionOthers <- assays(se)$counts[other_bins,idxCol] / assays(se)$offsets[other_bins,idxCol]


        if(is.nan(currentProportionOutlier) | length(ProportionOthers) == 0){
            index[idxRow,idxCol] <- 0
        } else{

            if(index[idxRow,idxCol] == 1){
                deltaMax <- 1-currentProportionOutlier
            } else {
                deltaMax <- currentProportionOutlier
            }


            if(deltaMin > deltaMax){
                index[idxRow,idxCol] <- -index[idxRow,idxCol]

                if(index[idxRow,idxCol] == 1){
                    deltaMax <- 1-currentProportionOutlier
                } else {
                    deltaMax <- currentProportionOutlier
                }
            }

            DeltaProportionOutlier <- runif(n = 1, min = deltaMin, max = deltaMax)

            #outlier_count <- round((currentProportion + DeltaProportionOutlier*index[idxRow,idxCol])*(assays(se)$offsets[outlier_bin,idxCol]+2)-1,0)
            if(currentProportionOutlier != 1){
                DeltaOthers <- -DeltaProportionOutlier*index[idxRow,idxCol] * (ProportionOthers/(1-currentProportionOutlier))
            } else {
                DeltaOthers <- -DeltaProportionOutlier*index[idxRow,idxCol] / rep(1/length(other_bins), times = length(other_bins))

            }

            #other_counts <- round((ProportionOthers + DeltaOthers)*(assays(se)$offsets[other_bins,idxCol]+2)-1,0)
            if((currentProportionOutlier + DeltaProportionOutlier*index[idxRow,idxCol] < 0) ||
               ( ProportionOthers + DeltaOthers < 0)){
            }
            injectedCounts <- rmultinom(n = 1,
                                        size = assays(se)$offsets[other_bins,idxCol],
                                        prob = c(currentProportionOutlier + DeltaProportionOutlier*index[idxRow,idxCol],
                                                 ProportionOthers + DeltaOthers))
            rownames(injectedCounts) <- c(outlier_bin,other_bins)

            assays(se)$counts[rownames(injectedCounts),idxCol] <- injectedCounts
            binOutliers[outlier_bin,idxCol] <- index[idxRow,idxCol]
            locusOutliers[rownames(locusOutliers) == rowData(se)[outlier_bin,"locus"],idxCol] <- index[idxRow,idxCol]

        }
    }
    assay(se, 'binCorruptedCounts', withDimnames=FALSE) <- binOutliers
    assay(se, 'locusCorruptedCounts', withDimnames=FALSE) <- locusOutliers
    metadata(se)$locusCorruptedCounts <- index

    return(se)
}



# Adapted from OUTRIDER - Brechtmann et al.
.evalAucPRLoss <- function(se){
    # Calculate the precision-recall curve of all p-values, using corrupted
    # counts as true positives.
    scores <- -as.vector(assay(se, 'pValue'))
    labels <- as.vector(assay(se, 'trueCorruptions') != 0) + 0

    if(any(is.na(scores))){
        warning(sum(is.na(scores)), " P-values where NAs.")
        scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
    }
    pr <- pr.curve(scores, weights.class0=labels)
    return(max(0, pr$auc.integral, na.rm=TRUE))
}

# Adapted from OUTRIDER - Brechtmann et al.
.evalAucPRLossSplicing <- function(se){
    # Calculate the precision-recall curve of all p-values, using corrupted
    # counts as true positives.
    order <- match(rownames(metadata(se)$locusCorruptedCounts),rownames(metadata(se)$pValuesLocus))
    scores <- -as.vector( metadata(se)$pValuesLocus[order,])
    labels <- as.vector(metadata(se)$locusCorruptedCounts != 0) + 0

    if(any(is.na(scores))){
        warning(sum(is.na(scores)), " P-values where NAs.")
        scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
    }
    pr <- pr.curve(scores, weights.class0=labels)
    return(max(0, pr$auc.integral, na.rm=TRUE))
}


# Adapted from OUTRIDER - Brechtmann et al.
.getBestQ <- function(se){
    if('optimalEncDim' %in% names(metadata(se))){
        return(metadata(se)[['optimalEncDim']])
    }

    if('encDimTable' %in% names(metadata(se))){
        encTable <- metadata(se)[['encDimTable']]
        return(.getBestQDT(encTable, 'aucPR'))
    }
    # warning('Please find the optimal encoding dimension by running. ')
    return(NA_integer_)
}

# Adapted from OUTRIDER - Brechtmann et al.
.getBestQDT <- function(dt, usedEvalMethod='aucPR', digits=10){
    if('evalMethod' %in% colnames(dt)){
        testFun <- ifelse(all(dt[,"evalMethod" == usedEvalMethod]),
                which.max, which.min)
    } else {
    testFun <- which.max
    }

    dt[,"encodingDimension"][testFun(round(dt$evaluationLoss, digits))]
}
