#' @title calculating the offsets to perform aberrant expression or splicing
#' analysis.

#'

#' @description

#' calculateOffsets is used to calculate the offsets for aberrant expression or
#' splicing analysis. it can use ordinary offsets such as `edgeR` and `DESeq2`
#' offsets for aberrant expression analysis, or use the total counts aggregated
#' per gene to perform aberrant splicing analysis. It returns a
#' `SummarizedExperiment` which  includes an `offsets` matrix in the assays
#' slot.


#' @param se A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' When offsets are calculated for aberrant splicing (`AS`) or other
#' proportions (`proportion`), a rowData column is required which is used to
#' group features over which the sum is taken as offset.
#'
#' @param method method used to calculate the offsets. `TMM` uses edgeR's
#' trimmed mean of M-values, `geommean` uses DESeq2's size factors, `unit` uses
#' a unit matrix of offsets, `AS` and `proportion` use a rowData column to
#' aggregate feature counts (e.g. calculating the total bin counts per gene)
#' for aberrant splicing or other usage analyses.
#'
#' @param aggregation column of rowData used to aggregate feature counts. Only
#' used when `method` is `AS` or `proportion`. If not specified, the `locus`
#' column is used when available, followed by the `symbol` column.
#'
#' @param zeroCountOffsets When using the `AS` or `proportion` method, it is
#' possible that some offsets are 0. As offsets need to be a strict positive
#' value, these offsets are changed to a strict positive number (defined by
#' `zeroOffsets`, with a default of 1). The corresponding count can also be
#' changed, with a default change to 1. Note that, although it is possible,
#' it is not recommended to specify `zeroCountOffsets` greater than
#' `zeroOffsets`.
#'
#' @param zeroOffsets When using the `AS` or `proportion` method, it is
#' possible that some offsets are 0. As offsets need to be a strict positive
#' value, these offsets are changed to a strict positive number (defined by
#' `zeroOffsets`, with a default of 1). The corresponding count can also be
#' changed, with a default change to 1. Note that, although it is possible,
#' it is not recommended to specify `zeroCountOffsets` greater than
#' `zeroOffsets`.

#' @param mergeGeneASpli logical value. When ASpli is used to obtain feature
#' counts, one can aggregate feature counts to calculate offsets based on
#' both gene-level or ASpli-cluster. When `mergeGeneASpli` is TRUE, one uses
#' gene-level aggregation, except where no annotated gene is available, and
#' use ASpli-cluster aggregation for these features.

#' @param filterna logical value to filter NA locus values in rowData.
#' These are not necessarily features without annotated gene.
#'
#' @param saveall logical value to save all aggregated counts and offsets.
#'  Default is FALSE to save memory. If TRUE, also counts and offsets of
#'  gene-level aggregation and ASpli-cluster aggregation are saved, instead of
#'  only using the final used counts and offsets.


#' @return An updated `SummarizedExperiment` instance, now including an
#' `offsets` matrix in the assays slot. When `method` is `AS` or `proportion`,
#' an `realCounts` matrix is available in the assays slot, which correspond to
#' the original counts, while the `counts` matrix in the assays slot is updated
#' with to `zeroCountOffsets` when the calculated offset is equal to 0. When
#' `saveall` is TRUE, also matrices `Locuscounts`, `Locusoffsets`, `ASplicounts`
#' and `ASplioffsets` are saved in the assays slot, corresponding to the
#' adapted counts and offsets using the `locus` rowData column, and using the
#' `symbol` rowData column respectively.
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
#' @import IRanges
#' @import S4Vectors
#' @importFrom dplyr `%>%` group_by summarise across mutate
#' @importFrom rrcov PcaHubert
#' @importFrom limma lmFit strsplit2
#' @importFrom data.table data.table .N
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats model.matrix p.adjust pnbinom  pnorm  qnbinom rlnorm rmultinom runif
#' @examples
#'
#' data(saseRExample, package = "saseR")
#'
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
#' @export
#'
#'
calculateOffsets <- function(se,
                             method,
                             aggregation = NULL,
                             zeroCountOffsets = 1,
                             zeroOffsets = 1,
                             mergeGeneASpli = TRUE,
                             filterna = TRUE,
                             saveall = FALSE){
    if(!exists("method")){
        stop("Please select a method to calculate offsets.
             Examples are TMM, geommean, unit or proportion.")
    }
    if(method == "TMM"){ #edgeR TMM offsets

        DGE <- edgeR::DGEList(counts=counts(se))
        DGE <- edgeR::calcNormFactors(DGE)

        assay(se, 'offsets', withDimnames = FALSE) <- matrix(rep(c(DGE$samples$lib.size *
                                               DGE$samples$norm.factors),
                                         each = nrow(se)),
                                     ncol = ncol(se))
    } else if (method == "geommean"){ # DESeq2 size factors

        dse <- DESeqDataSet(se = se, design = getDesign(se))
        dse <- DESeq2::estimateSizeFactors(dse)
        assay(se,'offsets', withDimnames = FALSE) <-
            matrix(rep(exp(colData(dse)$sizeFactor),
                each = nrow(se)), ncol = ncol(se))

    } else if (method == "unit"){ # all offsets are set to 1

        assays(se)$offsets <- matrix(1, ncol = ncol(se), nrow = nrow(se))

    } else if (method %in% c("AS", "proportion")){ # Using the aggregation column to define proportions
        if(is.null(aggregation)){
            if(is.null(se@metadata$DataType))
                aggregation <- "locus"
            else if(se@metadata$DataType == "junctions"){
                aggregation <- "symbol"
            } else {
                aggregation <- "locus"
            }
        }
        if(!(aggregation %in% colnames(rowData(se)))){
            stop("Aggregation column not found in rowData.")
        }
        if(zeroOffsets <= 0){
            stop("zeroOffsets must be larger than 0.")
        }

        message(paste("Using ", aggregation, " column in rowData to calculate offsets.", sep = ""))
        offsets <- aggregate(assays(se)$counts,
                             by = list("locus" =
                                       as.factor(rowData(se)[,aggregation])),
                             FUN = sum)

        offsets <- merge(list("locus" = rowData(se)[,aggregation]),
                         offsets,
                         by = "locus")

        order <- match(rowData(se)[,aggregation], offsets[,"locus"])
        offsets <- offsets[order,] %>% mutate("locus" = NULL)

        colnames(offsets) <- colnames(se)
        rownames(offsets) <- rownames(se)


        assays(se)$realCounts <- assays(se)$counts

        assays(se)$counts[offsets == 0] <- zeroCountOffsets
        offsets[offsets == 0] <- zeroOffsets

        assays(se)$offsets <- as.matrix(offsets)

        rowData(se)$locus <- rowData(se)[,aggregation]

        if("ASpliCluster" %in% colnames(rowData(se))){
            if(se@metadata$DataType == "junctions" &
               aggregation == "symbol" &
               mergeGeneASpli == TRUE &
               sum(rowData(se)$symbol == "-" & (!is.na(rowData(se)$ASpliCluster))) != 0){
                if("ASpliCluster" %in% colnames(rowData(se))){

                    offsets <- aggregate(assays(se)$realCounts,
                                         by = list("locus" =
                                                       as.factor(rowData(se)[,"ASpliCluster"])),
                                         FUN = sum)

                    offsets <- merge(list("locus" = rowData(se)[,"ASpliCluster"]),
                                     offsets,
                                     by = "locus")

                    order <- match(rowData(se)[,"ASpliCluster"], offsets[,"locus"])
                    offsets <- offsets[order,] %>% mutate("locus" = NULL)

                    colnames(offsets) <- colnames(se)
                    rownames(offsets) <- rownames(se)

                    assays(se)$ASplicounts <- assays(se)$counts
                    assays(se)$ASplicounts[offsets == 0] <- zeroCountOffsets
                    offsets[offsets == 0] <- zeroOffsets

                    assays(se)$ASplioffsets <- as.matrix(offsets)

                    assays(se)$Locuscounts <- assays(se)$counts
                    assays(se)$Locusoffsets <- assays(se)$offsets

                    assays(se)$counts[rowData(se)$symbol == "-" &
                                          (!is.na(rowData(se)$ASpliCluster))] <- assays(se)$ASplicounts[rowData(se)$symbol == "-" &
                                                                                                            (!is.na(rowData(se)$ASpliCluster))]
                    assays(se)$offsets[rowData(se)$symbol == "-" &
                                          (!is.na(rowData(se)$ASpliCluster))] <- assays(se)$ASplioffsets[rowData(se)$symbol == "-" &
                                                                                                             (!is.na(rowData(se)$ASpliCluster))]
                    rowData(se)$locus <- rowData(se)$symbol
                    levels(rowData(se)$locus) <- c(levels(rowData(se)$locus),c(paste("ASpliCluster_", rowData(se)$ASpliCluster[rowData(se)$symbol == "-" &
                                                                                                                                 (!is.na(rowData(se)$ASpliCluster))], sep = "")))
                    rowData(se)$locus[rowData(se)$symbol == "-" &
                                          (!is.na(rowData(se)$ASpliCluster))] <- paste("ASpliCluster_", rowData(se)$ASpliCluster[rowData(se)$symbol == "-" &
                                                                                                                                     (!is.na(rowData(se)$ASpliCluster))], sep = "")
                    rowData(se)$BooleanASpliCluster <- rowData(se)$symbol == "-" & (!is.na(rowData(se)$ASpliCluster))

                    message(paste("Changed ", sum(rowData(se)$symbol == "-" &
                                                      (!is.na(rowData(se)$ASpliCluster))), " rows to ASpli cluster offsets.", sep = ""))
                    if(saveall == FALSE){ #remove output for memory saving
                        assays(se)$Locuscounts <- NULL
                        assays(se)$Locusoffsets <- NULL
                        assays(se)$ASplicounts <- NULL
                        assays(se)$ASplioffsets <- NULL


                    }
                }
            }
        }
        if(filterna == TRUE){
            filter <- (!is.na(rowData(se)$locus) & rowData(se)$locus != "-" & rowData(se)$locus != "noHit")
            se <- se[filter, ]
            message(paste(sum(!filter), " features filtered due to NA, '-', or 'noHit' locus.", sep = ""))
        }
    } else {
        stop("Offset calculation method not found.")
    }


    return(se)


}
