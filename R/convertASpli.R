
#' @title Converting SummarizedExperiment with different ASpli-counts to
#' gene, bin or junction level counts.

#'

#' @description

#' convertASpli is used to obtain gene, bin or junction level counts from
#' a SummarizedExperiment object containing these three in the metadata slot.


#' @param ASpliSE An SummarizedExperiment object obtained by the
#' `BamtoAspliCounts` function, which contains gene, bin and junction level
#' counts in the metadata slot.
#'
#' @param type character vector specifying the counts to be extracted. Can be
#' `gene`, `bin` and `junction`.
#'
#' @param filter logical value specifying to filter bin counts based on the
#' type of bin. Regions that contain both intronic and exonic regions are
#' filtered as the intron and exon bins are already defined in other features.
#'
#' @param ... extra parameters for filtering.
#'

#' @return A `SummarizedExperiment` instance, representing the `gene`, `bin` or
#' `junction` counts as specified.
#'
#' @import ASpli
#' @import edgeR
#' @import S4Vectors
#' @import MASS
#' @import pracma
#' @import SummarizedExperiment
#' @import precrec
#' @import PRROC
#' @import BiocGenerics
#' @import methods
#' @import GenomicRanges
#' @import DESeq2
#' @import IRanges
#' @importFrom igraph graph_from_edgelist clusters
#' @importFrom rrcov PcaHubert
#' @importFrom limma lmFit strsplit2
#' @importFrom data.table data.table .N
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats model.matrix p.adjust pnbinom  pnorm  qnbinom rlnorm
#' rmultinom runif
#' @examples
#' data(saseRExample, package = "saseR")
#'
#' SEgenes <- convertASpli(ASpliSE, type = "gene")
#' SEbins <- convertASpli(ASpliSE, type = "bin")
#' SEjunctions <- convertASpli(ASpliSE, type = "junction")
#'
#' @export
#'
#'


convertASpli <- function(ASpliSE, type="none", filter = TRUE, ...){
    if(!(type %in% c("gene","bin","junction"))){
        message("Give type of counts that you want to analyse,
                e.g. 'gene', 'bin', 'junction'")}

    else if(type == "gene"){
        se <- .convertGene(ASpliSE)
        se <- se[rowSums(is.na(counts(se))) == 0, ]
    }
    else if(type == "bin"){
        se <- .convertBin(ASpliSE)
        if(filter == TRUE){
            se <- .FilterBinCounts(se, ...)
        }
        se <- se[rowSums(is.na(counts(se))) == 0, ]

    }
    else if(type == "junction"){
        se <-.convertJunction(ASpliSE)
        se <- se[rowSums(is.na(counts(se))) == 0, ]

    }
    return(se)
}




.convertGene <- function(ASpliSE){

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            "counts" =
               as.matrix(ASpliSE@metadata$geneCounts[,ASpliSE@colData$Names])),
        metadata = list(
            "DataType" = "gene"),
        rowData =
            DataFrame(
                ASpliSE@metadata$geneCounts[,
                            !(colnames(ASpliSE@metadata$geneCounts) %in%
                                                  ASpliSE@colData$Names)]),
        colData = DataFrame("Names" = ASpliSE@colData$Names)
    )

    return(se)

}


.convertBin <- function(ASpliSE){

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            "counts" =
                as.matrix(ASpliSE@metadata$binCounts[,ASpliSE@colData$Names])),
        metadata = list(
            "DataType" = "bin"),
        rowData =
            DataFrame(
                ASpliSE@metadata$binCounts[,
                                !(colnames(ASpliSE@metadata$binCounts) %in%
                                                 ASpliSE@colData$Names)]),
        colData = DataFrame("Names" = ASpliSE@colData$Names)
    )

    return(se)

}

.convertJunction <- function(ASpliSE){
    junctions <- metadata(ASpliSE)$junctionCounts
    number_of_samples <- length(colData(ASpliSE)$Names)

    start_J1 <- grep("StartHit", colnames(junctions)) + 1
    start_J2 <- grep("EndHit", colnames(junctions)) + 1
    start_J3 <- grep(ASpliSE@colData$Names[1],
                                  colnames(junctions))[1]
    end_J3 <- start_J3 + number_of_samples - 1

    J1 <- as.character(junctions$StartHit)
    J2 <- as.character(junctions$EndHit)
    J3 <- rownames(junctions)

    clusters <- .ASplimakeClusters(J1, J2, J3, bStrongFilter = FALSE)
    clustercounts <- .ASplimakeCountDataWithClusters(
        junctions[names(clusters$membership),start_J3:end_J3],
        clusters)

    rowData <- DataFrame(
        ASpliSE@metadata$junctionCounts[,
                            !(colnames(ASpliSE@metadata$junctionCounts) %in%
                                              ASpliSE@colData$Names)],
        "ASpliCluster" = NA)

    rowData$ASpliCluster[match(rownames(clustercounts), rownames(rowData))] <-
        clustercounts$locus

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(
            "counts" = as.matrix(
                ASpliSE@metadata$junctionCounts[,ASpliSE@colData$Names])),
        metadata = list(
            "DataType" = "junctions"),
        rowData = rowData,
        colData = DataFrame("Names" = ASpliSE@colData$Names)
    )

    return(se)

}


.FilterBinCounts <- function(se,
                            ignoreExternal = FALSE,
                            ignoreIo = TRUE,
                            ignoreI = FALSE){
    se <- se[ ! ignoreExternal | rowData(se)$event != "external" ,]
    se <- se[ ! ignoreIo | rowData(se)$feature != "Io" ,]
    se <- se[ ! ignoreI | rowData(se)$feature != "I" ,]

    return(se)
}

