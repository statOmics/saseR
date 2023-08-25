
#' @title Converting BAM files to an ASpli-SummarizedExperiment, which contains
#' gene, bin and junction counts.

#'

#' @description

#' BamtoAspliCounts is used to convert BAM files to a SummarizedExperiment
#' that contains gene, bin and junction level counts in the metadata slot.
#' This function is adapted from the ASpli package. More information can be
#' found in the corresponding package and paper.
#'
#' @param features An object of class ASpliFeatures, obtained by using the
#' binGenome function of ASpli.

#' @param targets A dataframe containing sample, bam and experimental factors
#' columns.
#'
#'
#' @param minReadLength Minimum read length of sequenced library. It is used
#' for computing E1I and IE2 read summarization. Make sure this number is
#' smaller than the maximum read length in every bam file, otherwise no E1I or
#' IE2 will be found
#'
#' @param maxISize Maximum intron expected size. Junctions longer than this
#' size will be dicarded.
#'
#' @param libType Defines how reads will be treated according their sequencing
#' library type (paired (PE, default) or single end (SE)).
#'
#' @param strandMode controls the behavior of the strand getter. It indicates
#' how the strand of a pair should be inferred from the strand of the first
#' and last alignments in the pair. 0: strand of the pair is always *.
#' 1: strand of the pair is strand of its first alignment.
#' This mode should be used when the paired-end data was generated using one of
#'  the following stranded protocols: Directional Illumina (Ligation), Standard
#'  SOLiD. 2: strand of the pair is strand of its last alignment.
#'  This mode should be used when the paired-end data was generated using one
#'  of the following stranded protocols: dUTP, NSR, NNSR, Illumina stranded
#'  TruSeq PE protocol. For more information see ?strandMode

#' @param minAnchor Minimum percentage of read that should be aligned to an
#' exon-intron boundary.
#'
#'
#' @param threshold Minimum number of reads supporting junctions. Default=5.
#'
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.
#'
#'
#' @return A `SummarizedExperiment` instance, containing `gene`, `bin` and
#' `junction` counts in the metadata slot.
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
#' @importFrom stats aggregate median model.matrix p.adjust pnbinom pnorm qnbinom rlnorm rmultinom runif
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
#'
#'
# Adapted from ASpli - Mancini et al.
BamtoAspliCounts <- function(
                           features,
                           targets,
                           minReadLength = 100,
                           maxISize = 50000,
                           libType="PE",
                           strandMode=0,
                           minAnchor = 0,
                           threshold = 5,
                           BPPARAM = bpparam()){


    counts <- .gbCounts(features = features,
                           targets = targets,
                           minReadLength = minReadLength,
                           maxISize = maxISize,
                           minAnchor = minAnchor,
                           libType=libType,
                           strandMode=strandMode,
                           BPPARAM = BPPARAM)


    jcounts <- ASpli::jCounts(counts= counts,
                               features = features,
                               minReadLength = minReadLength,
                               threshold = 5,
                               minAnchor = minAnchor,
                               libType=libType,
                               strandMode=strandMode)




    ASpliSE <- SummarizedExperiment::SummarizedExperiment(metadata =
                                                               list("geneCounts" = counts@gene.counts,
                                                                   "binCounts" = counts@exon.intron.counts,
                                                                    "junctionCounts" = jcounts@junctionsPJU,
                                                                    "IRCounts" = jcounts@junctionsPIR),
                                                           colData = list("Names" = rownames(targets)))
    return(ASpliSE)

}




# Adapted from ASpli - Mancini et al.
.gbCounts <- function( features, targets,  minReadLength,
                      maxISize, minAnchor = 10,
                      libType="SE",
                      strandMode=0,
                      BPPARAM = bpparam()) {

        counts <- .readCounts( features = features,
                          bam = NULL,
                          targets = targets,
                          readLength = minReadLength,
                          maxISize = maxISize,
                          minAnchor = minAnchor,
                          libType=libType,
                          strandMode=strandMode,
                          BPPARAM = BPPARAM)

    counts@.ASpliVersion = "2" #Marks ASpliCounts object with the ASpli update 2.0.0

    return(counts)
}

.readCounts <- function( features, bam, targets, cores = 1,
                        readLength,
                        maxISize,
                        minAnchor = 10,
                        libType=libType,
                        strandMode=strandMode,
                        BPPARAM = bpparam()) {


    if(!is.null(bam)){
        .Deprecated("gbCounts")
    }
    minReadLength <- readLength
    cores <- 1 #Allways use 1 core.

    #Create result object
    counts <- new(Class="ASpliCounts")
    counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.

    #Generates sample names in case there arent any
    targets <- ASpli:::.generateSamplesNames(targets)
    counts@targets <- ASpli:::.condenseTargetsConditions(targets) #ACH
    group                  <- counts@targets$condition
    counts@condition.order <- levels(factor( group, unique( group ), ordered = TRUE ))

    #Minimal anchors
    minAnchor <- if ( ! is.null(minAnchor) ) minAnchor else 10
    minA <- round( minAnchor * minReadLength / 100 )
    ptm <- proc.time()
    if(is.null(bam)) {
        ntargets <- nrow(targets)
    }else{
        ntargets <- 1
    }


    counts_list <- bplapply(X = c(1:ntargets), BPPARAM = BPPARAM,
                            FUN = function(i){
                                .bamToCounts(counts = counts,
                                            features = features,
                                            targets = targets[i,],
                                            minReadLength = minReadLength,
                                            maxISize = maxISize,
                                            minA = minA,
                                            libType=libType,
                                            strandMode=strandMode)})

    for (target in c(1:ntargets)){

        if(ncol(counts@gene.counts) == 0){
            counts@gene.counts <- counts_list[[target]]$gene.hits
        } else{
            counts@gene.counts <- cbind(counts@gene.counts,
                                        ASpli:::.extractCountColumns(counts_list[[target]]$gene.hits, targets[target, ]))
            colnames(counts@gene.counts)[ncol(counts@gene.counts)] <- rownames(targets)[target]
        }


        if(ncol(counts@exon.intron.counts) == 0){
            counts@exon.intron.counts <- counts_list[[target]]$exons.hits
        } else{
            counts@exon.intron.counts <- cbind(counts@exon.intron.counts,
                                               ASpli:::.extractCountColumns(counts_list[[target]]$exons.hits, targets[target, ]))
            colnames(counts@exon.intron.counts)[ncol(counts@exon.intron.counts)] <- rownames(targets)[target]
        }

        if(ncol(counts@e1i.counts) == 0){
            counts@e1i.counts <- counts_list[[target]]$e1i.hits
        }else{
            counts@e1i.counts <- cbind(counts@e1i.counts,
                                       ASpli:::.extractCountColumns(counts_list[[target]]$e1i.hits, targets[target, ]))
            colnames(counts@e1i.counts)[ncol(counts@e1i.counts)] <- rownames(targets)[target]
        }

        if(ncol(counts@ie2.counts) == 0){
            counts@ie2.counts <- counts_list[[target]]$ie2.hits
        }else{
            counts@ie2.counts <- cbind(counts@ie2.counts,
                                       ASpli:::.extractCountColumns(counts_list[[target]]$ie2.hits, targets[target, ]))
            colnames(counts@ie2.counts)[ncol(counts@ie2.counts)] <- rownames(targets)[target]
        }

        if(ncol(counts@junction.counts) == 0){
            counts@junction.counts <- counts_list[[target]]$junction.hits
            junction.hits <- counts_list[[target]]$junction.hits
        }else{
            dt1                    <- data.table(counts@junction.counts, keep.rownames = TRUE)
            dt2                    <- data.table(ASpli:::.extractCountColumns(counts_list[[target]]$junction.hits,
                                                                              targets[target, ]),
                                                 keep.rownames = T)
            dt3                    <- data.frame(merge(dt1, dt2, by="rn", all.x=T, all.y=TRUE))
            junction.hits <- counts_list[[target]]$junction.hits
            for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
                dt3[, s]           <- as.character(dt3[, s])
                junction.hits[, s] <- as.character(junction.hits[, s])
            }
            rownames(dt3)          <- dt3[, "rn"]
            dt3                    <- dt3[, -1]
            dt3[dt2$rn, 1:8]       <- ASpli:::.extractDataColumns(junction.hits, targets[target, ])
            counts@junction.counts <- dt3
            counts@junction.counts[is.na(counts@junction.counts)] <- 0
        }

        if(length(grep("NA", rownames(counts@junction.counts))) > 0){
            print(target)
            break
        }
        if(length(grep("NA", rownames(junction.hits ))) > 0){
            print(target)
        }
        gc()


    }



    for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
        counts@junction.counts[, s] <- as.factor(counts@junction.counts[, s])
    }
    colnames(counts@junction.counts)[9:ncol(counts@junction.counts)] <- rownames(targets)
    junctions.order <- sort(rownames(counts@junction.counts))
    junctions.order <- strsplit2(junctions.order, "[.]")
    junctions.order <- GRanges(seqnames=junctions.order[, 1], IRanges(start=as.numeric(junctions.order[, 2]), end=as.numeric(junctions.order[, 3])))
    junctions.order <- sort(junctions.order)
    junctions.order <- paste(junctions.order@seqnames, junctions.order@ranges@start, (junctions.order@ranges@start+junctions.order@ranges@width-1), sep=".")
    counts@junction.counts <- counts@junction.counts[junctions.order, ]

    # Create result object
    counts <- ASpli::rds( counts, targets )
    gc()
    return(counts)


}

.bamToCounts <- function(  counts = counts,
                          features = features,
                          targets = targets[i,],
                          minReadLength = minReadLength,
                          maxISize = maxISize,
                          minA = minA,
                          libType=libType,
                          strandMode=strandMode){

    bam <- ASpli::loadBAM(targets, cores = NULL,
                           libType=libType,
                           strandMode=strandMode)

    # Count Genes
    gene.hits <- ASpli:::.counterGenes( bam, ASpli::featuresg( features ))
    counts@gene.counts <- gene.hits

    # Count exons
    bins <- ASpli::featuresb( features )
    exons.hits <- ASpli:::.counterBin( bam, bins, gene.hits)
    counts@exon.intron.counts <- exons.hits



    # Count introns
    introns <- c( bins[ mcols(bins)$feature == "I" ],
                  bins[ mcols(bins)$feature == "Io"],
                  bins[ mcols(bins)$eventJ  == "IR"])

    # Count exon1 - intron regions
    e1i <- introns
    start( e1i ) <- start( introns ) - ( minReadLength - minA )
    end( e1i )   <- start( introns ) + ( minReadLength - minA )
    e1i.hits     <- ASpli:::.counterJbin(bam, e1i, gene.hits, minReadLength)
    counts@e1i.counts <- e1i.hits


    # Count intron - exon2 regions
    ie2 <- introns
    start( ie2 ) <- end( introns ) - ( minReadLength - minA )
    end( ie2 )   <- end( introns ) + ( minReadLength - minA )
    ie2.hits     <- ASpli:::.counterJbin( bam, ie2, gene.hits, minReadLength )
    counts@ie2.counts <- ie2.hits

    # Count junctions
    junction.hits    <- ASpli:::.counterJunctions( features, bam, maxISize )
    counts@junction.counts <- junction.hits


    return(list("gene.hits" = as.data.frame(gene.hits),
                "exons.hits" = as.data.frame(exons.hits),
                "e1i.hits" = as.data.frame(e1i.hits),
                "ie2.hits" = as.data.frame(ie2.hits),
                "junction.hits" = as.data.frame(junction.hits)))

}



jcounts <- function( counts,
                           targets,
                           features,
                           bam,
                           readLength,
                           threshold = 5,
                           cores = 1,
                           minAnchor = 10,
                           libType=libType,
                           strandMode=strandMode) {

        if(!.hasSlot(counts, ".ASpliVersion")){
            counts@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.
        }
        if(counts@.ASpliVersion == "1"){
            #Version conflict
            if(is.null(bam)){
                stop("Counts object is ASpli v1 but no bam was loaded. Please see vignette for new pipeline.")
            }
            .Deprecated("jCounts")
        }else{
            targets <- counts@targets
        }
        minReadLength <- readLength
        cores <- 1
        libType=libType
        strandMode=strandMode

        as  <- new(Class = "ASpliAS")
        as@.ASpliVersion = "1" #Last version before 2.0.0 was 1.14.0.
        as@targets <- targets

        #df0 <- countsj(counts)[ countsj(counts)$multipleHit == "-", ]
        #df0 <- df0[ df0$gene != "noHit" , ]

        targets <- ASpli:::.condenseTargetsConditions( targets )

        jcounts <- ASpli:::.filterJunctionBySample( df0=df0,
                                            targets=targets,
                                            threshold=threshold )

        # Junctions PSI:
        junctionsPSI    <- ASpli:::.junctionsPSI_SUM( df0, targets )
        as@junctionsPJU <- junctionsPSI
        message("Junctions PJU completed")

        # Junctions PIR:
        if(is.null(bam)) {
            ntargets <- nrow(targets)
            for(target in 1:ntargets){

                if(ntargets > 1){
                    #Load bam from current target
                    #agrego el libType y StrandMode
                    bam <- loadBAM(targets[target, ], cores = NULL,
                                   libType=libType, strandMode=strandMode)
                    junctionsPIR <- .junctionsDiscover( df=jcounts,
                                                        minReadLength=minReadLength,
                                                        targets=targets[target, ],
                                                        features=features,
                                                        minAnchor = minAnchor,
                                                        bam=bam)

                }



                if(ncol(as@junctionsPIR) == 0){
                    as@junctionsPIR <- junctionsPIR
                }else{
                    as@junctionsPIR <- cbind(as@junctionsPIR, junctionsPIR[, 3:6])
                }
            }
            junctions.order <- c(1, 2,
                                 seq(from=3, to=ncol(as@junctionsPIR), by=4),
                                 seq(from=4, to=ncol(as@junctionsPIR), by=4),
                                 seq(from=5, to=ncol(as@junctionsPIR), by=4))
            as@junctionsPIR <- as@junctionsPIR[, junctions.order]
            colnames(as@junctionsPIR)[c(-2, -1)] <- rep(rownames(targets), times=3)
            inicio_j1 <- 3
            inicio_j2 <- inicio_j1+nrow(targets)
            inicio_j3 <- inicio_j2+nrow(targets)
            j1 <- .sumByCond( as@junctionsPIR[, inicio_j1:(inicio_j1+nrow(targets)-1)],     targets )
            j2 <- .sumByCond( as@junctionsPIR[, inicio_j2:(inicio_j2+nrow(targets)-1)],     targets )
            j3 <- .sumByCond( as@junctionsPIR[, inicio_j3:(inicio_j3+nrow(targets)-1)],     targets )
            pirValues <- ( j1 + j2 ) / ( j1 + j2 + 2 * j3 )
            as@junctionsPIR <- cbind(as@junctionsPIR, pirValues)
        }else{

            junctionsPIR <- .junctionsDiscover( df=jcounts,
                                                minReadLength=minReadLength,
                                                targets=targets,
                                                features=features,
                                                minAnchor = minAnchor,
                                                bam=bam)
            as@junctionsPIR <- junctionsPIR
        }
        message("Junctions PIR completed")

        jranges <- .createGRangesExpJunctions( rownames( jcounts ) )

        # : refactor this code to other functions
        # ---------------------------------------------------------------------- #
        # Get all bins that are intronic or are associated to a Intron retention
        # event
        ic <- rbind( countsb(counts)[countsb(counts)$feature == "I",],
                     countsb(counts)[countsb(counts)$feature == "Io",],
                     countsb(counts)[countsb(counts)$event   == "IR*",],
                     countsb(counts)[countsb(counts)$event   == "IR",])
        # Get A GRanges object for intron bins, ordered by ic
        intranges <- featuresb(features)[ rownames(ic) ]

        # get exclusion junction counts, and make and index to ordered by ic
        dfe1e2 <- .e1e2JPIR( intranges, jcounts, targets )
        colnames(dfe1e2)[c(-1,-2)] <- rownames(targets)
        indexOrder <- match( dfe1e2$jbin, rownames( ic ) )

        # Get counts of inclusion junctions
        e1i <- .extractCountColumns( countse1i( counts ), targets )[ rownames(ic) ,]
        ie2 <- .extractCountColumns( countsie2( counts ), targets )[ rownames(ic) ,]

        j3 <- data.frame( matrix( NA,
                                  nrow =  nrow( e1i ),
                                  ncol =  length( targets$condition ) ),
                          stringsAsFactors = FALSE )
        colnames( j3 ) <- colnames( e1i )

        j3bin <- rep( NA , nrow( j3 ) )
        j3bin[ indexOrder ] <- rownames( dfe1e2 )
        j3[ indexOrder, ] <- .extractCountColumns( dfe1e2, targets )

        # Sum exclusion and inclusion counts by condition
        sumE1i <- .sumByCond( e1i, targets )
        sumIe2 <- .sumByCond( ie2, targets )
        sumJ3  <- .sumByCond( j3,  targets )

        # Calculates pir
        pirValues <- ( sumE1i + sumIe2 ) / ( sumE1i + sumIe2 + 2 * sumJ3 )

        # Creates result object
        result <- cbind(
            data.frame( event = ic$event ),
            data.frame( J1 = paste( rownames( e1i ), "E1I", sep="_") ),
            e1i,
            data.frame( J2 = paste( rownames( ie2 ), "IE2", sep="_") ),
            ie2,
            data.frame( J3 = j3bin ),
            j3,
            pirValues )


        message("Junctions IR PIR completed")

        as@irPIR <- result
        # ---------------------------------------------------------------------- #

        # ---------------------------------------------------------------------- #
        # Get all exons, except those that are associated to a intron retention
        # event
        ec <- countsb(counts)[countsb(counts)$feature == "E",]
        ec <- ec[ec$event != "IR",]
        ec <- ec[ec$event != "IR*",]

        exranges <- featuresb( features )[ rownames( ec ) ]

        fillAndReorderBy <- function( df , orderNames ) {
            indexOrder <- match( rownames( df ) , orderNames )
            result <- data.frame(
                matrix(
                    NA,
                    nrow = length( orderNames ),
                    ncol = ncol( df ) ) )
            result[ indexOrder, ] <- df
            colnames( result ) <- colnames( df )
            rownames( result ) <- orderNames
            return( result )
        }

        dfstart  <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'start' )
        dfstart  <- fillAndReorderBy( dfstart , rownames( ec ) )
        dfend    <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'end' )
        dfend    <- fillAndReorderBy( dfend , rownames( ec ) )
        dfwithin <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'within' )
        dfwithin <- fillAndReorderBy( dfwithin , rownames( ec ) )

        events   <- mcols( exranges ) $ event
        # ---------------------------------------------------------------------- #

        # ---------------------------------------------------------------------- #
        # Get the subset of previosly selected exons and gets only those associated
        # with an alternative splicing site usage event
        getAlternativeSS <- function( df, events ) {
            rbind(
                df[ events == "Alt3ss", ],
                df[ events == "Alt5ss", ],
                df[ events == "Alt3ss*", ],
                df[ events == "Alt5ss*", ] )
        }

        altJ1 <- getAlternativeSS( dfstart , events )
        altJ2 <- getAlternativeSS( dfend , events )
        altJ3 <- getAlternativeSS( dfwithin , events )
        colnames(altJ1)[-ncol(altJ1)] <- rownames(targets)
        colnames(altJ2)[-ncol(altJ2)] <- rownames(targets)
        colnames(altJ3)[-ncol(altJ3)] <- rownames(targets)

        sumAltJ1 <- .sumByCond( .extractCountColumns( altJ1, targets ), targets )
        sumAltJ1[is.na(sumAltJ1)] <- 0
        sumAltJ2 <- .sumByCond( .extractCountColumns( altJ2, targets ), targets )
        sumAltJ2[is.na(sumAltJ2)] <- 0
        sumAltJ3 <- .sumByCond( .extractCountColumns( altJ3, targets ), targets )
        sumAltJ3[is.na(sumAltJ3)] <- 0

        altPsiValues <- ( sumAltJ1 + sumAltJ2 ) / ( sumAltJ1 + sumAltJ2 + sumAltJ3 )

        result <- cbind(
            data.frame( event = mcols( exranges[ rownames( altJ1) ] )$ event ),
            data.frame( J1 = altJ1$overlappedSubjectNames ),
            .extractCountColumns( altJ1, targets ),
            data.frame( J2 = altJ2$overlappedSubjectNames ),
            .extractCountColumns( altJ2, targets ),
            data.frame( J3 = altJ3$overlappedSubjectNames ),
            .extractCountColumns( altJ3, targets ),
            altPsiValues )

        message("Junctions AltSS PSI completed")
        altPSI( as ) <- result
        # ---------------------------------------------------------------------- #

        # ---------------------------------------------------------------------- #
        # Get the subset of previosly selected exons and gets only those associated
        # with an exon skipping event and those not assigned to any splice event.
        getES <- function( df, events ) {
            rbind(
                df[ events == "ES", ],
                df[ events == "-", ],
                df[ events == "ES*", ] )
        }

        esJ1 <- getES( dfstart , events )
        esJ2 <- getES( dfend , events )
        esJ3 <- getES( dfwithin , events )
        colnames(esJ1)[-ncol(esJ1)] <- rownames(targets)
        colnames(esJ2)[-ncol(esJ2)] <- rownames(targets)
        colnames(esJ3)[-ncol(esJ3)] <- rownames(targets)

        sumEsJ1 <- .sumByCond( .extractCountColumns( esJ1, targets ), targets )
        sumEsJ1[is.na(sumEsJ1)] <- 0
        sumEsJ2 <- .sumByCond( .extractCountColumns( esJ2, targets ), targets )
        sumEsJ2[is.na(sumEsJ2)] <- 0
        sumEsJ3 <- .sumByCond( .extractCountColumns( esJ3, targets ), targets )
        sumEsJ3[is.na(sumEsJ3)] <- 0

        esPsiValues <- ( sumEsJ1 + sumEsJ2 ) / ( sumEsJ1 + sumEsJ2 + 2 * sumEsJ3 )

        result <- cbind(
            data.frame( event = mcols( exranges[ rownames( esJ1) ] )$ event ),
            data.frame( J1 = esJ1$overlappedSubjectNames ),
            .extractCountColumns( esJ1, targets ),
            data.frame( J2 = esJ2$overlappedSubjectNames ),
            .extractCountColumns( esJ2, targets ),
            data.frame( J3 = esJ3$overlappedSubjectNames ),
            .extractCountColumns( esJ3, targets ),
            esPsiValues )

        message("Junctions ES PSI completed")

        esPSI( as ) <- result
        # ---------------------------------------------------------------------- #

        # TODO: joint podria ser un getter, pero no es necesario mantener toda
        # esta data repetida
        joint( as ) <- rbind( altPSI( as ), esPSI( as ), irPIR( as ) )

        return( as )

    }


