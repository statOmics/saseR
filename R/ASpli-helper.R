
# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASpligenerateSamplesNames <- function (targets, collapse = "_")
{
    is.sequential <- function(x) {
        all(diff(x) == diff(x)[1])
    }
    my.make.unique <- function(s, sep = "_") {
        tab <- unique(s)
        tab <- setNames(rep(1, length(tab)), tab)
        sapply(s, function(ss) {
            sss <- paste(ss, tab[ss], sep = sep)
            tab[ss] <<- tab[ss] + 1
            return(sss)
        })
    }
    r <- suppressWarnings(as.numeric(rownames(targets)))
    if (all(!is.na(r))) {
        if (is.sequential(r)) {
            if (!"condition" %in% colnames(targets)) {
                targets <- .ASplicondenseTargetsConditions(targets,
                                                      collapse)
            }
            rownames(targets) <- my.make.unique(targets$condition,
                                                sep = collapse)
        }
    }
    return(targets)
}




# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASpliextractCountColumns <- function (aDataframe, targets)
{
    result <- aDataframe[, match(row.names(targets), colnames(aDataframe)),
                         F]
    colnames(result) <- as.character(row.names(targets))
    return(result)
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.
.ASplicounterJbin <- function (reads, feature, genes, l)
{
    cores = 1
    ungapped <- mclapply(reads, mc.cores = cores, function(x) {
        x[njunc(x) == 0, ]
    })
    hits <- mclapply(ungapped, mc.cores = cores, function(x) {
        co <- countOverlaps(feature, x, ignore.strand = TRUE,
                            minoverlap = l)
        gc()
        return(co)
    })
    jbinToGeneIndex <- match(feature@elementMetadata$locus, rownames(genes))
    gene_coordinates <- genes$gene_coordinates[jbinToGeneIndex]
    result <- data.frame(feature@elementMetadata[c("event", "locus",
                                                   "locus_overlap", "symbol")], gene_coordinates, as.data.frame(feature@ranges))
    for (i in 1:length(hits)) {
        result[, names(hits)[i]] <- hits[[i]]
    }
    result$names <- NULL
    colnames(result)[1:8] <- c("event", "locus", "locus_overlap",
                               "symbol", "gene_coordinates", "start", "end", "length")
    rownames(result) <- names(hits[[1]])
    return(result)
}



# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASplicondenseTargetsConditions <- function (targets, collapse = "_")
{
    if (!"condition" %in% colnames(targets)) {
        for (i in 2:ncol(targets)) {
            if (!is.character(targets[, i])) {
                targets[, i] <- as.character(targets[, i])
            }
        }
        targets <- data.frame(targets, condition = apply(targets[,
                                                                 -1, drop = FALSE], 1, paste, collapse = collapse),
                              stringsAsFactors = FALSE)
    }
    return(targets)
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASplicounterGenes <- function (reads, feature)
{
    cores = 1
    hits <- mclapply(reads, mc.cores = cores, function(x) {
        co <- countOverlaps(feature, x, ignore.strand = TRUE)
        gc()
        return(co)
    })
    effectiveGeneLength <- sum(width(feature))
    geneStarts <- sapply(start(feature), min)
    geneEnds <- sapply(end(feature), max)
    geneWidths <- geneEnds - geneStarts + 1
    strand <- as.character(unlist(runValue(strand(feature))))
    result <- data.frame(symbol = feature@elementMetadata$symbol,
                         locus_overlap = feature@elementMetadata$locus_overlap,
                         gene_coordinates = feature@elementMetadata$gene_coordinates,
                         start = geneStarts, end = geneEnds, length = geneWidths,
                         effective_length = effectiveGeneLength)
    for (i in 1:length(hits)) {
        result[, names(hits)[i]] <- hits[[i]]
    }
    rownames(result) <- names(hits[[1]])
    return(result)
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASplicounterBin <- function (reads, feature, genes)
{
    cores = 1
    hits <- mclapply(reads, mc.cores = cores, function(x) {
        co <- countOverlaps(feature, x, ignore.strand = TRUE)
        gc()
        return(co)
    })
    binToGeneIndex <- match(feature@elementMetadata$locus, rownames(genes))
    geneCoordinates <- genes$gene_coordinates[binToGeneIndex]
    result <- data.frame(feature@elementMetadata[c("feature",
                                                   "event", "locus", "locus_overlap", "symbol")], geneCoordinates,
                         as.data.frame(feature@ranges)[, c("start", "end", "width")])
    for (i in 1:length(hits)) {
        result[, names(hits)[i]] <- hits[[i]]
    }
    colnames(result)[1:9] <- c("feature", "event", "locus", "locus_overlap",
                               "symbol", "gene_coordinates", "start", "end", "length")
    rownames(result) <- names(hits[[1]])
    return(result)
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASpliextractDataColumns <- function(aDataframe, targets){
    result <- aDataframe[, -match(row.names(targets), colnames(aDataframe))]
    return(result)
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.

.ASplicounterJunctions <- function (features, bam, maxISize)
{
    cores = 1
    ujunctions <- mclapply(bam, mc.cores = cores, function(x) {
        junctions <- unlist(junctions(x))
        strand(junctions) <- "*"
        start(junctions) <- start(junctions) - 1
        end(junctions) <- end(junctions) + 1
        ujunctions <- unique(junctions)
        gc()
        return(ujunctions)
    })
    jranges <- unique(unlist(GRangesList(unlist(ujunctions))))
    maxWidth <- maxISize + 2
    jranges <- jranges[width(jranges) <= maxISize]
    fcoord <- paste(seqnames(jranges), start(jranges), end(jranges),
                    sep = ".")
    jranges@ranges@NAMES <- fcoord
    jcounts <- mclapply(bam, mc.cores = cores, function(x) {
        junctions <- unlist(junctions(x))
        strand(junctions) <- "*"
        start(junctions) <- start(junctions) - 1
        end(junctions) <- end(junctions) + 1
        count <- countMatches(jranges, junctions)
        jc <- data.frame(row.names = names(jranges), count)
        gc()
        return(jc)
    })
    df <- do.call("cbind", jcounts)
    colnames(df) <- names(jcounts)
    jranges <- .ASpliovBinJunction(features, jranges)
    jrdf <- data.frame(as.data.frame(jranges@elementMetadata$hitBin),
                       as.data.frame(jranges@elementMetadata$hitGen), as.data.frame(jranges@elementMetadata$hitGenStrand),
                       as.data.frame(jranges@elementMetadata$undef), as.data.frame(jranges@elementMetadata$symbol),
                       as.data.frame(jranges@elementMetadata$gene_coordinates),
                       as.data.frame(jranges@elementMetadata$bin_spanned), as.data.frame(jranges@elementMetadata$j_within_bin),
                       row.names = names(jranges))
    colnames(jrdf) <- c("junction", "gene", "strand", "multipleHit",
                        "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")
    rownames(jrdf) <- names(jranges)
    aa <- merge(jrdf, df, by.x = "row.names", by.y = "row.names",
                sort = FALSE)
    rnames <- paste(start(jranges) - 1, end(jranges) + 1, sep = "-")
    rownames(aa) <- fcoord
    aa$Row.names <- NULL
    return(aa)
}

.ASpliovBinJunction <- function (features, jranges)
{
    genes <- .ASpligetGeneGRanges(features)
    ovGene <- .ASpligetJunctionOverlapGeneData(jranges, genes)
    ambiguos <- ovGene[[1]]
    hitGen <- ovGene[[2]]
    hitGenStrand <- ovGene[[3]]
    gene_coordinates <- ovGene[[4]]
    symbol <- ovGene[[5]]
    hitBin <- .ASpligetJunctionMatchingBins(features, jranges)
    exonsBins <- featuresb(features)[featuresb(features)@elementMetadata$feature ==
                                         "E", ]
    span <- .ASpligetJunctionSpanningBins(jranges, exonsBins)
    j_within_bin <- .ASpligetJunctionWithinBins(jranges, exonsBins)
    mcols(jranges) <- append(mcols(jranges), DataFrame(hitBin = hitBin,
                                                       hitGen = hitGen, hitGenStrand = hitGenStrand, gene_coordinates = gene_coordinates,
                                                       undef = ambiguos, bin_spanned = span, j_within_bin = j_within_bin,
                                                       symbol = symbol))
    return(jranges)
}


.ASpligetGeneGRanges <- function (aspliFeatures)
{
    feature <- featuresg(aspliFeatures)
    aggregate_first <- function(data, by) {
        d = b = NULL
        data <- data.table(d = data, b = by)
        ans <- data[, list(A = first(d)), by = b]
        return(ans$A)
    }
    unlistedFeatures <- unlist(feature)
    geneAndChr <- paste(names(unlistedFeatures), as.character(seqnames(unlistedFeatures)),
                        sep = "_")
    geneChr <- aggregate_first(as.data.frame(seqnames(unlistedFeatures))[,
                                                                         1], by = geneAndChr)
    geneStarts <- aggregate(as.data.frame(start(unlistedFeatures)),
                            by = list(geneAndChr), FUN = min)[, 2]
    geneEnds <- aggregate(as.data.frame(end(unlistedFeatures)),
                          by = list(geneAndChr), FUN = max)[, 2]
    strand <- aggregate_first(as.data.frame(strand(unlistedFeatures))[,
                                                                      1], by = geneAndChr)
    geneCoordinates <- rep(feature@elementMetadata$gene_coordinates,
                           table(names(unlistedFeatures)))
    geneCoordinates <- aggregate_first(geneCoordinates, by = geneAndChr)
    symbols <- rep(feature@elementMetadata$symbol, table(names(unlistedFeatures)))
    symbols <- aggregate_first(symbols, by = geneAndChr)
    geneNames <- aggregate_first(as.data.frame(names(unlistedFeatures),
                                               stringsAsFactors = FALSE)[, 1], by = geneAndChr)
    genes <- GRanges(seqnames = geneChr, strand = strand, ranges = IRanges(geneStarts,
                                                                           geneEnds), gene_coordinates = geneCoordinates, symbol = symbols)
    mcols(genes)$names <- geneNames
    return(genes)
}


.ASpligetJunctionOverlapGeneData <- function (jranges, genes)
{
    hitGen <- rep("-", length(jranges))
    hitGenStrand <- rep("*", length(jranges))
    gene_coordinates <- rep("-", length(jranges))
    ambiguos <- rep("-", length(jranges))
    symbol <- rep("-", length(jranges))
    overGene <- findOverlaps(jranges, genes, type = "within")
    overGeneDF <- as.data.frame(overGene)
    posJrange <- overGeneDF$queryHits
    posGene <- overGeneDF$subjectHits
    overGeneDF$queryHits <- names(jranges)[as.numeric(overGeneDF$queryHits)]
    overGeneDF$subjectHits <- mcols(genes)$names[as.numeric(overGeneDF$subjectHits)]
    table <- table(overGeneDF$queryHits)
    if (nrow(overGeneDF) > 0) {
        ttG <- data.frame(aggregate(subjectHits ~ queryHits,
                                    data = overGeneDF, paste, collapse = ";"))
    }
    else {
        ttG <- data.frame(names = character(0))
        for (i in 1:ncol(overGeneDF)) {
            ttG[, i + 1] <- integer(0)
        }
        colnames(ttG)[2:ncol(ttG)] <- colnames(overGeneDF)
    }
    dd0 <- match(ttG$queryHits, names(jranges))
    hitGen[dd0] <- ttG$subjectHits
    dd <- match(ttG$queryHits, names(table))
    ttG$undef <- table[dd]
    ttG$tag <- rep("-", nrow(ttG))
    ttG$tag[ttG$undef > 1] <- "yes"
    ambiguos[dd0] <- ttG$tag
    hitGen[posJrange] <- mcols(genes)$names[posGene]
    hitGen[-posJrange] <- "noHit"
    hitGenStrand[posJrange] <- as.character(strand(genes)[posGene])
    gene_coordinates[posJrange] <- mcols(genes)$gene_coordinates[posGene]
    symbol[posJrange] <- as.character(genes@elementMetadata$symbol[posGene])
    return(list(ambiguos, hitGen, hitGenStrand, gene_coordinates,
                symbol))
}

.ASpligetJunctionMatchingBins <- function (features, jranges)
{
    hitBin <- rep("-", length(jranges))
    annJunctions <- featuresj(features)
    overJ <- findOverlaps(jranges, annJunctions, type = "equal")
    overJDF <- as.data.frame(overJ)
    namesJ <- as.numeric(overJDF[, 1])
    namesAnnJ <- as.numeric(overJDF[, 2])
    hitBin[namesJ] <- names(annJunctions[namesAnnJ])
    hitBin[-namesJ] <- "noHit"
    return(hitBin)
}


.ASpligetJunctionSpanningBins <- function (jranges, exonsBins)
{
    over <- findOverlaps(jranges, exonsBins)
    overDF <- as.data.frame(over)
    namesJ <- as.numeric(overDF[, 1])
    overDF[, 1] <- names(jranges[namesJ])
    namesBins <- as.numeric(overDF[, 2])
    overDF[, 2] <- names(exonsBins[namesBins])
    if (nrow(overDF) > 0) {
        tt <- data.frame(aggregate(subjectHits ~ queryHits, data = overDF,
                                   paste, collapse = ";"))
    }
    else {
        tt <- data.frame(names = character(0))
        for (i in 1:ncol(overDF)) {
            tt[, i + 1] <- integer(0)
        }
        colnames(tt)[2:ncol(tt)] <- colnames(overDF)
    }
    span <- rep("-", length(jranges))
    te <- match(names(jranges), tt$queryHits)
    span[!is.na(te)] <- tt$subjectHits[te][!is.na(te)]
    return(span)
}


.ASpligetJunctionWithinBins <- function (jranges, exonsBins) {
    j_within_bin <- rep("-", length(jranges))
    overJunctionWithinBins <- findOverlaps(jranges, exonsBins,
                                           type = "within")
    if (length(overJunctionWithinBins) > 0) {
        overJunctionWithinBinsDF <- as.data.frame(overJunctionWithinBins)
        namesJ <- as.numeric(overJunctionWithinBinsDF[, 1])
        namesB <- as.numeric(overJunctionWithinBinsDF[, 2])
        overJunctionWithinBinsDF[, 1] <- names(jranges[namesJ])
        overJunctionWithinBinsDF[, 2] <- names(exonsBins[namesB])
        agtt <- data.frame(aggregate(subjectHits ~ queryHits,
                                     data = overJunctionWithinBinsDF, paste, collapse = ";"))
        tw <- match(names(jranges), agtt$queryHits)
        j_within_bin[!is.na(tw)] <- agtt$subjectHits[tw][!is.na(tw)]
    }
    return(j_within_bin)
}

# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.
.ASplijunctionsPSI_SUM <- function (df, targets)
{
    getBoundarySharedData <- function(allJunctionCounts, junctionsGRanges,
                                      boundaryType = "start") {
        junctionNames <- rownames(allJunctionCounts)
        j.boundary <- findOverlaps(jranges, drop.self = TRUE,
                                   drop.redundant = FALSE, type = boundaryType)
        queryNames <- names(jranges[queryHits(j.boundary)])
        subjectNames <- names(jranges[subjectHits(j.boundary)])
        sharedNames <- data.frame(aggregate(subjectNames ~ queryNames,
                                            FUN = paste, collapse = ";"))
        junctionCounts <- data.frame(names = queryNames, .ASpliextractCountColumns(allJunctionCounts[subjectNames,
        ], targets), row.names = NULL)
        sharedRowSum <- data.frame(aggregate(. ~ names, data = junctionCounts,
                                             sum))
        colnames(sharedRowSum)[-1] <- rownames(targets)
        rownames(sharedRowSum) <- sharedRowSum$names
        sharedRowSum <- .ASpliextractCountColumns(sharedRowSum, targets)
        sharedRowSumOrdered <- data.frame(row.names = junctionNames,
                                          matrix(NA, nrow = nrow(allJunctionCounts), ncol = ncol(sharedRowSum)))
        colnames(sharedRowSumOrdered) <- colnames(sharedRowSum)
        orderIndex <- match(row.names(sharedRowSum), junctionNames)
        sharedRowSumOrdered[orderIndex, ] <- sharedRowSum
        sharedRowSumOrdered[is.na(sharedRowSumOrdered)] <- 0
        orderIndex <- match(sharedNames$queryNames, junctionNames)
        sharedNamesOrdered <- as.character(rep("-", nrow(sharedRowSumOrdered)))
        sharedNamesOrdered[orderIndex] <- sharedNames$subjectNames
        j1 <- .ASplisumByCond(.ASpliextractCountColumns(allJunctionCounts,
                                              targets), targets)
        j2 <- .ASplisumByCond(sharedRowSumOrdered, targets)
        jratioResult <- j1/(j1 + j2)
        colnames(jratioResult) <- paste(colnames(jratioResult),
                                        boundaryType, sep = ".")
        result <- data.frame(sharedNamesOrdered, sharedRowSumOrdered,
                             jratioResult, check.names = FALSE)
        colnames(result)[1] <- if (boundaryType == "start")
            "StartHit"
        else "EndHit"
        return(result)
    }
    junctionNames <- rownames(df)
    jranges <- sort(.ASplicreateGRangesExpJunctions(junctionNames))
    sharedStartData <- getBoundarySharedData(df, jranges, "start")
    sharedEndData <- getBoundarySharedData(df, jranges, "end")
    result <- do.call(cbind, list(df, sharedStartData, sharedEndData))
    return(result)
}


.ASplisumByCond <- function (countDf, targets)
{
    countDf[is.na(countDf)] <- 0
    uniqueConditions <- unique(targets$condition)
    nConditions <- length(uniqueConditions)
    result <- matrix(data = 0, nrow = nrow(countDf), ncol = nConditions)
    for (i in 1:nConditions) {
        result[, i] <- rowSums(countDf[, targets$condition ==
                                           uniqueConditions[i], drop = FALSE])
    }
    colnames(result) <- uniqueConditions
    return(result)
}

# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.
.ASplifilterJunctionBySample <- function (df0, targets, threshold)
{
    vecMin <- function(a, b) {
        return(pmin(a, b))
    }
    cropped <- .ASpliextractCountColumns(df0, targets)
    uniqueConditions <- unique(targets$condition)
    filter <- matrix(Inf, ncol = length(uniqueConditions), nrow = nrow(cropped))
    for (i in 1:length(uniqueConditions)) {
        byCond <- cropped[, targets$condition == uniqueConditions[i],
                          drop = FALSE]
        for (j in 1:ncol(byCond)) {
            filter[, i] <- vecMin(filter[, i], byCond[, j])
        }
    }
    filter <- rowSums(filter > threshold) > 0
    return(df0[filter, ])
}


# This function is copied from the ASpli package - Mancini et al. as it is not
# exported.
.ASplicreateGRangesExpJunctions <- function(jnames){
    split <- data.frame(matrix(unlist(strsplit(jnames, "[.]")),
                               byrow = TRUE, ncol = 3), row.names = jnames, stringsAsFactors = FALSE)
    colnames(split) <- c("chrom", "start", "end")
    split$start <- as.numeric(split$start)
    split$end <- as.numeric(split$end)
    split$chrom <- as.character(split$chrom)
    jranges <- GRanges(seqnames = split$chrom, ranges = IRanges(start = split$start,
                                                                end = split$end, names = jnames))
    return(jranges)
}


.ASpligetJPSIByOverlap <- function (jranges, exonGRanges, jcounts, targets, overlapType)
{
    jranges.edge <- jranges
    if (overlapType == "start") {
        start(jranges.edge) <- end(jranges.edge)
    }
    else if (overlapType == "end") {
        end(jranges.edge) <- start(jranges.edge)
    }
    overlapped <- findOverlaps(exonGRanges, jranges.edge, type = overlapType)
    overlappedQueryNames <- names(exonGRanges[queryHits(overlapped)])
    overlappedSubjectNames <- names(jranges.edge[subjectHits(overlapped)])
    df1 <- data.frame(names = overlappedQueryNames, .ASpliextractCountColumns(jcounts[subjectHits(overlapped),
    ], targets), row.names = NULL)
    result <- merge(x = data.frame(aggregate(. ~ names, data = df1,
                                             sum)), y = data.frame(aggregate(overlappedSubjectNames ~
                                                                                 overlappedQueryNames, FUN = paste, collapse = ";")),
                    by.x = 1, by.y = 1)
    rownames(result) <- result$names
    result[, 1] <- NULL
    return(result)
}



.ASplimakeClusters <- function (J1, J2, J3, bStrongFilter = FALSE)
{
    junctions_of_interest <- !is.na(J3) & J3 != "-"
    J1 <- J1[junctions_of_interest]
    J2 <- J2[junctions_of_interest]
    J3 <- J3[junctions_of_interest]
    JJ3 <- strsplit(J3, ";")
    main_junctions <- unlist(lapply(JJ3, function(s) {
        return(s[1])
    }))
    nJJ3 <- unlist(lapply(JJ3, length)) - 1
    JJ3 <- unlist(lapply(JJ3[nJJ3 > 0], function(s) {
        return(s[-1])
    }))
    resMed <- cbind(rep(main_junctions, times = nJJ3), JJ3)
    junctions_of_interest <- J1 != "-" & !is.na(J1)
    J1 <- J1[junctions_of_interest]
    J1 <- strsplit(J1, ";")
    nJ1 <- unlist(lapply(J1, length))
    J1 <- unlist(J1)
    resStart <- cbind(rep(main_junctions[junctions_of_interest],
                          times = nJ1), J1)
    junctions_of_interest <- J2 != "-" & !is.na(J2)
    J2 <- J2[junctions_of_interest]
    J2 <- strsplit(J2, ";")
    nJ2 <- unlist(lapply(J2, length))
    J2 <- unlist(J2)
    resEnd <- cbind(rep(main_junctions[junctions_of_interest],
                        times = nJ2), J2)
    if (bStrongFilter) {
        resStart <- resStart[which(resStart[, 1] %in% J3 & resStart[,
                                                                    2] %in% J3), ]
        resEnd <- resEnd[which(resEnd[, 1] %in% J3 & resEnd[,
                                                            2] %in% J3), ]
    }
    g <- graph_from_edgelist(rbind(resEnd, resStart, resMed),
                             directed = FALSE)
    gclus <- clusters(g)
    return(gclus)
}


.ASplimakeCountDataWithClusters <- function (countData, clusters)
{
    countData <- countData[names(clusters$membership), ]
    countData <- cbind(length = "", countData)
    countData <- cbind(end = "", countData)
    countData <- cbind(start = "", countData)
    countData <- cbind(gene_coordinates = "", countData)
    countData <- cbind(symbol = "", countData)
    countData <- cbind(locus_overlap = "", countData)
    countData <- cbind(locus = clusters$membership, countData)
    countData <- cbind(event = "", countData)
    countData <- cbind(feature = "", countData)
    countData <- countData[order(rownames(countData)), ]
    return(countData)
}

