# This code is copied from the ASpli package, with a minor change at line 248, 
# which gives errors in the original package.

setGeneric ( name = "binGenome", 
             def = function( genome, geneSymbols = NULL, 
                             logTo = "ASpli_binFeatures.log", cores = 1 ) standardGeneric( "binGenome" ) )

setMethod(
    f = "binGenome",
    signature = "TxDb",
    definition = function ( genome, geneSymbols = NULL, logTo = "ASpli_binFeatures.log", cores = 1) {
        
        #Normalize seqnames. If . present in name, changes it to _ and warns the user
        if(length(grep("[.]", seqlevels(genome))) > 0){
            seqlevels(genome) <- gsub("[.]", "_", seqlevels(genome))
            warning("Some seqnames had a '.' present in their names. ASpli had to normalize them using '_'.")
        }
        
        features <- new( Class = "ASpliFeatures" )

        if ( is.null( geneSymbols ) ) {
            # Recupera los nombres de los genes
            geneSymbols <- data.frame( names( transcriptsBy(genome) ), stringsAsFactors = FALSE)
            row.names(geneSymbols) <- names( transcriptsBy(genome) )
            colnames(geneSymbols)  <- "symbol" 
        }
        if ( ! is.null ( logTo ) ) {
            sink( file = logTo )
        }
        genes.by.exons <- .createGRangesGenes( genome, geneSymbols ) 
        
        msg <- paste( "* Number of extracted Genes =" , length( genes.by.exons ) )
        message( msg )
        if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
        
        
        exon.bins <- .createGRangesExons(genome, geneSymbols)
        #add locus_overlap
        index <- match(exon.bins@elementMetadata$locus, names(genes.by.exons))
        locus_overlap <- rep("-", length(exon.bins))
        locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
        mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(locus_overlap=locus_overlap))
        
        msg <- paste( "* Number of extracted Exon Bins =", length( exon.bins ) )
        message( msg )
        if ( sink.number() > 0 ) cat( msg, "\n" ) 
        
        
        intron.tot <- .createGRangesIntrons(genome, geneSymbols)
        #add locus_overlap
        index <- match(intron.tot@elementMetadata$locus, names(genes.by.exons))
        locus_overlap <- rep("-", length(intron.tot))
        locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
        mcols(intron.tot) <- append( mcols(intron.tot), 
                                     DataFrame(locus_overlap=locus_overlap) )
        
        msg <- paste( "* Number of extracted intron bins =", length( intron.tot ) )
        message( msg )
        if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
        
        transcripts <- .createGRangesTranscripts(genome)
        
        msg <- paste( "* Number of extracted trascripts =" , length( unlist( transcripts ) ) )
        message( msg )
        if ( sink.number() > 0 ) cat( paste( msg,"\n") )
        
        junctions <- .createGRangesJunctions( genome ) 
        
        #add locus_overlap
        index <- match( junctions@elementMetadata$locus, names( genes.by.exons ) )
        locus_overlap <- rep( "-", length( junctions ) )
        locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[ index ]
        mcols(junctions) <- append( mcols( junctions ), 
                                    DataFrame( locus_overlap = locus_overlap ) )
        
        msg <- paste("* Number of extracted junctions =", length(junctions) )
        message( msg)
        if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
        
        intron.bins <- intron.tot[ intron.tot@elementMetadata$feature == "I"  ]
        intron.orig <- intron.tot[ intron.tot@elementMetadata$feature == "Io" ]
        
        class  <- rep("fullI", length( intron.orig ))
        mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( class=class ))
        event  <- rep("-", length( intron.orig ))
        eventJ <- rep("-", length( intron.orig ))
        mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( event=event ))
        mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( eventJ=eventJ ))
        
        exons <- exons( genome )
        
        exons.introns <- .findAsBin( exons, exon.bins, intron.bins, transcripts, 
                                     junctions, logTo )
        fullT <- c(exons.introns,intron.orig)
        
        transcriptExons <- exonsBy(genome, use.names=TRUE)
        
        tByGene    <- transcriptsBy(genome, by = "gene")
        
        a         <- unlist(tByGene)
        #geneNames <- DataFrame(gene=names(a),tx=mcols(a)$tx_name)
        geneNames <- DataFrame(gene=names(a))
        rownames(geneNames)<-mcols(a)$tx_name
        mcols( transcriptExons ) <- append ( mcols( transcriptExons ), geneNames[names(transcriptExons), , drop=FALSE] )

        features@genes <- genes.by.exons
        features@bins <- fullT
        features@junctions <- junctions
        features@transcriptExons <- transcriptExons
        
        #TODO hacer el disjoin directo con fullT y tratar de machear los rangos que no coinciden en el disjoin con 
        #los originales de fullT. Esto tiene el problema del overlap de locus y de perder la referencia de los bines
        #de donde salio el disjoin.
        
        # Start the clock!
        #ptm <- proc.time()
        
        locus <- unique(fullT$locus)
        fullTcorregido <- fullT[fullT$locus %in% locus]
        dtFullTcorregido <- as.data.table(fullTcorregido)
        lFullTcorregidos <- split(dtFullTcorregido, dtFullTcorregido$locus)
        
        message( "Correcting Io ends, this might take a while..." )
        #pb <- txtProgressBar(min=1,max=length(lFullTcorregidos),style=3)
        #i <- 0
        fullTcorregido <- pbmcapply::pbmclapply(lFullTcorregidos, mc.cores = cores, function(l){
            #print(l$locus[1])
            lIo    <- l[l$feature == "Io", ]
            lSinIo <- l[l$feature != "Io", ]
            rangos <- GRanges(seqnames = lSinIo$seqnames, ranges = IRanges(lSinIo$start, lSinIo$end, lSinIo$width), strand = lSinIo$strand)
            dis    <- disjoin(rangos)
            #Si los elementos disjuntos son menos que los rangos, es porque hay elementos que no deberian estar (ej AT1G02010). Los remuevo
            #Si no, simplemente ajusto los rangos
            if(length(dis) != length(rangos)){
                #Si los elementos disjuntos son mayores que los rangos, es porque no sabemos que pasa (ej AT2G31902).
                if(length(dis) > length(rangos)){
                    #print(l$locus[1])
                }else{
                    ol <- data.frame(findOverlaps(rangos, dis, type = "equal"))
                    lSinIo <- lSinIo[ol$queryHits]
                }
            }else{
                #Si es la misma cantidad simplemente alcanza con reemplazar los rangos anteriores (ej AT1G05860).
                dis <- data.frame(ranges(dis))        
                lSinIo$start <- dis$start
                lSinIo$end   <- dis$end
                lSinIo$width <- dis$width
            }
            return(rbind(lSinIo, lIo))
        })
        
        fullTcorregido <- rbindlist(fullTcorregido)
        fullTcorregido <- GRanges(seqnames      = fullTcorregido$seqnames, 
                                  ranges        = IRanges(start = fullTcorregido$start, end = fullTcorregido$end, width = fullTcorregido$width),
                                  strand        = fullTcorregido$strand, 
                                  locus         = fullTcorregido$locus, 
                                  bin           = fullTcorregido$bin, 
                                  feature       = fullTcorregido$feature, 
                                  symbol        = fullTcorregido$symbol, 
                                  locus_overlap = fullTcorregido$locus_overlap,
                                  class         = fullTcorregido$class, 
                                  event         = fullTcorregido$event, 
                                  eventJ        = fullTcorregido$eventJ)
        # Stop the clock
        #proc.time() - ptm
        names(fullTcorregido) <- paste0(fullTcorregido$locus, ":", fullTcorregido$bin)
        features@bins <- fullTcorregido
        message( "Genome binning completed" )
        return( features ) 
    })




.createGRangesGenes <- function( genomeTxDb, geneSymbols ) {
  
  # -------------------------------------------------------------------------- #
  # Se generan los rangos de los bins que corresponden a los bins exonicos.
  # Pimero se extraen los exones anotados para cada gen y despues se entrecruzan
  # los exones del mismo gen ( con la funcion disjoin ) dando los rangos de los
  # bins
  exons                  <- exonsBy( genomeTxDb, by="gen" ) #extract exons by gene
  exons.by.gene.disjoint <- disjoin( exons ) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrega datos de la coordenadas genomicas, el solapamiento de genes y 
  # simbolos de los genes
  geneCoordinates <- .createGRangesGenes.getGeneCoordinates( 
      exons.by.gene.disjoint , genomeTxDb )
  locusOverlap    <- .createGRangesGenes.getLocusOverlap( 
      exons.by.gene.disjoint )
  symbol          <- geneSymbols[ names( exons.by.gene.disjoint ), ]
  
  metadata        <- DataFrame( gene_coordinates = geneCoordinates,
      locus_overlap    = locusOverlap,
      symbol           = symbol ) 
  
  mcols( exons.by.gene.disjoint ) <- append(
      mcols( exons.by.gene.disjoint ), metadata )
  # -------------------------------------------------------------------------- #
  
  return( exons.by.gene.disjoint )
}

.createGRangesGenes.getGeneCoordinates <- function ( binsGRangesList , aTxDb) {
  
  geneChr    <- as.character( seqnames( genes( aTxDb ) ) )
  geneStarts <- sapply( start( binsGRangesList ), min )
  geneEnds   <- sapply( end(   binsGRangesList ), max )
  
  geneCoordinates <- paste0( geneChr, ':', geneStarts, '-', geneEnds )
  
  return( geneCoordinates )
}

# deprecated
.createGRangesGenes.getGeneCoordinates.old <- function ( binsGRangesList ) {
  
  chromosome                <- sapply( seqnames( binsGRangesList ),unique )
  chromosome                <- as.character( chromosome )
  
  geneStarts             <- sapply( start( binsGRangesList ), min )
  geneEnds               <- sapply( end( binsGRangesList ), max )
  
  browser                <- paste( chromosome, geneStarts, sep=":" ) 
  
  geneCoordinates <- paste( browser, geneEnds, sep="-" )
  
  return (geneCoordinates)
  
}

.createGRangesGenes.getLocusOverlap <- function ( binsGRangesList ) {
  
  # -------------------------------------------------------------------------- #
  # Busca los rangos que se solapan, debido a que se espera que los rangos de 
  # binsGRangesList sean disjuntos ( es decir, no se solapan para un mismo gen),
  # los solapamientos que se encuentran deben ser entre rangos de genes 
  # diferentes.
  #Modificaciones hechas por AR 20201128
  if(F){
    if ( packageVersion("IRanges") < 2.6 ) {
      locus <- findOverlaps( binsGRangesList, ignoreSelf=TRUE, ignore.strand = T ) 
    } else {
      locus <- findOverlaps( binsGRangesList, drop.self=TRUE, ignore.strand = T ) 
      #locus <- findOverlaps( binsGRangesList[c(i, j)], drop.self=TRUE, ignore.strand = T ) 
    }
    
    
    # -------------------------------------------------------------------------- #
    
    locusOverlap <- rep("-", length( binsGRangesList ) )  
    
    # -------------------------------------------------------------------------- #
    # Si hay genes que se solapan, se recuperan todos los nombres de los genes
    if ( length( locus ) > 0 ) {
      names      <- names( binsGRangesList[ to( locus ) ] )
      geneIndex  <- from ( locus ) 
      aggregated <- data.frame( aggregate( names ~ geneIndex, 
                                           FUN = paste, collapse = ';' ) )
      locusOverlap[ aggregated$geneIndex ] <- aggregated$names
    }
    # -------------------------------------------------------------------------- #
  }else{
    starts   <- aggregate(start ~ group_name, data.frame(binsGRangesList), min)
    ends     <- aggregate(end ~ group_name, data.frame(binsGRangesList), max)
    seqnames <- aggregate(seqnames ~ group_name, data.frame(binsGRangesList), unique)
    strands  <- aggregate(strand ~ group_name, data.frame(binsGRangesList), unique)
    rangos   <- merge(strands, merge(seqnames, merge(starts, ends, by = "group_name"), by = "group_name"), by = "group_name")
    grangos  <- makeGRangesListFromDataFrame(rangos, names.field = "group_name")
    names(grangos) <- rangos$group_name

    locus <- findOverlaps( grangos, drop.self=TRUE, ignore.strand = T ) 
    
    # -------------------------------------------------------------------------- #
    
    locusOverlap <- rep("-", length( grangos ) )  
    
    # -------------------------------------------------------------------------- #
    # Si hay genes que se solapan, se recuperan todos los nombres de los genes
    if ( length( locus ) > 0 ) {
      names      <- names( grangos[ to( locus ) ] )
      geneIndex  <- from ( locus ) 
      aggregated <- data.frame( aggregate( names ~ geneIndex, 
                                           FUN = paste, collapse = ';' ) )
      locusOverlap[ aggregated$geneIndex ] <- aggregated$names
    }
    # -------------------------------------------------------------------------- #    
    
  }
  
  
  
  return ( locusOverlap )
  
}
.createGRangesExons <- function( aTxDb, geneSymbols ) {
  
  # -------------------------------------------------------------------------- #
  # Recupera los exones de los genes con mas de un exon
  exonsM <- .createGRangesExons.getExonsFromMultiExonicGenes( aTxDb )
  # -------------------------------------------------------------------------- # 
  
  # -------------------------------------------------------------------------- #
  # Genera los bins para los exones seleccionados
  exon.bins <- .createGRangesExons.getExonBins ( exonsM )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Genera nombres para los bins creados
  exonData <-  .createGRangesExons.getExonBinNames ( exon.bins )
  exon.bins.names          <- exonData [[1]] 
  exon.gene.names.repeated <- exonData [[2]] 
  exon.bins.id             <- exonData [[3]] 
  # Setea los nombres de los bins
  names( exon.bins ) <- exon.bins.names
  # -------------------------------------------------------------------------- #  
  
  # -------------------------------------------------------------------------- #
  # Agrega metadatos
  feature.exon    <- rep("E", length(exon.bins.id))
  geneSymbolIndex <- match( exon.gene.names.repeated, rownames( geneSymbols ) )  
  symbol          <- geneSymbols[geneSymbolIndex, ]
  
  metadata <- DataFrame( locus   = exon.gene.names.repeated, 
      bin     = exon.bins.id,
      feature = feature.exon,
      symbol  = symbol )
  
  mcols( exon.bins ) <- append ( mcols( exon.bins ), metadata )
  # -------------------------------------------------------------------------- #
  
  return (exon.bins)
}

.createGRangesExons.getExonBinNames <- function ( aGRanges ) {
  
  # -------------------------------------------------------------------------- #
  # Recupera el numero de exones para cada gen.
  nExonsByGene <- table( names( aGRanges ) )          
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Genera un lista que contiene el numero de orden de cada exon
  exon.bins.num <- unlist( lapply( nExonsByGene, seq ) )
  # Genera identificadores de orden para los bins
  exonBinIds  <- sprintf('E%03d', exon.bins.num)
  # Genera los nombre de los bins
  exonGeneNames <- rep( names(nExonsByGene), nExonsByGene) 
  exonBinNames <- paste(exonGeneNames, exonBinIds ,sep=":")
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Devuelve los nombre de los bins, genes y los identificadores de orden de los
  # bins
  return ( list( exonBinNames, exonGeneNames, exonBinIds) )
  # -------------------------------------------------------------------------- #
  
}

.createGRangesExons.getExonsFromMultiExonicGenes <- function ( aTxDb ) {
  
  # Extrae todos los exones
  exons <- exonsBy(aTxDb, by="gen" ) #extract exons by gen 
  # se queda con los GRanges de los genes que tienen mas de un exon
  multipleExons <- exons[ lengths( exons )  > 1 ] 
  
  return ( multipleExons )
  
}

.createGRangesExons.getExonBins <- function ( aGRangesList ) {
  
  # Genera los segmentos que corresponden a los bins
  exon.by.gene.disjoint <- disjoin( aGRangesList )
  
  # Convierte el GRangeList a un GRanges.
  # Se mezclan los GRanges de todos lo genes, es una sola lista
  exon.by.gene.disjoint.unlist <- unlist( exon.by.gene.disjoint )
  
  # Se eliminan GRanges duplicados, esto es porque pueden tener origen en 
  # distintos genes.
  exon.bins <- exon.by.gene.disjoint.unlist[ ! duplicated( exon.by.gene.disjoint.unlist ) ] 
  
  return( exon.bins )
  
}

# ---------------------------------------------------------------------------- #
# .createGRangesIntrons
.createGRangesIntrons <- function( aTxDb, geneSymbols ) {
  
  # -------------------------------------------------------------------------- #
  # Recupera los intrones de la anotacion por transcripto
  introns <- intronsByTranscript( aTxDb, use.names=TRUE )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Convierte el GRangesList de los intrones por transcripto en un GRanges de
  # todos los intrones solamente.
  introns.ulst <- unlist( introns )
  # Elimina intrones duplicados
  introns.no.dups <- introns.ulst[ ! duplicated( introns.ulst ) ]
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Recupera los nombres de los genes de cada transcripto
  introns.tx.names <- names( introns.no.dups ) 
  suppressMessages ( 
      introns.map <- ensembldb::select( 
          aTxDb, 
          keys    = introns.tx.names, 
          keytype = 'TXNAME', 
          columns = 'GENEID' )
  )
  # Elimina posibles datos con nombres de gen no valido.
  introns.map <- introns.map[ ! is.na( introns.map$GENEID ) , ]
  introns.no.dups <-  introns.no.dups[ ! is.na( introns.map$GENEID ) , ]
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrupa los intrones por gen. Devuelve un GRangesList
  introns.by.gene <- split( introns.no.dups, introns.map$GENEID ) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Recupera todos los rangos del GrangesList en un unico GRanges.
  nIntronsByGene      <- lengths( introns.by.gene )
  intronOrig          <- unlist( introns.by.gene )
  geneNames           <- rep( names( introns.by.gene ), nIntronsByGene )
  names( intronOrig ) <- geneNames
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Procesa los intrones originales
  intronOrigIndex   <- sprintf( 'Io%03d', unlist( lapply( nIntronsByGene, seq ) ) )
  
  names(ranges(intronOrig)) <- paste( geneNames, intronOrigIndex, sep=":" ) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrega metadatos
  feature.o.intron  <- rep( "Io", length( intronOrigIndex ) )
  
  metadata <- DataFrame( locus   = geneNames ,
      bin     = intronOrigIndex,  
      feature = feature.o.intron )
  
  mcols( intronOrig ) <- append( mcols(intronOrig) , metadata)
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Procesa los bins intronicos.
  # Obtiene los rangos de los bins.
  introns.by.gene.disjoint <- disjoin( introns.by.gene )
  # Recupera el numero de intrones por gen
  nIntronsByGene           <- lengths( introns.by.gene.disjoint )
  # Recupera los nombres del gen para cada intron
  geneNames                <- rep (names(introns.by.gene.disjoint),nIntronsByGene )
  # Junta todos los rangos de GRangeList de los intrones en un unido GRanges.
  intron.bins              <- unlist( introns.by.gene.disjoint )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Genera los identificadores de cada intron
  intron.bins.num          <- sprintf( 'I%03d', 
      unlist( lapply( nIntronsByGene, seq ) ) )
  # Agrega los nombes de los genes a cada intron
  intron.bin.names         <- paste( geneNames , 
      intron.bins.num, sep=":" ) 
  names( intron.bins )     <- intron.bin.names  
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrega metadatos
  
  feature.intron <- rep("I", length(intron.bins.num) )
  
  
  metadata <- DataFrame( locus   = geneNames,       
      bin     = intron.bins.num,  #add metadata    
      feature = feature.intron  )  #add metadata 
  
  mcols(intron.bins) <- append( mcols( intron.bins ), metadata )                    
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Junta los intrones originales y los bins intronicos.
  # Si hay un intron originales que es igual a uno intronico, se elimina.
  # La eliminacion se hace implicitamente con la funcion unique. Se eliminan
  # los originales porque estan despues que los bins intronicos.
  intrones.totales  <- c( intron.bins, intronOrig )
  intron.tot.u      <- sort( unique( intrones.totales ) )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrega los simbolos de los genes
  symbolIndex <- match( mcols(intron.tot.u)[, 'locus'], rownames( geneSymbols ) )
  symbol <- geneSymbols[symbolIndex, ]
  mcols(intron.tot.u) <- append( mcols( intron.tot.u ), DataFrame( symbol = symbol ) )
  # -------------------------------------------------------------------------- #
  
  return (intron.tot.u)
}
# Fin de createGRangesIntrons
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .createGRangesTranscripts
.createGRangesTranscripts <- function(genome) {
  transcripts <- transcriptsBy(genome)
  return(transcripts)
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .createGRangesJunctions
.createGRangesJunctions <- function( aTxDb ) {
  
  # -------------------------------------------------------------------------- #
  # Recupera los intrones por transcripto
  introns.no.dups <- .createGRangesJunctions.getUniqueIntrons( aTxDb )
  # Recupera los nombres de los transcriptos para cada intron
  introns.tx.names <- names( introns.no.dups )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Agrupa los intrones por gen
  introns.by.gene <- .createGRangesJunctions.intronsByGene( aTxDb, 
      introns.tx.names , introns.no.dups )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Genera las junturas a partir de los intrones 
  nIntronsByGene <- lengths( introns.by.gene )
  geneNames      <- rep( names( introns.by.gene ) , nIntronsByGene )
  junctions <- .createGRangesJunctions.getJunctions( introns.by.gene, nIntronsByGene, geneNames  )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Modifica el rangos de las junturas
  start( junctions ) <- start( junctions ) - 1
  end( junctions )   <- end( junctions ) + 1
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Agrega metadatos
  mcols( junctions) <- append( mcols( junctions ), DataFrame( locus = geneNames ) ) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Ordena las junturas
  junctions <- sort( junctions )
  # -------------------------------------------------------------------------- #
  
  return (junctions)
  
}

.createGRangesJunctions.getUniqueIntrons <- function ( aTxDb ) {
  
  introns <- intronsByTranscript( aTxDb, use.names=TRUE )
  
  introns.ulst <- unlist( introns )          
  
  introns.no.dups <- introns.ulst[ !duplicated( introns.ulst ) ]
  
  return ( introns.no.dups  ) 
}

.createGRangesJunctions.intronsByGene <- function ( aTxDb, transcripts, introns.no.dups ) {
  suppressMessages(
      introns.map <- ensembldb::select( aTxDb, 
          keys    = transcripts, 
          keytype = 'TXNAME', 
          columns = 'GENEID' )
  )
  introns.map <- introns.map[ !is.na( introns.map$GENEID ),  ]
  introns.no.dups <- introns.no.dups[ ! is.na( introns.map$GENEID ) ]
  # Agrupa los intrones por gen
  introns.by.gene <- split( introns.no.dups, introns.map$GENEID )
  
  return ( introns.by.gene )
}

.createGRangesJunctions.getJunctions <- function ( intronsByGene, nIntronsByGene, geneNames  ) {
  junctions          <- unlist( intronsByGene ) 
  introns.num        <- sprintf( 'J%03d',  unlist( lapply( nIntronsByGene, seq ) ) )
  junction.names     <- paste( geneNames, introns.num, sep = ":" ) 
  names( junctions ) <- junction.names
  return( junctions )
}
# Fin de .createGRangesJunctions
# ---------------------------------------------------------------------------- #

.shiftVector <- function ( x, positions,  dir="right", default = NA ) {
  invert    <- positions < 0
  positions <- if (invert) -positions else positions
  positions <- positions %% length(x)
  
  dir <- tolower( dir )
  dir <- if (dir == 'right') 0 else if (dir == 'left') 1 else NA
  dir <- if (invert) dir <- 1 - dir else dir
  
  if ( dir == 0 ) {
    return( c( rep( default, times=positions ), x[1:(length(x)-positions )] ) )
  } else 
  if ( dir == 1 ) {
    return( c( x[(positions+1):length(x)], rep( default, times=positions )) )
  }
  warning("Direction no recognized. Return unmodified value.")
  return (x)
}


# ---------------------------------------------------------------------------- #
# .findAsBin 
.findAsBin <- function( exons, exon.bins, intron.bins, transcripts, junctions, logTo = NULL) {
    
    # -------------------------------------------------------------------------- #
    # Get all alternative bins. I.e. those that are shared between exons bins and
    # intron bins.
    as.bins <- findOverlaps( exon.bins, intron.bins, type = c( "equal" ) ) 
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Add 'as' tag to alternative bins
    exon.as   <- rep( "-", length( exon.bins ) ) 
    intron.as <- rep( "-", length( intron.bins ) ) 
    exon.as[ from( as.bins ) ] <- "Undefined AS" 
    intron.as[ to( as.bins ) ] <- "Undefined AS" 
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Looks for terminal bins.
    # A bin has one 5' and 3' terminal bin for each transcript.
    transcripts.unlist  <- unlist( transcripts )
    find.external.start <- findOverlaps( exon.bins, transcripts.unlist, type=c( "start" ) ) 
    find.external.end   <- findOverlaps( exon.bins, transcripts.unlist, type=c( "end" ) )
    # Add the 'terminal' tag to a bin. If a bin has a 'as' tag already, the new 
    # 'tag' overwrites the older one.
    exon.as[ find.external.start@from ] <- "external" 
    exon.as[ find.external.end@from ] <- "external" 
    # -------------------------------------------------------------------------- #
    
    
    # -------------------------------------------------------------------------- #
    # Add 'as' and 'external' tags to introns
    mcols(intron.bins) <- append( mcols( intron.bins ), DataFrame(class=intron.as)) 
    # and exons.
    mcols(exon.bins)   <- append( mcols( exon.bins ), DataFrame(class=exon.as)) 
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Join all bins
    bins <- unique ( sort( append( exon.bins, intron.bins ) ) ) 
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Assign a splice event to each alternative bin y putative junctions if 
    # possible.
    auxdf <- data.frame( feature = as.character( mcols( bins )$feature), 
                         class   = as.character( mcols( bins )$class),
                         strand  = as.character( strand( bins) ) )
    # Preassign default values            
    auxdf$events <- as.character(auxdf$class)
    auxdf$eventsJ <- as.character(auxdf$class)
    # -------------------------------------------------------------------------- #
    
    
    # -------------------------------------------------------------------------- #
    # Add class information of neighbor bins.
    auxdf$classPrev <- .shiftVector(as.character(auxdf$class),1, default = '-')
    auxdf$classNext <- .shiftVector(as.character(auxdf$class),-1, default = '-')
    isAS <- auxdf$class == 'Undefined AS'
    neighbourIsAS <- auxdf$classPrev =='Undefined AS' | auxdf$classNext == 'Undefined AS'
    # -------------------------------------------------------------------------- #
    
    
    # -------------------------------------------------------------------------- #
    # Looks for bins that overlap perfectly with a junction
    ji <- junctions
    start(ji) <- start(junctions)+1
    end(ji) <- end(junctions)-1
    
    over <- findOverlaps( bins, ji, type="equal")
    binOverlapsJunction <- rep( FALSE ,  nrow(auxdf) ) 
    binOverlapsJunction[queryHits(over)] <- TRUE
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Looks for bins that overlap perfectly with an exon
    over <- findOverlaps( bins, exons, type="equal")
    binOverlapsExon <- rep( FALSE ,  nrow(auxdf) ) 
    binOverlapsExon[queryHits(over)] <- TRUE
    # -------------------------------------------------------------------------- #
    
    
    # -------------------------------------------------------------------------- #
    # Looks for bins that overlaps with the start or the end of a junction 
    jranges.start <- junctions
    start(jranges.start) <- end(jranges.start)
    st <- findOverlaps(type = "start", bins, jranges.start )
    binOverlapsJunctionStart <- rep( FALSE ,  nrow(auxdf) )
    binOverlapsJunctionStart[ queryHits( st ) ] <- TRUE
    
    jranges.end <- junctions
    end(jranges.end) <- start(jranges.end)
    end <- findOverlaps( type="end", bins, jranges.end )
    binOverlapsJunctionEnd <- rep( FALSE ,  nrow(auxdf) )
    binOverlapsJunctionEnd[ queryHits( end ) ] <- TRUE
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    isPlusStrand <- auxdf$strand == '+'
    isMinusStrand <- auxdf$strand == '-'
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Add information about the neighbor bins
    auxdf$featurePrev <- .shiftVector( as.character( auxdf$feature ),  1 , default = "-")
    auxdf$featureNext <- .shiftVector( as.character( auxdf$feature ), -1 , default = "-")
    neighborsAreExons <- auxdf$featurePrev == "E" & auxdf$featureNext == "E" 
    neighborsAreIntrons <- auxdf$featurePrev == "I" & auxdf$featureNext == "I" 
    neighborsAreExonAndIntron <- auxdf$featurePrev == "E" & auxdf$featureNext == "I" 
    neighborsAreIntronAndExon <- auxdf$featurePrev == "I" & auxdf$featureNext == "E" 
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Assign a splicing event type to bins and junctions.
    # If at least one of the neighbor bins is AS (i.e.), the classification is 
    # made using annotated junctions. Otherwise is made using exon or intron 
    # annotation.
    
    auxdf[ isAS & neighbourIsAS & binOverlapsJunction , 'events' ] <- "IR*"
    #auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & binOverlapsJunctionEnd, 'events' ] <- "ES*"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsExon, 'events' ] <- "ES*"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & ! binOverlapsJunctionEnd & isPlusStrand, 'events' ] <- "Alt3ss*" 
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & ! binOverlapsJunctionEnd & isMinusStrand, 'events' ] <- "Alt5ss*"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  ! binOverlapsJunctionStart & binOverlapsJunctionEnd & isPlusStrand, 'events' ] <- "Alt5ss*" 
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  ! binOverlapsJunctionStart & binOverlapsJunctionEnd & isMinusStrand, 'events' ] <- "Alt3ss*" 
    
    auxdf[ isAS & neighbourIsAS & binOverlapsJunction , 'eventsJ' ] <- "IR"
    #auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & binOverlapsJunctionEnd, 'eventsJ' ] <- "ES"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsExon, 'eventsJ' ] <- "ES"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & ! binOverlapsJunctionEnd & isPlusStrand, 'eventsJ' ] <- "Alt3ss" 
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  binOverlapsJunctionStart & ! binOverlapsJunctionEnd & isMinusStrand, 'eventsJ' ] <- "Alt5ss"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  ! binOverlapsJunctionStart & binOverlapsJunctionEnd & isPlusStrand, 'eventsJ' ] <- "Alt5ss"
    auxdf[ isAS & neighbourIsAS & ! binOverlapsJunction &  ! binOverlapsJunctionStart & binOverlapsJunctionEnd & isMinusStrand, 'eventsJ' ] <- "Alt3ss" 
    
    auxdf[ isAS & !neighbourIsAS & neighborsAreExons , 'events' ] <- "IR"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntrons , 'events' ] <- "ES"
    auxdf[ isAS & !neighbourIsAS & neighborsAreExonAndIntron & isPlusStrand, 'events' ] <- "Alt5ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreExonAndIntron & isMinusStrand, 'events' ] <- "Alt3ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntronAndExon & isPlusStrand, 'events' ] <- "Alt3ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntronAndExon & isMinusStrand, 'events' ] <- "Alt5ss"
    
    auxdf[ isAS & !neighbourIsAS & neighborsAreExons , 'eventsJ' ] <- "IR"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntrons , 'eventsJ' ] <- "ES"
    auxdf[ isAS & !neighbourIsAS & neighborsAreExonAndIntron & isPlusStrand, 'eventsJ' ] <- "Alt5ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreExonAndIntron & isMinusStrand, 'eventsJ' ] <- "Alt3ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntronAndExon & isPlusStrand, 'eventsJ' ] <- "Alt3ss"
    auxdf[ isAS & !neighbourIsAS & neighborsAreIntronAndExon & isMinusStrand, 'eventsJ' ] <- "Alt5ss"
    # -------------------------------------------------------------------------- #
    
    # -------------------------------------------------------------------------- #
    # Count events by type.
    AsNotExternal <- sum( intron.as == "Undefined AS" )
    totAS <- sum( auxdf$class == "Undefined AS" )
    
    ES     <- sum( auxdf$events == "ES" )
    IR     <- sum( auxdf$events == "IR" )
    Alt5ss <- sum( auxdf$events == "Alt5ss" )
    Alt3ss <- sum( auxdf$events == "Alt3ss" )
    
    mult       <- sum( isAS & neighbourIsAS )
    multES     <- sum( auxdf$events == "ES*" )
    multIR     <- sum( auxdf$events == "IR*" )
    multAlt5ss <- sum( auxdf$events == "Alt5ss*" )
    multAlt3ss <- sum( auxdf$events == "Alt3ss*" )
    # -------------------------------------------------------------------------- #
    
    
    # -------------------------------------------------------------------------- #
    # Append junction metadata to bins 
    mcols(bins) <- append(mcols(bins), DataFrame(event=auxdf$events))
    mcols(bins) <- append(mcols(bins), DataFrame(eventJ=auxdf$eventsJ))
    # -------------------------------------------------------------------------- #
    
    exportTextMessage <- function (totAS, AsNotExternal, ES, IR, Alt5ss, Alt3ss, 
                                   mult, multES,multIR, multAlt5ss, multAlt3ss   ) {
        
        paste0("* Number of AS bins (not include external) = ", totAS, "\n",
               "* Number of AS bins (include external) = ", AsNotExternal,"\n",
               "* Classified as: \n",
               "\tES bins = ", ES, "\t(",round(ES/totAS*100), "%)\n", 
               "\tIR bins = " , IR,"\t(",round(IR/totAS*100), "%)\n",
               "\tAlt5'ss bins = ", Alt5ss, "\t(",round(Alt5ss/totAS*100), "%)\n", 
               "\tAlt3'ss bins = ", Alt3ss,"\t(",round(Alt3ss/totAS*100), "%)\n",
               "\tMultiple AS bins = ", mult, "\t","(",round(mult/totAS*100), "%)\n", 
               "\tclassified as:\n",
               "\t\t\t", "ES bins = ", multES, "\t(",round(multES/mult*100), "%)\n", 
               "\t\t\t","IR bins = " , multIR,"\t(",round(multIR/mult*100), "%)\n",
               "\t\t\t","Alt5'ss bins = ", multAlt5ss, "\t(",round(multAlt5ss/mult*100), "%)\n", 
               "\t\t\t","Alt3'ss bins = ", multAlt3ss,"\t(",round(multAlt3ss/mult*100), "%)\n" )
    }
    
    textMsg <- exportTextMessage( totAS, AsNotExternal, ES, IR, Alt5ss, Alt3ss,
                                  mult, multES,multIR, multAlt5ss, multAlt3ss )
    message(textMsg)
    
    if ( sink.number() > 0 ) { 
        cat( textMsg )
        sink()
    }
    return(bins)
}