

setGeneric("getDesign", function(object) standardGeneric("getDesign"))
setGeneric("getSizeFactors", function(object) standardGeneric("getSizeFactors"))
setGeneric("counts", function(object) standardGeneric("counts"))
setGeneric("metadata", function(object) standardGeneric("metadata"))



setMethod("getDesign",
          signature = "SummarizedExperiment",
          definition = function(object) object@metadata$design
)

setMethod("getSizeFactors",
          signature = "SummarizedExperiment",
          definition = function(object) object@colData$sizeFactor
)

setMethod("counts",
          signature = "SummarizedExperiment",
          definition = function(object) assays(object)[["counts"]]
)

setMethod("metadata",
          signature = "SummarizedExperiment",
          definition = function(object) object@metadata
)
