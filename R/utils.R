logit2 <- function(x) { log2(x) - log2(1-x) }

ilogit2 <- function(x) { 2^(x) / (1+2^(x)) }

.default.450k.annotation <- "ilmn.v1.2"

.show.ExpressionSet <- function(object) {
    cat(class(object), " (storageMode: ", storageMode(object), ")\n", sep = "")
    cat("assayData:", paste(dim(object)[[1]], "features,", dim(object)[[2]], "samples"), "\n")
    cat("  element names:",
        paste(assayDataElementNames(object), collapse=", "), "\n")
    Biobase:::.showAnnotatedDataFrame(phenoData(object),
                                      labels=list(object="phenoData"))
}
    
.show.annotation <- function(annotation, indent = "  ") {
    cat("Annotation\n")
    if(length(annotation) == 1) {
        cat(sprintf("%sarray: %s\n", indent, annotation))
    } else {
        sapply(seq(along = annotation), function(ii) {
            cat(sprintf("%s%s: %s\n", indent, names(annotation)[ii], annotation[ii]))
        })
    }
}

