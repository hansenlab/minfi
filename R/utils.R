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

.show.preprocessMethod <- function(preprocessMethod) {
    if(length(preprocessMethod) == 3 && is.null(names(preprocessMethod)))
        names(preprocessMethod) <- c("rg.norm", "minfi", "manifest")
    if(length(preprocessMethod) == 0)
        preprocessMethod <- c(rg.norm = "unknown", minfi = "unknown", manifest = "unknown")
    cat("Preprocessing\n")
    cat(sprintf("  Method: %s\n  minfi version: %s\n  Manifest version: %s\n",
                preprocessMethod["rg.norm"], preprocessMethod["minfi"],
                preprocessMethod["manifest"]))
}

.getManifestString <- function(annotation) {
    if(length(annotation) == 1)
        return(paste0(annotation, "manifest"))
    if("array" %in% names(annotation))
        return(paste0(annotation["array"], "manifest"))
    stop("unable to get the manifest string for this object")
}

.getAnnotationString <- function(annotation) {
    if(length(annotation) == 1)
        return(paste0(annotation, "annotation"))
    if(all(c("array", "annotation") %in% names(annotation)))
        return(paste0(annotation["array"], "annotation.", annotation["annotation"]))
    stop("unable to get the annotation string for this object")
}

.betaFromMethUnmeth <- function(Meth, Unmeth, object, offset = 0,
                                betaThreshold = 0, minZero = TRUE) {
    stopifnot(offset >= 0)
    stopifnot(betaThreshold >= 0 & betaThreshold <= 0.5)
    if(minZero) {
        Meth <- pmax(Meth, 0)
        Unmeth <- pmax(Unmeth, 0)
    }
    beta <- Meth / (Meth + Unmeth + offset)
    if(betaThreshold > 0) {
        beta <- pmin(pmax(beta, betaThreshold), 1-betaThreshold)
    }
    beta
}

