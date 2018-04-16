utils::globalVariables(c("channel"))

logit2 <- function(x) log2(x) - log2(1 - x)

ilogit2 <- function(x) 2^x / (1 + 2^x)

.default.27k.annotation  <- "ilmn12.hg19"
.default.450k.annotation <- "ilmn12.hg19"
.default.epic.annotation <- "ilm10b2.hg19"
.metharray.types <- c("IlluminaHumanMethylation450k",
                      "IlluminaHumanMethylationEPIC",
                      "IlluminaHumanMethylation27k")

.seqnames.order.all <- c(paste0("chr", c(1:22, "X", "Y")), "multi", "unmapped")
.seqnames.order <- paste0("chr", c(1:22, "X", "Y"))

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

.show.availableAnnotation <- function(object, indent = "  ") {
    available <- .availableAnnotation(object)
    cat("Available annotation\n")
    sapply(available$names, function(xx) {
        cat(sprintf("%s%s\n", indent, xx))
    })
    cat("Defaults\n")
    sapply(available$defaults, function(xx) {
        cat(sprintf("%s%s\n", indent, xx))
    })
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
    if(length(annotation) == 1) {
        if(annotation == "Unknown")
            stop("Cannot get Manifest object for an 'Unknown' array")
        return(paste0(annotation, "manifest"))
    }
    if("array" %in% names(annotation)) {
        if(annotation["array"] == "Unknown")
            stop("Cannot get Manifest object for an 'Unknown' array")
        return(paste0(annotation["array"], "manifest"))
    }
    stop("unable to get the manifest string for this object")
}

.getAnnotationString <- function(annotation) {
    if(length(annotation) == 1) {
        if(annotation == "Unknown")
            stop("Cannot get Annotation object for an 'Unknown' array")
        return(sprintf("%sanno", annotation))
    }
    if(all(c("array", "annotation") %in% names(annotation))) {
        if(annotation["array"] == "Unknown")
            stop("Cannot get Annotation object for an 'Unknown' array")
        return(sprintf("%sanno.%s", annotation["array"], annotation["annotation"]))
    }
    stop("unable to get the annotation string for this object")
}

.betaFromMethUnmeth <- function(Meth, Unmeth, object, offset = 0,
                                betaThreshold = 0, minZero = TRUE) {
    stopifnot(offset >= 0)
    stopifnot(betaThreshold >= 0 & betaThreshold <= 0.5)
    if (minZero) {
        Meth <- pmax2(Meth, 0)
        Unmeth <- pmax2(Unmeth, 0)
    }
    beta <- Meth / (Meth + Unmeth + offset)
    if (betaThreshold > 0) {
        beta <- pmin2(pmax2(beta, betaThreshold), 1 - betaThreshold)
    }
    beta
}

.checkAssayNames <- function(object, names) {
    nms <- names(assays(object, withDimnames = FALSE))
    if(!all(names %in% nms))
        return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                       class(object), paste0(names, collapse = ", ")))
    else
        NULL
}

.digestMatrix <- function(mat, digits = 6) {
    content <- sprintf(paste0("%.", digits, "f"), mat)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", digits)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(c(content, rownames(mat), colnames(mat)))
}

.digestVector <- function(vec, digits = 6) {
    content <- sprintf(paste0("%.", digits, "f"), vec)
    ## Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", digits)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(content)
}

.isGenomicOrStop <- function(object) {
    if(!is(object, "GenomicMethylSet") && !is(object, "GenomicRatioSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'GenomicMethylSet' or 'GenomicRatioSet'",
                     class(object)))
}

.isMethylOrStop <- function(object) {
    if(!is(object, "MethylSet") && !is(object, "GenomicMethylSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'MethylSet' or 'GenomicMethylSet'",
                     class(object)))
}

.isGenomicMethylOrStop <- function(object) {
    if(!is(object, "GenomicMethylSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'GenomicMethylSet'",
                     class(object)))
}

.isMethylOrRatioOrStop <- function(object) {
    if(!is(object, "MethylSet") && !is(object, "GenomicMethylSet") &&
       !is(object, "RatioSet") && !is(object, "GenomicRatioSet"))
        stop(sprintf("object is of class '%s', but needs to be of class '[Genomic]MethylSet' or '[Genomic]RatioSet",
                     class(object)))
}

.isMethylOrRatio <- function(object) {
    (is(object, "MethylSet") || is(object, "GenomicMethylSet") ||
        is(object, "RatioSet") || is(object, "GenomicRatioSet"))
}


.isRGOrStop <- function(object) {
    if(!is(object, "RGChannelSet"))
        stop(sprintf("object is of class '%s', but needs to be of class 'RGChannelSet' or 'RGChannelSetExtended'",
                     class(object)))
}

.is27k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation27k"
}

.is450k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation450k"
}

.isEPIC <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationEPIC"
}

.harmonizeSex <- function(vector) {
    ## FIXME: not done
    stop("function not done")
    validMale <- c("M", "MALE")
    validFemale <- c("F", "FEMALE")
    ## validUnknown <- c("U", "Unknown")
    if(is.factor(vector))
        vector <- as.character(vector)
    if(!is.character(vector))
        stop("[.harmonizeSet] argument 'vector' needs to be either a character or a factor")
    vector <- toupper(vector)
    vector[vector %in% validMale] <- "M"
    vector[vector %in% validFemale] <- "F"
    if(any(! vector %in% c("M", "F")))
        stop("[.harmonizeSet] could not harmonize the vector argument to be either 'M' or 'F'")
    vector
}

.harmonizeDataFrames <- function(x, y) {
    stopifnot(is(x, "DataFrame"))
    stopifnot(is(y, "DataFrame"))
    x.only <- setdiff(names(x), names(y))
    y.only <- setdiff(names(y), names(x))
    if(length(x.only) > 0) {
        df.add <- x[1,x.only]
        is.na(df.add[1,]) <- TRUE
        y <- cbind(y, df.add)
    }
    if(length(y.only) > 0) {
        df.add <- y[1,y.only]
        is.na(df.add[1,]) <- TRUE
        x <- cbind(x, df.add)
    }
    list(x = x, y = y[, names(x)])
}

getMethSignal <- function(object, what = c("Beta", "M"), ...) {
    what <- match.arg(what)
    switch(what,
           "Beta" = getBeta(object, ...),
           "M" = getM(object, ...)
           )
}

.pDataAdd <- function(object, df) {
    stopifnot(is(df, "data.frame") || is(df, "DataFrame"))
    pD <- colData(object)
    if(any(names(df) %in% names(pD))) {
        alreadyPresent <- intersect(names(df), names(pD))
        warning(sprintf("replacing the following columns in colData(object): %s",
                        paste(alreadyPresent, collapse = ", ")))
        pD[, alreadyPresent] <- df[, alreadyPresent]
        df <- df[, ! names(df) %in% alreadyPresent]
    }
    if(ncol(df) > 0) {
        ## Work around for bug in cbind(DataFrame, DataFrame)
        rownam <- rownames(pD)
        pD <- cbind(pD, df)
        rownames(pD) <- rownam
    }
    colData(object) <- pD
    object
}

.pDataFix <- function(df) {
    characterColumns <- c("Slide", "Array", "Sample_Name", "Basename", "SampleID")
    for(col in characterColumns) {
        if(col %in% names(df))
            df[[col]] <- as.character(df[[col]])
    }
    df
}

.NA_type <- function(type) {
    c(vector(type), NA)
}
