# Internal functions -----------------------------------------------------------

.getAnnotationFromOutType <- function(outType = c(
    "IlluminaHumanMethylation450k",
    "IlluminaHumanMethylationEPIC",
    "IlluminaHumanMethylation27k")) {
    outType <- match.arg(outType)
    if (outType == "IlluminaHumanMethylation450k") {
        anno <-  .default.450k.annotation
    } else if (outType == "IlluminaHumanMethylation27k") {
        anno <- .default.27k.annotation
    } else {
        anno <- .default.epic.annotation
    }
    c(array = outType, annotation = anno)
}

.getLociFromOutType <- function(outType = c("IlluminaHumanMethylation450k",
                                            "IlluminaHumanMethylationEPIC",
                                            "IlluminaHumanMethylation27k")) {
    outType <- match.arg(outType)
    manifest <- getManifest(outType)
    probesI <- getProbeInfo(manifest, type = "I")$Name
    probesII <- getProbeInfo(manifest, type = "II")$Name
    probes <- c(probesI, probesII)
    probes
}

# Convert the rgSet into the outType array.
.convertArray_450k_epic <- function(rgSet,
                                    outType = c("IlluminaHumanMethylation450k",
                                                "IlluminaHumanMethylationEPIC"),
                                    verbose = verbose) {

    outType <- match.arg(outType)
    .isRGOrStop(rgSet)
    stopifnot(.is450k(rgSet) || .isEPIC(rgSet))

    array <- annotation(rgSet)["array"]
    if (array == outType) stop("'rgSet' already in the 'outType' array type.")
    manifest1 <- getManifest(outType)
    manifest2 <- getManifest(rgSet)

    keepAddresses <- list(
        I = NULL,
        II = NULL,
        SnpI = NULL,
        SnpII = NULL,
        Control = NULL)

    # Probes of Type I
    probes1 <- getProbeInfo(manifest1, type = "I")
    probes2 <- getProbeInfo(manifest2, type = "I")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name), ]
    probes2 <- probes2[match(commonNames, probes2$Name), ]
    stopifnot(all(probes1$Color == probes2$Color))
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqB))
    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressA, probes2$AddressB)
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$I <- unname(translate)

    # Probes of Type II
    probes1 <- getProbeInfo(manifest1, type = "II")
    probes2 <- getProbeInfo(manifest2, type = "II")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$II <- unname(translate)

    # Probes of Type SnpI
    probes1 <- getProbeInfo(manifest1, type = "SnpI")
    probes2 <- getProbeInfo(manifest2, type = "SnpI")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqB))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqA))
    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressA, probes2$AddressB)
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$SnpI <- unname(translate)

    # Probes of Type SnpII
    probes1 <- getProbeInfo(manifest1, type = "SnpII")
    probes2 <- getProbeInfo(manifest2, type = "SnpII")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    # Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$SnpII <- unname(translate)

    # Probes of Type Control
    probes1 <- getProbeInfo(manifest1, type = "Control")
    probes2 <- getProbeInfo(manifest2, type = "Control")
    commonAddress <- intersect(probes1$Address, probes2$Address)
    probes1 <- probes1[match(commonAddress, probes1$Address),]
    probes2 <- probes2[match(commonAddress, probes2$Address),]
    keepAddresses$Control <- unname(probes1$Address)

    # Update rgSet
    keepAddresses <- do.call("c", keepAddresses)
    keepAddresses <- keepAddresses[keepAddresses %in% rownames(rgSet)]
    rgSet  <- rgSet[keepAddresses, ]
    annotation(rgSet) <- .getAnnotationFromOutType(outType)
    rgSet
}

# Exported methods -------------------------------------------------------------

setMethod(
    "combineArrays",
    signature(object1 = "RGChannelSet", object2 = "RGChannelSet"),
    function(object1, object2,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array1 <- annotation(object1)[["array"]]
        array2 <- annotation(object2)[["array"]]
        if (array1 == array2) outType <- array1
        if (array1 == "IlluminaHumanMethylation27k" ||
            array2 == "IlluminaHumanMethylation27k") {
            stop("27k arrays cannot be combined at the RGChannelSet level.")
        }
        object1 <- convertArray(object1, outType = outType, verbose = verbose)
        object2 <- convertArray(object2, outType = outType, verbose = verbose)
        features1 <- rownames(object1)
        features2 <- rownames(object2)
        features  <- intersect(features1, features2)
        object1  <- object1[features, ]
        object2  <- object2[features, ]
        rgSet <- combine(object1, object2)
        rgSet$ArrayTypes <- rep(
            x = c(array1, array2),
            times = c(ncol(object1), ncol(object2)))
        rgSet
    }
)

setMethod("combineArrays",
          signature(object1 = "MethylSet", object2 = "MethylSet"),
          function(object1, object2,
                   outType = c("IlluminaHumanMethylation450k",
                               "IlluminaHumanMethylationEPIC",
                               "IlluminaHumanMethylation27k"),
                   verbose = TRUE) {

              outType <- match.arg(outType)
              array1 <- annotation(object1)["array"]
              array2 <- annotation(object2)["array"]
              if (array1 == array2) outType <- array1
              object1 <- convertArray(
                  object = object1,
                  outType = outType,
                  verbose = verbose)
              object2 <- convertArray(
                  object = object2,
                  outType = outType,
                  verbose = verbose)
              common.features <- intersect(rownames(object1), rownames(object2))
              object1 <- object1[common.features, ]
              object2 <- object2[common.features, ]
              Mset <- combine(object1, object2)
              Mset$ArrayTypes <- rep(
                  x = c(array1, array2),
                  times = c(ncol(object1), ncol(object2)))
              Mset
          }
)

setMethod("combineArrays",
          signature(object1 = "RatioSet", object2 = "RatioSet"),
          function(object1, object2,
                   outType = c("IlluminaHumanMethylation450k",
                               "IlluminaHumanMethylationEPIC",
                               "IlluminaHumanMethylation27k"),
                   verbose = TRUE) {

              outType <- match.arg(outType)
              array1 <- annotation(object1)["array"]
              array2 <- annotation(object2)["array"]
              if (array1 == array2) outType <- array1
              object1 <- convertArray(
                  object = object1,
                  outType = outType,
                  verbose = verbose)
              object2 <- convertArray(
                  object = object2,
                  outType = outType,
                  verbose = verbose)
              common.features <- intersect(rownames(object1), rownames(object2))
              object1 <- object1[common.features, ]
              object2 <- object2[common.features, ]
              Rset <- combine(object1, object2)
              Rset$ArrayTypes <- rep(
                  x = c(array1, array2),
                  times = c(ncol(object1), ncol(object2)))
              Rset
          }
)

setMethod(
    "combineArrays",
    signature(object1 = "GenomicRatioSet", object2 = "GenomicRatioSet"),
    function(object1, object2,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array1 <- annotation(object1)["array"]
        array2 <- annotation(object2)["array"]
        if (array1 == array2) outType <- array1
        object1 <- convertArray(object1, outType = outType, verbose = verbose)
        object2 <- convertArray(object2, outType = outType, verbose = verbose)
        colData1 <- colData(object1)
        colData2 <- colData(object2)
        colData1$ArrayTypes <- array1
        colData2$ArrayTypes <- array2
        colData1 <- colData(object1)
        colData2 <- colData(object2)
        by <- c("row.names", intersect(names(colData1), names(colData2)))
        colData.merged <- merge(colData1, colData2, all = TRUE, by = by)
        colData(object1) <- colData.merged[match(
            x = colnames(object1),
            table = colData.merged[, "Row.names"]), ]
        colData(object2) <- colData.merged[match(
            x = colnames(object2),
            table = colData.merged[, "Row.names"]), ]
        gr.common <- intersect(granges(object1), granges(object2))
        object1 <- sort(subsetByOverlaps(object1, gr.common))
        object2 <- sort(subsetByOverlaps(object2, gr.common))
        GRset <- cbind(object1, object2)
        colnames(GRset) <- GRset$Row.names
        GRset$Row.names <- NULL
        GRset
    }
)

setMethod(
    "combineArrays",
    signature(object1 = "GenomicMethylSet", object2 = "GenomicMethylSet"),
    function(object1, object2,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array1 <- annotation(object1)["array"]
        array2 <- annotation(object2)["array"]
        if (array1 == array2) outType <- array1
        object1 <- convertArray(object1, outType = outType, verbose = verbose)
        object2 <- convertArray(object2, outType = outType, verbose = verbose)
        colData1 <- colData(object1)
        colData2 <- colData(object2)
        colData1$ArrayTypes <- array1
        colData2$ArrayTypes <- array2
        by <- c("row.names", intersect(names(colData1), names(colData2)))
        colData.merged <- merge(colData1, colData2, all = TRUE, by = by)
        colData(object1) <- colData.merged[match(
            x = colnames(object1),
            table = colData.merged[, "Row.names"]), ]
        colData(object2) <- colData.merged[match(
            x = colnames(object2),
            table = colData.merged[, "Row.names"]), ]
        gr.common <- intersect(granges(object1), granges(object2))
        object1 <- sort(subsetByOverlaps(object1, gr.common))
        object2 <- sort(subsetByOverlaps(object2, gr.common))
        GMset <- cbind(object1, object2)
        colnames(GMset) <- GMset$Row.names
        GMset$Row.names <- NULL
        GMset
    }
)

setMethod(
    "convertArray",
    signature(object = "RGChannelSet"),
    function(object,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array <- annotation(object)[["array"]]
        if (array == outType) return(object)
        if (verbose) message(sprintf("[convertArray] Casting as %s", outType))
        .convertArray_450k_epic(
            rgSet = object,
            outType = outType,
            verbose = verbose)
    }
)

setMethod(
    "convertArray",
    signature(object = "MethylSet"),
    function(object,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array <- annotation(object)[["array"]]
        if (array == outType) return(object)
        if (verbose) message(sprintf("[convertArray] Casting as %s", outType))
        common.features <- intersect(
            x = rownames(object),
            y = .getLociFromOutType(outType))
        object <- object[common.features, ]
        annotation(object) <- .getAnnotationFromOutType(outType)
        object
    }
)

setMethod(
    "convertArray",
    signature(object = "RatioSet"),
    function(object,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array <- annotation(object)[["array"]]
        if (array == outType) return(object)
        if (verbose) message(sprintf("[convertArray] Casting as %s", outType))
        common.features <- intersect(
            x = rownames(object),
            y = .getLociFromOutType(outType))
        object <- object[common.features, ]
        annotation(object) <- .getAnnotationFromOutType(outType)
        object
    }
)

setMethod(
    "convertArray",
    signature(object = "GenomicMethylSet"),
    function(object,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {


        outType <- match.arg(outType)
        array <- annotation(object)[["array"]]
        if (array == outType) return(object)
        if (verbose) message(sprintf("[convertArray] Casting as %s", outType))
        common.features <- intersect(
            x = rownames(object),
            y = .getLociFromOutType(outType))
        object <- object[common.features, ]
        object@annotation <- .getAnnotationFromOutType(outType)
        object
    }
)

setMethod(
    "convertArray",
    signature(object = "GenomicRatioSet"),
    function(object,
             outType = c("IlluminaHumanMethylation450k",
                         "IlluminaHumanMethylationEPIC",
                         "IlluminaHumanMethylation27k"),
             verbose = TRUE) {

        outType <- match.arg(outType)
        array <- annotation(object)[["array"]]
        if (array == outType) return(object)
        if (verbose) message(sprintf("[convertArray] Casting as %s", outType))
        common.features <- intersect(
            x = rownames(object),
            y = .getLociFromOutType(outType))
        object <- object[common.features, ]
        object@annotation <- .getAnnotationFromOutType(outType)
        object
    }
)

# TODOs ------------------------------------------------------------------------

# TODO: Lots of duplicated code; DRY
