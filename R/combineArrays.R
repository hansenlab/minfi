.checkCombineAnnotation <- function(object1, object2, outType) {
    array1 <- annotation(object1)["array"]
    array2 <- annotation(object2)["array"]
    stopifnot(array1 %in% .metharray.types)
    stopifnot(array2 %in% .metharray.types)
    if(! array1 %in% outType && ! array2 %in% outType) {
        stop("`combineArrayTypes` requires that one of the two arrays being combined is of the same kind as `outType`")
    }
    if(array1 == array2) {
        stop("The two objects 'object1' and 'object2' are the same array type; the function `combineArrayTypes` is for combining different types of arrays. Have a look at 'cbind'")
    }
    if(array1 == outType) {
        outAnno <- annotation(object1)
    } else if(array2 == outType) {
        outAnno <- annotation(object2)
    }
    list(annotation = outAnno)
}

.getAnnotationFromOutType <- function(outType = c("IlluminaHumanMethylation450k", 
        "IlluminaHumanMethylationEPIC", 
        "IlluminaHumanMethylation27k")){
    outType <- match.arg(outType)
    array <- outType
    if (outType=="IlluminaHumanMethylation450k"){
        anno <-  .default.450k.annotation
    } else if (outType=="IlluminaHumanMethylation27k"){
        anno <-  .default.27k.annotation
    } else {
        anno <-  .default.epic.annotation
    }
    c(array =  outType, annotation = anno)
}


.getProbesFromOutType <- function(outType = c("IlluminaHumanMethylation450k", 
        "IlluminaHumanMethylationEPIC", 
        "IlluminaHumanMethylation27k")){
    outType <- match.arg(outType )
    manifest <- getManifest(outType)
    probesI  <- getProbeInfo(manifest, type="I")$Name
    probesII  <- getProbeInfo(manifest, type="II")$Name
    probes <- c(probesI, probesII)
    probes
}

# Combine array functions:
setMethod("combineArrays",
          signature(object1 = "RGChannelSet", object2 = "RGChannelSet"),
          function(object1, object2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC"), verbose = TRUE) {
    outType <- match.arg(outType)
    array1 <- annotation(object1)[["array"]]
    array2 <- annotation(object2)[["array"]]
    if (array1=="IlluminaHumanMethylation27k" | array2=="IlluminaHumanMethylation27k"){
        stop("27k arrays cannot be combined at the RGChannelSet level.")
    }
    object1 <- convertArray(object1, outType = outType)
    object2 <- convertArray(object2, outType = outType)
    features1 <- rownames(object1)
    features2 <- rownames(object2)
    features  <- intersect(features1, features2)
    object1  <- object1[features,]
    object2  <- object2[features,]
    rgSet <- combine(object1, object2)
    rgSet
})


setMethod("combineArrays",
          signature(object1 = "MethylSet", object2 = "MethylSet"),
          function(object1, object2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"), verbose = TRUE) {
    outType <- match.arg(outType)
    object1 <- convertArray(object1, outType = outType)
    object2 <- convertArray(object2, outType = outType)
    common.features <- intersect(rownames(object1), rownames(object2))
    object1 <- object1[common.features,]
    object2 <- object2[common.features,]
    Mset <- combine(object1, object2)
    Mset
})


setMethod("combineArrays",
          signature(object1 = "RatioSet", object2 = "RatioSet"),
          function(object1, object2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"), verbose = TRUE) {
    outType <- match.arg(outType)
    object1 <- convertArray(object1, outType = outType)
    object2 <- convertArray(object2, outType = outType)
    common.features <- intersect(rownames(object1), rownames(object2))
    object1 <- object1[common.features,]
    object2 <- object2[common.features,]
    Rset <- combine(object1, object2)
    Rset
})


setMethod("combineArrays",
          signature(object1 = "GenomicRatioSet", object2 = "GenomicRatioSet"),
          function(object1, object2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"), verbose = TRUE) {
    outType <- match.arg(outType)
    object1 <- convertArray(object1, outType = outType)
    object2 <- convertArray(object2, outType = outType)

    colData1$ArrayTypes <- annotation(object1)["array"]
    colData2$ArrayTypes <- annotation(object2)["array"]
    colData1 <- colData(object1)
    colData2 <- colData(object2)
    by <- c("row.names", intersect(names(colData1), names(colData2)))
    colData.merged <- merge(colData1, colData2, all = TRUE, by = by)
    colData(object1) <- colData.merged[match(colnames(object1), colData.merged[, "Row.names"]), ]
    colData(object2) <- colData.merged[match(colnames(object2), colData.merged[, "Row.names"]), ]
    gr.common <- intersect(granges(object1), granges(object2))
    object1 <- sort(subsetByOverlaps(object1, gr.common))
    object2 <- sort(subsetByOverlaps(object2, gr.common))
    GRset <- cbind(object1, object2)
    colnames(GRset) <- GRset$Row.names
    GRset$Row.names <- NULL
    GRset
})


setMethod("combineArrays",
          signature(object1 = "GenomicMethylSet", object2 = "GenomicMethylSet"),
          function(object1, object2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"), verbose = TRUE) {
    outType <- match.arg(outType)
    object1 <- convertArray(object1, outType = outType)
    object2 <- convertArray(object2, outType = outType)
    
    colData1 <- colData(object1)
    colData2 <- colData(object2)
    colData1$ArrayTypes <- annotation(object1)["array"]
    colData2$ArrayTypes <- annotation(object2)["array"]
    by <- c("row.names", intersect(names(colData1), names(colData2)))
    colData.merged <- merge(colData1, colData2, all = TRUE, by = by)
    colData(object1) <- colData.merged[match(colnames(object1), colData.merged[, "Row.names"]), ]
    colData(object2) <- colData.merged[match(colnames(object2), colData.merged[, "Row.names"]), ]
    gr.common <- intersect(granges(object1), granges(object2))
    object1 <- sort(subsetByOverlaps(object1, gr.common))
    object2 <- sort(subsetByOverlaps(object2, gr.common))
    GMset <- cbind(object1, object2)
    colnames(GMset) <- GMset$Row.names
    GMset$Row.names <- NULL
    GMset
})

setMethod("convertArray",
          signature(object = "RGChannelSet"),
          function(object, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC"), verbose = TRUE) {
    outType <- match.arg(outType)
    array <- annotation(object)[["array"]]
    if(array == outType){
        rgSet <- object
        #stop("The 'object' is already in the 'outType' array type and the function `convertArrayTypes` is for converting to a different types of arrays.") 
    } else {
        if(verbose) {
            message(sprintf("[convertArray] Casting as %s", outType))
        }
        rgSet <- .convertArray_450k_epic(rgSet = object, outType = outType, verbose = verbose)
    }
    rgSet
})


setMethod("convertArray",
        signature(object = "MethylSet"),
        function(object, outType = c("IlluminaHumanMethylation450k", 
            "IlluminaHumanMethylationEPIC", 
            "IlluminaHumanMethylation27k"), verbose = TRUE){
    outType <- match.arg(outType)
    array <- annotation(object)[["array"]]
    if (array != outType){
        common.features <- intersect(rownames(object), .getProbesFromOutType(outType))
        object <- object[common.features,]
        annotation(object) <- .getAnnotationFromOutType(outType)
    }
    object
})


setMethod("convertArray",
        signature(object = "RatioSet"),
        function(object, outType = c("IlluminaHumanMethylation450k", 
            "IlluminaHumanMethylationEPIC", 
            "IlluminaHumanMethylation27k"), verbose = TRUE){
    outType <- match.arg(outType)
    array <- annotation(object)[["array"]]
    if (array != outType){
        common.features <- intersect(rownames(object), .getProbesFromOutType(outType))
        object <- object[common.features,]
        annotation(object) <- .getAnnotationFromOutType(outType)
    }
    object
})


setMethod("convertArray",
        signature(object = "GenomicMethylSet"),
        function(object, outType = c("IlluminaHumanMethylation450k", 
            "IlluminaHumanMethylationEPIC", 
            "IlluminaHumanMethylation27k"), verbose = TRUE){
    outType <- match.arg(outType)
    array <- annotation(object)[["array"]]
    if (array != outType){
        common.features <- intersect(rownames(object), .getProbesFromOutType(outType))
        object <- object[common.features,]
        object@annotation <- .getAnnotationFromOutType(outType)
    }
    object
})


setMethod("convertArray",
        signature(object = "GenomicRatioSet"),
        function(object, outType = c("IlluminaHumanMethylation450k",
         "IlluminaHumanMethylationEPIC", 
         "IlluminaHumanMethylation27k"), verbose = TRUE){
    outType <- match.arg(outType)
    array <- annotation(object)[["array"]]
    if (array != outType){
        common.features <- intersect(rownames(object), .getProbesFromOutType(outType))
        object <- object[common.features,]
        object@annotation <- .getAnnotationFromOutType(outType)
    }
    object
})





# Convert the rgSet into the outType array. 
.convertArray_450k_epic <- function(rgSet, outType=c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC"),
                                         verbose = verbose) {
    outType <- match.arg(outType)
    .isRGOrStop(rgSet)
    stopifnot(.is450k(rgSet) | .isEPIC(rgSet))
    
    array <- annotation(rgSet)["array"]
    if (array==outType){
        stop("'rgSet' already in the 'outType' array type.")
    }
    manifest1 <- getManifest(outType)
    manifest2 <- getManifest(rgSet)

    keepAddresses <- list(I = NULL, II = NULL, SnpI = NULL,
                          SnpII = NULL, Control = NULL)
    ## Probes of Type I
    probes1 <- getProbeInfo(manifest1, type = "I")
    probes2 <- getProbeInfo(manifest2, type = "I")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$Color == probes2$Color))
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqB))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressA, probes2$AddressB)
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$I <- unname(translate)
        
    ## Probes of Type II
    probes1 <- getProbeInfo(manifest1, type = "II")
    probes2 <- getProbeInfo(manifest2, type = "II")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$II <- unname(translate)
    
    ## Probes of Type SnpII
    probes1 <- getProbeInfo(manifest1, type = "SnpI")
    probes2 <- getProbeInfo(manifest2, type = "SnpI")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqB))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressB, probes2$AddressA)
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$SnpI <- unname(translate)

    ## Probes of Type SnpII
    probes1 <- getProbeInfo(manifest1, type = "SnpII")
    probes2 <- getProbeInfo(manifest2, type = "SnpII")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet) %in% names(translate))
    rownames(rgSet)[wh] <- translate[rownames(rgSet)[wh]]
    keepAddresses$SnpII <- unname(translate)
    
    ## Probes of Type Control
    probes1 <- getProbeInfo(manifest1, type = "Control")
    probes2 <- getProbeInfo(manifest2, type = "Control")
    commonAddress <- intersect(probes1$Address, probes2$Address)
    probes1 <- probes1[match(commonAddress, probes1$Address),]
    probes2 <- probes2[match(commonAddress, probes2$Address),]
    keepAddresses$Control <- unname(probes1$Address)


    keepAddresses <- do.call("c", keepAddresses)
    keepAddresses <- keepAddresses[keepAddresses %in% rownames(rgSet)]
    rgSet  <- rgSet[keepAddresses,]
    annotation(rgSet) <- .getAnnotationFromOutType(outType)
    rgSet
}
