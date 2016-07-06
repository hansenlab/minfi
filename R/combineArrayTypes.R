combineArrayTypes <- function(rgSet1, rgSet2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC"), verbose=TRUE){
    .isRGOrStop(rgSet1)
    .isRGOrStop(rgSet2)
    outType <- match.arg(outType)
    sampleNames1 <- sampleNames(rgSet1)
    sampleNames2 <- sampleNames(rgSet2)
    if(annotation(rgSet1)["array"] == annotation(rgSet2)["array"])
        stop("The two objects 'rgSet1' and 'rgSet2' are the same array type; the function `combineArrayTypes` is for combining different types of arrays. Have a look at 'combine'") 
    if(length(intersect(sampleNames1, sampleNames2)) > 0) {
        stop("The two objects 'rgSet1' and 'rgSet2' must have different sample names.")
    }
    if(verbose) {
        message(sprintf("[combineArrayTypes] Casting as %s", annotation(rgSet1)["array"]))
    }
    if(outType == "IlluminaHumanMethylationEPIC") {
        if(.isEPIC(rgSet1) && .is450k(rgSet2)) {
            rgSet <- .combineArrayTypes_450k_epic(rgSet1 = rgSet1, rgSet2 = rgSet2, verbose = verbose)
        } else if(.isEPIC(rgSet2) && .is450k(rgSet1)) {
            rgSet <- .combineArrayTypes_450k_epic(rgSet1 = rgSet2, rgSet2 = rgSet1, verbose = verbose)
        } else {
            stop("Currently, 'combineArrayTypes' only supports combining 'IlluminaHumanMethylation450k' and 'IlluminaHumanMethylationEPIC' arrays.")
        }
    }
    if(outType == "IlluminaHumanMethylation450k") {
        if(.isEPIC(rgSet1) && .is450k(rgSet2)) {
            rgSet <- .combineArrayTypes_450k_epic(rgSet1 = rgSet2, rgSet2 = rgSet1, verbose = verbose)
        } else if(.isEPIC(rgSet2) && .is450k(rgSet1)) {
            rgSet <- .combineArrayTypes_450k_epic(rgSet1 = rgSet1, rgSet2 = rgSet2, verbose = verbose)
        } else {
            stop("Currently, 'combineArrayTypes' only supports combining 'IlluminaHumanMethylation450k' and 'IlluminaHumanMethylationEPIC' arrays.")
        }
    }
    rgSet
}

.combineArrayTypes_450k_epic <- function(rgSet1, rgSet2,
                                         verbose = verbose) {
    ## This function makes the output array equal to rgSet1
    .isRGOrStop(rgSet1)
    .isRGOrStop(rgSet2)
    stopifnot((.is450k(rgSet1) && .isEPIC(rgSet2)) || (.isEPIC(rgSet1) && .is450k(rgSet2)))
    keepAddresses <- list(I = NULL, II = NULL, SnpI = NULL,
                          SnpII = NULL, Control = NULL)
    ## Probes of Type I
    probes1 <- getProbeInfo(rgSet1, type = "I")
    probes2 <- getProbeInfo(rgSet2, type = "I")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$Color == probes2$Color))
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqB))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressA, probes2$AddressB)
    wh <- which(rownames(rgSet2) %in% names(translate))
    rownames(rgSet2)[wh] <- translate[rownames(rgSet2)[wh]]
    keepAddresses$I <- unname(translate)
        
    ## Probes of Type II
    probes1 <- getProbeInfo(rgSet1, type = "II")
    probes2 <- getProbeInfo(rgSet2, type = "II")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet2) %in% names(translate))
    rownames(rgSet2)[wh] <- translate[rownames(rgSet2)[wh]]
    keepAddresses$II <- unname(translate)
    
    ## Probes of Type SnpII
    probes1 <- getProbeInfo(rgSet1, type = "SnpI")
    probes2 <- getProbeInfo(rgSet2, type = "SnpI")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqB))
    stopifnot(all(probes1$ProbeSeqB == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- c(probes1$AddressA, probes1$AddressB)
    names(translate) <- c(probes2$AddressB, probes2$AddressA)
    wh <- which(rownames(rgSet2) %in% names(translate))
    rownames(rgSet2)[wh] <- translate[rownames(rgSet2)[wh]]
    keepAddresses$SnpI <- unname(translate)

    ## Probes of Type SnpII
    probes1 <- getProbeInfo(rgSet1, type = "SnpII")
    probes2 <- getProbeInfo(rgSet2, type = "SnpII")
    commonNames <- intersect(probes1$Name, probes2$Name)
    probes1 <- probes1[match(commonNames, probes1$Name),]
    probes2 <- probes2[match(commonNames, probes2$Name),]
    stopifnot(all(probes1$ProbeSeqA == probes2$ProbeSeqA))
    ## Translating rgSet2 addresses to rgSet1 addresses
    translate <- probes1$AddressA
    names(translate) <- probes2$AddressA
    wh <- which(rownames(rgSet2) %in% names(translate))
    rownames(rgSet2)[wh] <- translate[rownames(rgSet2)[wh]]
    keepAddresses$SnpII <- unname(translate)
    
    ## Probes of Type Control
    probes1 <- getProbeInfo(rgSet1, type = "Control")
    probes2 <- getProbeInfo(rgSet2, type = "Control")
    commonAddress <- intersect(probes1$Address, probes2$Address)
    probes1 <- probes1[match(commonAddress, probes1$Address),]
    probes2 <- probes2[match(commonAddress, probes2$Address),]
    ## Even with a common address, there are 11 probes with
    ## different Type / Color / ExtendedType
    ## Discussion with Illumina support reveals that these
    ## are actually the same probes
    ## wh <- sapply(colnames(probes1), function(nam) {
    ##     which(probes1[, nam] != probes2[, nam])
    ## })
    ## wh <- unique(unname(unlist(wh)))
    ## probes1 <- probes1[-wh,]
    ## probes2 <- probes2[-wh,]
    keepAddresses$Control <- unname(probes1$Address)
    
    keepAddresses <- do.call(c, keepAddresses)
    rgSet1 <- rgSet1[keepAddresses,]
    rgSet2 <- rgSet2[keepAddresses,]
    array1 <- annotation(rgSet1)["array"]
    array2 <- annotation(rgSet2)["array"]
    annotation(rgSet2) <- annotation(rgSet1)

    ## combine does not fill out missing stuff it seems
    rgSet <- combine(rgSet1, rgSet2)
    pData(rgSet)$ArrayTypes <- rep(c(array1, array2),
                                   times = c(ncol(rgSet1), ncol(rgSet2)))
    rgSet
}
    
