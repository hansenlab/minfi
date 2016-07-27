preprocessRaw <- function(rgSet) {
    .isRGOrStop(rgSet)
    locusNames <- getManifestInfo(rgSet, "locusNames")
    M <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                dimnames = list(locusNames, sampleNames(rgSet)))
    U <- M
    TypeII.Name <- getProbeInfo(rgSet, type = "II")$Name
    M[TypeII.Name,] <- getGreen(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
    U[TypeII.Name,] <- getRed(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    M[TypeI.Red$Name,] <- getRed(rgSet)[TypeI.Red$AddressB,]
    M[TypeI.Green$Name,] <- getGreen(rgSet)[TypeI.Green$AddressB,]
    U[TypeI.Red$Name,] <- getRed(rgSet)[TypeI.Red$AddressA,]
    U[TypeI.Green$Name,] <- getGreen(rgSet)[TypeI.Green$AddressA,]
    out <- MethylSet(Meth = M, Unmeth = U, phenoData = phenoData(rgSet),
                     annotation = annotation(rgSet))
    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.          
    ## packageVersion expects a string
    out@preprocessMethod <- c(rg.norm = "Raw (no normalization or bg correction)",
                              minfi = as.character(packageVersion("minfi")),
                              manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
                              #packageVersion(getManifest(rgSet)))
    out
}


normalize.illumina.control <- function(rgSet, reference=1) {
    ## This function returns an rgset, not a methylset
    ## code duplication
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)    
    
    if (.is450k(rgSet) || .isEPIC(rgSet)) {
        AT.controls <- getControlAddress(rgSet, controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(rgSet, controlType = c("NORM_C", "NORM_G"))
    }
    if (.is27k(rgSet)) {
        AT.controls <- getControlAddress(rgSet, controlType = "Normalization-Red")
        CG.controls <- getControlAddress(rgSet, controlType = "Normalization-Green")
    }
    Green.avg <- colMeans(Green[CG.controls, , drop = FALSE])
    Red.avg <- colMeans(Red[AT.controls, , drop = FALSE])
    ref <- (Green.avg + Red.avg)[reference]/2
    if(is.na(ref))
        stop("perhaps 'reference' refer to an array that is not present.")
    Green.factor <- ref/Green.avg
    Red.factor <- ref/Red.avg
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assayDataElement(rgSet, "Green") <- Green
    assayDataElement(rgSet, "Red") <- Red
    rgSet
}

bgcorrect.illumina <- function(rgSet) {
    .isRGOrStop(rgSet)
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)
    if (.is450k(rgSet) || .isEPIC(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
    }
    if (.is27k(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "Negative")
    }
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) sort(xx)[31])
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) sort(xx)[31])
    Green <- pmax(sweep(Green, 2, Green.bg), 0)
    Red <- pmax(sweep(Red, 2, Red.bg), 0)
    assayDataElement(rgSet, "Green") <- Green
    assayDataElement(rgSet, "Red") <- Red
    rgSet
}


preprocessIllumina <- function(rgSet, bg.correct = TRUE, normalize = c("controls", "no"),
                                reference = 1) {
    .isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if(normalize == "controls") {
        rgSet <- normalize.illumina.control(rgSet, reference = reference)
    }
    if(bg.correct) {
        rgSet <- bgcorrect.illumina(rgSet)
    }
    out <- preprocessRaw(rgSet)
    preprocess <- sprintf("Illumina, bg.correct = %s, normalize = %s, reference = %d",
                          bg.correct, normalize, reference)
    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.          
    ## packageVersion expects a string
    out@preprocessMethod <- c(rg.norm = preprocess,
                              minfi = as.character(packageVersion("minfi")),
                              manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
                              #packageVersion(getManifest(rgSet)))
    out
}


detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, sampleNames(rgSet)))

    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")   
    r <- getRed(rgSet)
    rBg <- r[controlIdx,,drop=FALSE]
    rMu <- matrixStats::colMedians(rBg)
    rSd <- matrixStats::colMads(rBg)

    g <- getGreen(rgSet)
    gBg <- g[controlIdx,,drop=FALSE]
    gMu <- matrixStats::colMedians(gBg)
    gSd <- matrixStats::colMads(gBg)

    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    for (i in 1:ncol(rgSet)) {   
        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
        detP[TypeI.Red$Name, i] <- 1-pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*2)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
        detP[TypeI.Green$Name, i] <- 1-pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*2)
        ## Type II
        intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
        detP[TypeII$Name, i] <- 1-pnorm(intensity, mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
    }
    detP
}
