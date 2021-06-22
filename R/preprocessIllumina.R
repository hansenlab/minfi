# Internal functions -----------------------------------------------------------

normalize.illumina.control <- function(rgSet, reference = 1) {
    # This function returns an rgset, not a methylset
    # code duplication
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)

    if (.is450k(rgSet) || .isEPIC(rgSet) || .isAllergy(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_C", "NORM_G"))
    }
    if (.is27k(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Red")
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Green")
    }

    Green.avg <- colMeans2(Green, rows = match(CG.controls, rownames(Green)))
    Red.avg <- colMeans2(Red, rows = match(AT.controls, rownames(Red)))
    ref <- (Green.avg + Red.avg)[reference] / 2
    if (is.na(ref)) {
        stop("perhaps 'reference' refer to an array that is not present.")
    }
    Green.factor <- ref / Green.avg
    Red.factor <- ref / Red.avg
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

bgcorrect.illumina <- function(rgSet) {
    .isRGOrStop(rgSet)
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)
    if (.is450k(rgSet) || .isEPIC(rgSet) || .isAllergy(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
    }
    if (.is27k(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "Negative")
    }
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Green <- pmax2(sweep(Green, 2, Green.bg), 0)
    Red <- pmax2(sweep(Red, 2, Red.bg), 0)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

# Exported functions -----------------------------------------------------------

# TODO: Document: This does not realize the result for a DelayedMatrix-backed
#       RGChannelSet{Extended}
preprocessIllumina <- function(rgSet, bg.correct = TRUE,
                               normalize = c("controls", "no"), reference = 1) {
    .isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if (normalize == "controls") {
        rgSet <- normalize.illumina.control(rgSet, reference = reference)
    }
    if (bg.correct) {
        rgSet <- bgcorrect.illumina(rgSet)
    }
    out <- preprocessRaw(rgSet)
    preprocess <- sprintf(
        "Illumina, bg.correct = %s, normalize = %s, reference = %d",
        bg.correct, normalize, reference)

    # TODO: The manifest package version is currently not updated since
    #       `packageVersion(getManifest(rgSet))` fails. `packageVersion()`
    #       expects a string
    out@preprocessMethod <- c(
        rg.norm = preprocess,
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}
