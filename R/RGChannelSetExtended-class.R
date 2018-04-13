# Exported classes -------------------------------------------------------------

setClass("RGChannelSetExtended", contains = "RGChannelSet")

# Validity methods -------------------------------------------------------------

setValidity("RGChannelSetExtended", function(object) {
    msg <- validMsg(
        NULL,
        .checkAssayNames(
            object, c("Red", "Green", "RedSD", "GreenSD", "NBeads")))
    if (is.null(msg)) TRUE else msg
})

# Exported functions -----------------------------------------------------------

RGChannelSetExtended <- function(Green = new("matrix"), Red = new("matrix"),
                                 GreenSD = new("matrix"), RedSD = new("matrix"),
                                 NBeads = new("matrix"), annotation = "",
                                 ...) {
    # Check rownames, colnames
    assays <- SimpleList(
        Green = Green,
        Red = Red,
        GreenSD = GreenSD,
        RedSD = RedSD,
        NBeads = NBeads)

    new("RGChannelSetExtended",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation
    )
}

getNBeads <- function(object) {
    if (!is(object, "RGChannelSetExtended")) {
        msg <- "object is of class '%s', but needs to be of class 'RGChannelSetExtended'"
        stop(sprintf(msg, class(object)))
    }
    assay(object, "NBeads")
}

# Exported methods -------------------------------------------------------------

setMethod(
    "updateObject",
    signature(object = "RGChannelSetExtended"),
    function(object, ..., verbose = FALSE) {
        if (verbose) message("updateObject(object = 'RGChannelSetExtended')")
        if ("assayData" %in% names(getObjectSlots(object))) {
            # This is an ExpressionSet based object
            object <- RGChannelSetExtended(
                Green = getObjectSlots(object)[["assayData"]][["Green"]],
                Red = getObjectSlots(object)[["assayData"]][["Red"]],
                GreenSD = getObjectSlots(object)[["assayData"]][["GreenSD"]],
                RedSD = getObjectSlots(object)[["assayData"]][["RedSD"]],
                NBeads = getObjectSlots(object)[["assayData"]][["NBeads"]],
                colData = getObjectSlots(
                    getObjectSlots(object)[["phenoData"]])[["data"]],
                annotation = getObjectSlots(object)[["annotation"]])
        }
        object
    }
)

setMethod(
    "coerce",
    signature(from = "RGChannelSetExtended", to = "RGChannelSet"),
    function(from, to){
        if (nrow(from) > 0 || ncol(from) > 0) {
            assays(from) <- assays(from)[c("Green", "Red")]
        }
        class(from) <- "RGChannelSet"
        from
    }
)
