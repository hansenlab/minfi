# Exported classes -------------------------------------------------------------

setClass(
    "GenomicMethylSet",
    representation(preprocessMethod = "character", annotation = "character"),
    contains = "RangedSummarizedExperiment"
)

# Validity methods -------------------------------------------------------------

setValidity("GenomicMethylSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Meth", "Unmeth")))
    if (class(rowRanges(object)) != "GRanges") {
        msg <- validMsg(
            msg,
            sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowRanges'",
                    class(object)))
    }
    if (is.null(msg)) TRUE else msg
})

# Exported functions -----------------------------------------------------------

GenomicMethylSet <- function(gr = GRanges(), Meth = new("matrix"),
                             Unmeth = new("matrix"), annotation = "",
                             preprocessMethod = "", ...) {
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    new("GenomicMethylSet",
        SummarizedExperiment(
            assays = assays,
            rowRanges = as(gr, "GRanges"),
            ...),
        annotation = annotation,
        preprocessMethod = preprocessMethod
    )
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "GenomicMethylSet"), function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setMethod("getMeth", signature(object = "GenomicMethylSet"), function(object) {
    assay(object, "Meth")
})

setMethod(
    "getUnmeth",
    signature(object = "GenomicMethylSet"),
    function(object) assay(object, "Unmeth")
)

setMethod(
    "getBeta",
    signature(object = "GenomicMethylSet"),
    function(object, type = "", offset = 0, betaThreshold = 0) {
        if (type == "Illumina") {
            offset <- 100
        }
        .betaFromMethUnmeth(
            Meth = getMeth(object),
            Unmeth = getUnmeth(object),
            offset = offset,
            betaThreshold = betaThreshold)
    }
)

setMethod(
    "getM",
    signature(object = "GenomicMethylSet"),
    function(object, type = "", ...) {
        if (type == "") return(log2(getMeth(object) / getUnmeth(object)))
        if (type == "beta" || type == "Beta") logit2(getBeta(object, ...))
    }
)

setMethod(
    "getCN",
    signature(object = "GenomicMethylSet"),
    function(object, ...) log2(getMeth(object) + getUnmeth(object))
)

setMethod(
    "mapToGenome",
    signature(object = "GenomicMethylSet"),
    function(object, ...) object
)

setMethod(
    "preprocessMethod",
    signature(object = "GenomicMethylSet"),
    function(object) object@preprocessMethod
)

setMethod(
    "annotation",
    signature(object = "GenomicMethylSet"),
    function(object) object@annotation
)

setReplaceMethod(
    "annotation",
    signature(object = "GenomicMethylSet"),
    function(object, value) {
        object@annotation <- value
        object
    }
)

setMethod(
    "ratioConvert",
    signature(object = "GenomicMethylSet"),
    function(object, what = c("beta", "M", "both"), keepCN = TRUE, ...) {
        what <- match.arg(what)
        if (what == "beta") {
            Beta <- getBeta(object, ...)
            M <- NULL
        }
        if (what == "M") {
            Beta <- NULL
            M <- getM(object, ...)
        }
        if (what == "both") {
            Beta <- getBeta(object, ...)
            M <- getM(object, ...)
        }
        if (keepCN) {
            CN <- getCN(object)
        } else {
            CN <- NULL
        }
        GenomicRatioSet(
            Beta = Beta,
            M = M,
            CN = CN,
            gr = granges(object),
            colData = colData(object),
            annotation = annotation(object),
            preprocessMethod = preprocessMethod(object),
            metadata = metadata(object))
    }
)

setMethod(
    "updateObject",
    signature(object = "GenomicMethylSet"),
    function(object, ..., verbose = FALSE) {
        if (object@annotation["annotation"] == "ilmn.v1.2") {
            object@annotation["annotation"] <- .default.450k.annotation
        }
        object
    }
)

setMethod(
    "combine",
    signature(x = "GenomicMethylSet", y = "GenomicMethylSet"),
    function(x, y, ...) {
        colDataFix <- .harmonizeDataFrames(.pDataFix(colData(x)),
                                           .pDataFix(colData(y)))
        colData(x) <- colDataFix$x
        colData(y) <- colDataFix$y
        cbind(x,y)
    }
)
