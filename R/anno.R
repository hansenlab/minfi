setClass("IlluminaMethylationManifest",
         representation(data = "environment",
                        annotation = "character"))

setValidity("IlluminaMethylationManifest", function(object) {
    msg <- NULL
    if(! "TypeI" %in% ls(object@data) || !is(object@data[["TypeI"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with TypeI probes", sep = "\n")
    if(! "TypeII" %in% ls(object@data) || !is(object@data[["TypeII"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with TypeII probes", sep = "\n")
    if(! "TypeControl" %in% ls(object@data) || !is(object@data[["TypeControl"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with Control probes", sep = "\n")
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationManifest", function(object) {
    cat("IlluminaMethylationManifest object\n")
    .show.annotation(object@annotation)
    cat("Number of type I probes:", nrow(object@data[["TypeI"]]), "\n")
    cat("Number of type II probes:", nrow(object@data[["TypeII"]]), "\n")
    cat("Number of control probes:", nrow(object@data[["TypeControl"]]), "\n")
})

IlluminaMethylationManifest <- function(TypeI = new("data.frame"), TypeII = new("data.frame"),
                        TypeControl = new("data.frame"), annotation = "") {
    data <- new.env(parent = emptyenv())
    data[["TypeI"]] <- TypeI
    data[["TypeII"]] <- TypeII
    data[["TypeControl"]] <- TypeControl
    lockEnvironment(data, bindings = TRUE)
    manifest <- new("IlluminaMethylationManifest", annotation = annotation, data = data)
    manifest
}

setClass("IlluminaMethylationAnnotation",
         representation(data = "environment",
                        annotation = "character"))

setValidity("IlluminaMethylationAnnotation", function(object) {
    msg <- NULL
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationAnnotation", function(object) {
    cat("IlluminaMethylationAnnotation object\n")
    .show.annotation(object@annotation)
})

IlluminaMethylationAnnotation <- function(listOfObjects,
                                          annotation = "") {
    data <- new.env(parent = emptyenv())
    stopifnot(all(c("Locations.hg18", "Locations.hg19") %in% names(listOfObjects)))
    for(nam in names(listOfObjects)) {
        cat(nam, "\n")
        assign(nam, listOfObjects[[nam]], envir = data)
    }
    lockEnvironment(data, bindings = TRUE)
    anno <- new("IlluminaMethylationAnnotation",
                annotation = annotation, data = data)
    anno
}

