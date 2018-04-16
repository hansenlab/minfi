# Exported generics ------------------------------------------------------------

setGeneric("getMeth", function(object) standardGeneric("getMeth"))

setGeneric("getUnmeth", function(object) standardGeneric("getUnmeth"))

setGeneric("getBeta", function(object, ...) standardGeneric("getBeta"))

setGeneric("getM", function(object, ...) standardGeneric("getM"))

setGeneric("getCN", function(object, ...) standardGeneric("getCN"))

setGeneric("mapToGenome", function(object, ...) standardGeneric("mapToGenome"))

setGeneric("getManifest", function(object) standardGeneric("getManifest"))

setGeneric(
    "ratioConvert",
    function(object, ...) standardGeneric("ratioConvert"))

setGeneric(
    "preprocessMethod",
    function(object) standardGeneric("preprocessMethod"))


setGeneric(
    "convertArray",
    function(object, ...) standardGeneric("convertArray"))

setGeneric(
    "combineArrays",
    function(object1, object2, ...) standardGeneric("combineArrays"))
