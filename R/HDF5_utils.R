#-------------------------------------------------------------------------------
# Is a DelayedMatrix object (or the assays of a SummarizedExperiment object)
# backed by a HDF5 file?
#

.getSeedClasses <- function(seed) {
    if (is(seed, "SeedBinder") || is(seed, "ConformableSeedCombiner")) {
        seeds <- seed@seeds
        return(lapply(seeds, .getSeedClasses))
    } else if (is(seed, "DelayedArray")) {
        # A DelayedArray can have another DelayedArray as a seed
        seed <- seed(seed)
        return(.getSeedClasses(seed))
    }
    else {
        class(seed)
    }
}

# NOTE: Returns TRUE if *any* assay is HDF5Array-backed and FALSE if *all*
#       assays are not HDF5Array-backed
.isHDF5ArrayBacked <- function(object) {
    if (is(object, "SummarizedExperiment")) {
        return(all(vapply(X = assays(object, withDimnames = FALSE),
                          FUN = .isHDF5ArrayBacked,
                          FUN.VALUE = logical(1L))))
    }
    if (is(object, "DelayedArray")) {
        seed <- seed(object)
        seed_classes <- .getSeedClasses(seed)
        is_hdf5_backed <- vapply(unlist(seed_classes, use.names = FALSE),
                                 extends, class2 = "HDF5ArraySeed",
                                 logical(1L))
        return(any(is_hdf5_backed))
    } else if (is.matrix(object)) {
        FALSE
    } else if (is.null(object)) {
        FALSE
    } else {
        stop("Don't know how to handle object of class ", class(object))
    }
}
