# TODO: type() for all RealizationSink subclasses
setMethod("type", "HDF5RealizationSink", function(x) {
    x@type
})
setMethod("type", "arrayRealizationSink", function(x) {
    DelayedArray::type(x@result_envir$result)
})
setMethod("type", "RleRealizationSink", function(x) {
    x@type
})
# TODO: dimnames() for all RealizationSink subclasses
setMethod("dimnames", "arrayRealizationSink", function(x) {
    dimnames(x@result_envir$result)
})

# NOTE: DelayedArray::blockApply() with the option to write the blocks to
#       'sink'. Useful, for example, to apply a function across column-blocks
#       of a DelayedMatrix and write these results to disk and then wrap
#       these in a DelayedMatrix.
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/10
blockApplyWithRealization <- function(x, FUN, ..., grid = NULL, sink = NULL,
                                      BPREDO = list(), BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)
    grid <- DelayedArray:::.normarg_grid(grid, x)
    nblock <- length(grid)
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        viewport <- grid[[b]]
        block <- DelayedArray:::extract_block(x, viewport)
        if (!is.array(block)) {
            block <- DelayedArray:::.as_array_or_matrix(block)
        }
        attr(block, "from_grid") <- grid
        attr(block, "block_id") <- b
        block_ans <- FUN(block, ...)
        # NOTE: This is the only part different from DelayedArray::blockApply()
        if (!is.null(sink)) {
            write_block_to_sink(block_ans, sink, viewport)
            block_ans <- NULL
        }
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}

# NOTE: A mapply()-like function for conformable arrays.
# NOTE: Different from DelayedArray:::block_Mapply(); designed to have an API more like
#       DelayedArray::blockArray()
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/11
blockMapply <- function(FUN, ..., grids = NULL, BPREDO = list(),
                        BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)
    dots <- unname(list(...))
    dims <- lapply(dots, dim)
    if (!all(vapply(dims, function(dim) all(dim == dims[[1L]]), logical(1L)))) {
        stop("non-conformable arrays")
    }
    if (is.null(grids)) {
        grids <- replicate(length(dots), NULL)
    }
    grids <- mapply(DelayedArray:::.normarg_grid, grids, dots, SIMPLIFY = FALSE,
                    USE.NAMES = FALSE)
    grid_dims <- lapply(grids, dim)
    if (!all(vapply(grid_dims, function(dim) all(dim == grid_dims[[1L]]), logical(1L)))) {
        stop("non-conformable grids")
    }
    stopifnot(length(dots) == length(grids))
    nblock <- length(grids[[1]])
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        viewports <- lapply(grids, function(grid) grid[[b]])
        blocks <- mapply(function(x, viewport) {
            block <- DelayedArray:::extract_block(x, viewport)
            if (!is.array(block)) {
                block <- DelayedArray:::.as_array_or_matrix(block)
            }
            block
        }, x = dots, viewport = viewports, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        block_ans <- do.call(FUN, blocks)
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
        block_ans
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}
