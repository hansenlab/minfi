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
        block <- extract_block(x, viewport)
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

.sweep_DelayedArray <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
    stopifnot(is(x, "DelayedArray"))
    dims <- dim(x)
    if (length(dims) > 2) {
        warning("Not yet implemented. Coercing 'x' to an array")
        return(sweep(as.array(x), MARGIN, STATS, FUN, check.margin, ...))
    }

    FUN <- match.fun(FUN)
    # TODO: Allow check.margin to be FALSE
    if (!check.margin) {
        stop("'check.margin' must be TRUE")
    }
    if (check.margin) {
        dimmargin <- dims[MARGIN]
        dimstats <- dim(STATS)
        lstats <- length(STATS)
        if (lstats > prod(dimmargin)) {
            warning("STATS is longer than the extent of 'dim(x)[MARGIN]'")
        } else if (is.null(dimstats)) {
            cumDim <- c(1L, cumprod(dimmargin))
            upper <- min(cumDim[cumDim >= lstats])
            lower <- max(cumDim[cumDim <= lstats])
            if (lstats && (upper %% lstats != 0L || lstats %% lower != 0L)) {
                # TODO: Make this a warning with appropriate consequent behaviour
                stop("STATS does not recycle exactly across MARGIN")
            }
        } else {
            dimmargin <- dimmargin[dimmargin > 1L]
            dimstats <- dimstats[dimstats > 1L]
            if (length(dimstats) != length(dimmargin) ||
                any(dimstats !=  dimmargin)) {
                # TODO: Make this a warning with appropriate consequent behaviour
                stop("length(STATS) or dim(STATS) do not match dim(x)[MARGIN]")
            }
        }
    }
    perm <- c(MARGIN, seq_along(dims)[-MARGIN])

    if (length(STATS) == 1) {
        rle <- Rle(STATS, prod(dims))
        A <- RleArray(rle, dims[perm])
    } else {
        rle <- Rle(STATS, prod(dims[-MARGIN]))
        if (length(dims) > 2) {

        }
        if (MARGIN == 1) {
            # Works for MARGIN = 1 over matrix
            A <- aperm(RleArray(rle, rev(dims)), rev(perm))
        } else if (MARGIN == 2) {
            # Works for MARGIN = 2 over matrix
            A <- aperm(RleArray(rle, dims), perm)
        } else {
            stop("Crap")
        }
    }
    FUN(x, aperm(A, order(perm)), ...)
}

.sweep <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
    if (is.array(x)) {
        return(sweep(x, MARGIN, STATS, FUN, check.margin, ...))
    } else if (is(x, "DelayedArray")) {
        .sweep_DelayedArray(x, MARGIN, STATS, FUN, check.margin, ...)
    } else {
        stop("Expected 'x' to be an array or DelayedArray")
    }
}
