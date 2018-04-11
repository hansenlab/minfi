# ------------------------------------------------------------------------------
# Missing methods
#

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

# ------------------------------------------------------------------------------
# Helper functions for setting up ArrayGrid instances
#

colGrid <- function(x) {
    max_block_len <- max(nrow(x), DelayedArray:::get_max_block_length(type(x)))
    spacings <- DelayedArray:::get_spacings_for_linear_capped_length_blocks(
        refdim = dim(x),
        max_block_len = max_block_len)
    RegularArrayGrid(dim(x), spacings)
}

rowGrid <- function(x) {
    max_block_len <- max(ncol(x), DelayedArray:::get_max_block_length(type(x)))
    spacings <- DelayedArray:::get_spacings_for_hypercube_capped_length_blocks(
        refdim = dim(x),
        max_block_len = max_block_len)
    RegularArrayGrid(dim(x), spacings)
}

# ------------------------------------------------------------------------------
# Advanced block processing routines
#

# NOTE: DelayedArray::blockApply() with the option to write the blocks to
#       'sink'. Useful, for example, to apply a function across column-blocks
#       of a DelayedMatrix, write these results to disk, and then wrap
#       these in a DelayedMatrix.
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/10
blockApplyWithRealization <- function(x, FUN, ..., sink = NULL, x_grid = NULL,
                                      sink_grid = NULL, BPREDO = list(),
                                      BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)

    # Check conformable input_grids and output_grids
    x_grid <- DelayedArray:::.normarg_grid(x_grid, x)
    sink_grid <- DelayedArray:::.normarg_grid(sink_grid, sink)
    if (!identical(dim(x_grid), dim(sink_grid))) {
        stop("non-conformable 'x_grid' and 'sink_grid'")
    }

    # Loop over blocks of `x` and write to `sink`
    nblock <- length(x_grid)
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        x_viewport <- x_grid[[b]]
        sink_viewport <- sink_grid[[b]]
        block <- DelayedArray:::extract_block(x, x_viewport)
        if (!is.array(block)) {
            block <- DelayedArray:::.as_array_or_matrix(block)
        }
        attr(block, "from_grid") <- grid
        attr(block, "block_id") <- b
        block_ans <- FUN(block, ...)
        # NOTE: This is the only part different from DelayedArray::blockApply()
        if (!is.null(sink)) {
            write_block_to_sink(block_ans, sink, sink_viewport)
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
# NOTE: Different from DelayedArray:::block_Mapply(); designed to have an API
#       more like DelayedArray::blockArray()
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/11
blockMapply <- function(FUN, ..., MoreArgs = NULL, grids = NULL,
                        BPREDO = list(), BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)
    dots <- unname(list(...))
    dims <- lapply(dots, dim)
    all_same_dims <- all(
        vapply(dims, function(dim) all(dim == dims[[1L]]), logical(1L)))
    if (!all_same_dims) {
        stop("non-conformable arrays")
    }
    if (is.null(grids)) {
        grids <- replicate(length(dots), NULL)
    }
    grids <- mapply(
        FUN = DelayedArray:::.normarg_grid,
        grids,
        dots,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE)
    grid_dims <- lapply(grids, dim)
    all_same_grid_dims <- all(
        vapply(X = grid_dims,
               FUN = function(dim) all(dim == grid_dims[[1L]]),
               FUN.VALUE = logical(1L)))
    if (!all_same_grid_dims) {
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
        blocks <- mapply(
            FUN = function(x, grid, viewport) {
                block <- DelayedArray:::extract_block(x, viewport)
                if (!is.array(block)) {
                    block <- DelayedArray:::.as_array_or_matrix(block)
                }
                attr(block, "from_grid") <- grid
                attr(block, "block_id") <- b
                block
            },
            x = dots,
            grid = grids,
            viewport = viewports,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE)
        block_ans <- do.call(FUN, c(blocks, MoreArgs))
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
        block_ans
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}

# NOTE: blockMapply() with the option to write the blocks to multiple 'sinks'.
#       Useful, for example, to apply a function across column-blocks of
#       multiple DelayedMatrix objects, write these results to disk, and then
#       wrap these in a DelayedMatrix.
# NOTE: `input_grids`, `output_grids`, and `sinks` should all be lists
# TODO: See https://github.com/Bioconductor/DelayedArray/issues/11
blockMapplyWithRealization <- function(FUN, ..., MoreArgs = NULL, sinks = NULL,
                                       input_grids = NULL, output_grids = NULL,
                                       BPREDO = list(), BPPARAM = bpparam()) {
    FUN <- match.fun(FUN)
    # Check conformable inputs
    dots <- unname(list(...))
    dims <- lapply(dots, dim)
    all_same_dims <- all(
        vapply(dims, function(dim) all(dim == dims[[1L]]), logical(1L)))
    if (!all_same_dims) {
        stop("non-conformable arrays")
    }
    # Check valid `sinks`
    stopifnot(is.null(sinks) || is.list(sinks))
    # Check conformable input_grids and output_grids
    if (is.null(input_grids)) {
        input_grids <- replicate(length(dots), NULL)
    } else {
        stopifnot(is.list(input_grids))
    }
    input_grids <- mapply(
        FUN = DelayedArray:::.normarg_grid,
        input_grids,
        dots,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE)
    if (is.null(output_grids)) {
        output_grids <- replicate(length(sinks), NULL)
    } else {
        stopifnot(is.list(output_grids))
    }
    output_grids <- mapply(
        FUN = DelayedArray:::.normarg_grid,
        output_grids,
        sinks,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE)
    grid_dims <- lapply(c(input_grids, output_grids), dim)
    all_same_grid_dims <- all(
        vapply(X = grid_dims,
               FUN = function(dim) all(dim == grid_dims[[1L]]),
               FUN.VALUE = logical(1L)))
    if (!all_same_grid_dims) {
        stop("non-conformable 'input_grids' and 'output_grids'")
    }
    stopifnot(length(dots) == length(input_grids),
              length(sinks) == length(output_grids))

    # Loop over blocks of `dots` and write to `sinks`
    nblock <- length(input_grids[[1]])
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        input_viewports <- lapply(input_grids, function(grid) grid[[b]])
        output_viewports <- lapply(output_grids, function(grid) grid[[b]])
        blocks <- mapply(
            FUN = function(x, grid, viewport) {
                block <- DelayedArray:::extract_block(x, viewport)
                if (!is.array(block)) {
                    block <- DelayedArray:::.as_array_or_matrix(block)
                }
                attr(block, "from_grid") <- grid
                attr(block, "block_id") <- b
                block
            },
            x = dots,
            grid = input_grids,
            viewport = input_viewports,
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE)
        block_ans <- do.call(FUN, c(blocks, MoreArgs))
        if (!is.list(block_ans)) {
            block_ans <- list(block_ans)
        }
        # NOTE: This is the only part different from blockMapply()
        if (!is.null(sinks)) {
            mapply(function(ans, sink, viewport) {
                write_block_to_sink(ans, sink, viewport)
            }, ans = block_ans, sink = sinks, viewport = output_viewports,
            SIMPLIFY = FALSE, USE.NAMES = FALSE)
            block_ans <- NULL
        }
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
        block_ans
    },
    BPREDO = BPREDO,
    BPPARAM = BPPARAM)
}
