# TODO: type() for all RealizationSink subclasses
setMethod("type", "HDF5RealizationSink", function(x) {
    x@type
})
setMethod("type", "arrayRealizationSink", function(x) {
    x@result_envir$result
})
# TODO: dimnames() for all RealizationSink subclasses
setMethod("dimnames", "arrayRealizationSink", function(x) {
    dimnames(x@result_envir$result)
})

# NOTE: Conceptually, this is `sink[i, ] <- value`
# TODO: Formalise as `[<-`,RealizationSink-method
# TODO: Generalize to arbitrary i, j, `...`. This is a low priority and will
#       also be tricky.
subassignRowsToRealizationSink <- function(sink, i, value) {
    # Set up ArrayGrid over 'sink' and 'value'.
    # We're going to walk over columns of 'value', which imposes two
    # requirements:
    # 1. Need to increase the block length so each block is made of at least
    #    one column.
    # 2. The grid over 'sink' must have the same dim as the grid over 'value'
    sink_type <- type(sink)
    sink_dim <- dim(sink)
    max_block_len <- max(DelayedArray:::get_max_block_length(sink_type),
                         sink_dim[[1L]])

    value_grid <- defaultGrid(value, max_block_len)
    sink_grid <- RegularArrayGrid(refdim = sink_dim,
                                  spacings = c(sink_dim[[1L]],
                                               sink_dim[[2L]] /
                                                   length(value_grid)))
    stopifnot(dim(sink_grid) == dim(value_grid))
    nblock <- length(sink_grid)

    # TODO: It feels like 'sink_temp_array' shouldn't be necessary
    # Set up temporary 'Array' concrete subclass with appropriate dimensions
    # and type
    sink_dimnames <- dimnames(sink)
    # TODO: Check with minfi whether this should be filled with 0 or NA
    sink_temp_array <- RleArray(Rle(vector(sink_type, 1L), prod(sink_dim)),
                                dim = sink_dim,
                                dimnames = sink_dimnames)

    # Loop over blocks and write to sink
    bplapply(seq_len(nblock), function(b, i, sink_grid, value_grid, value,
                                       sink_temp_array) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }

        # Prepare 'sink_block'
        sink_viewport <- sink_grid[[b]]
        sink_block <- DelayedArray:::extract_block(sink_temp_array,
                                                   sink_viewport)
        if (!is.array(sink_block)) {
            sink_block <- DelayedArray:::.as_array_or_matrix(sink_block)
        }
        attr(sink_block, "from_grid") <- sink_grid
        attr(sink_block, "block_id") <- b

        # Prepare 'value_block'
        value_viewport <- value_grid[[b]]
        value_block <- DelayedArray:::extract_block(value, value_viewport)
        if (!is.array(value_block)) {
            value_block <- DelayedArray:::.as_array_or_matrix(value_block)
        }
        attr(value_block, "from_grid") <- value_grid
        attr(value_block, "block_id") <- b

        # Update 'sink_block' with 'value_block'
        sink_block[i, ] <- value_block

        # Write updated 'sink_block' to 'sink'
        write_block_to_sink(sink_block, sink, sink_viewport)
        sink_block <- NULL
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
        sink_block
    }, i = i, sink_grid = sink_grid, value_grid = value_grid, value = value,
    sink_temp_array = sink_temp_array)
}
