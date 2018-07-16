## ToDo
##  - handle non-HDF5 backends; ie. both DelayedArray and non-DelayedArray
##  - code sharing between DelayedArray and non-DelayedArray
##  - dealing with the case of arrays with unequal dimension
##  - set chunkdim somehow in arguments, perhaps compression level as well
##  - include output dir in arguments
##  - return with right annotation
##  - return a fully saved SE
##  - progress bar; verbosity

read.metharray2 <- function(basenames, extended = FALSE, verbose = FALSE,
                            nArrays = 1, BPPARAM = bpparam("SerialParam") ) {
    ## This function assumes that all IDATs are identical.
    ## This is not a safe assumption; will be relaxed later.

    ## Setting up file names
    basenames <- sub("_Grn\\.idat.*", "", basenames)
    basenames <- sub("_Red\\.idat.*", "", basenames)
    stopifnot(!anyDuplicated(basenames))
    G.files <- paste(basenames, "_Grn.idat", sep = "")
    names(G.files) <- basename(basenames)
    these.dont.exists <- !file.exists(G.files)
    if (any(these.dont.exists)) {
        G.files[these.dont.exists] <- paste0(G.files[these.dont.exists], ".gz")
    }
    if (!all(file.exists(G.files))) {
        noexits <- sub("\\.gz", "", G.files[!file.exists(G.files)])
        stop("The following specified files do not exist:",
             paste(noexits, collapse = ", "))
    }
    R.files <- paste(basenames, "_Red.idat", sep = "")
    names(R.files) <- basename(basenames)
    these.dont.exists <- !file.exists(R.files)
    if (any(these.dont.exists)) {
        R.files[these.dont.exists] <- paste0(R.files[these.dont.exists], ".gz")
    }
    if (!all(file.exists(R.files))) {
        noexits <- R.files[!file.exists(G.files)]
        stop("The following specified files do not exist:",
             paste(noexits, collapse = ", "))
    }
    ans_ncol <- length(G.files)
    first.file <- readIDAT(G.files[1])
    rownames.check <- rownames(first.file$Quants)
    ans_nrow <- length(rownames.check)
    grid <- RegularArrayGrid(c(ans_nrow, ans_ncol),
                             c(ans_nrow, nArrays)) 
    readAndCheck <- function(b, files, grid, sink, sink_lock, rownames.check) {
        cat(".")
        viewport <- grid[[b]]
        col_range <- ranges(viewport)[2L]
        col_indices <- seq(from = start(col_range), to = end(col_range))
        Quants <- lapply(files[col_indices], function(xx) {
            readIDAT(xx)[["Quants"]]
        })
        sapply(Quants, function(xx) stopifnot(identical(rownames(xx), rownames.check)))
        Mean <- do.call(cbind, lapply(Quants, function(xx) xx[,"Mean", drop=FALSE]))
        write_block(x = sink, viewport = viewport, block = Mean)
    }
    Gmean_sink <- HDF5RealizationSink(
        dim = c(ans_nrow, ans_ncol),
        ## NOTE: Never allow dimnames.
        dimnames = NULL,
        type = "integer",
        name = "Green",
        chunkdim = c(ans_nrow, 1))
    on.exit(close(Gmean_sink), add = TRUE)
    sink_lock <- ipcid()
    on.exit(ipcremove(sink_lock), add = TRUE)
    Gmean <- bplapply(seq_along(grid),
                  FUN = readAndCheck,
                  files = G.files,
                  grid = grid, 
                  sink = Gmean_sink,
                  sink_lock = sink_lock,
                  rownames.check = rownames.check,
                  BPPARAM = BPPARAM)
    ipcremove(sink_lock)
    Gmean <- as(Gmean_sink, "DelayedArray")
    Rmean_sink <- HDF5RealizationSink(
        dim = c(ans_nrow, ans_ncol),
        ## NOTE: Never allow dimnames.
        dimnames = NULL,
        type = "integer",
        name = "Red",
        chunkdim = c(ans_nrow, 1))
    on.exit(close(Rmean_sink), add = TRUE)
    Rmean <- bplapply(seq_along(grid),
                  FUN = readAndCheck,
                  files = R.files,
                  grid = grid, 
                  sink = Rmean_sink,
                  sink_lock = sink_lock,
                  rownames.check = rownames.check,
                  BPPARAM = BPPARAM)
    Rmean <- as(Rmean_sink, "DelayedArray")
    out <- RGChannelSet(Red = Rmean, Green = Gmean)
    rownames(out) <- rownames.check
    colnames(out) <- names(G.files)
    colData(out)$filename
    ## TODO: set annotation slot
    out
}

