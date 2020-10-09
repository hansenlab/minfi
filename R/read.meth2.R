## ToDo
##  - handle non-HDF5 backends; ie. both DelayedArray and non-DelayedArray
##  - code sharing between DelayedArray and non-DelayedArray
##  - DONE dealing with the case of arrays with unequal dimension 
##  - SEMI-DONE set chunkdim somehow in arguments, perhaps compression level as well
##  - DONE include output dir in arguments
##  - return with right annotation
##  - POSTPONED TO WAIT ON PETE return a fully saved SE
##  - POSTPONED TO WAIT ON PETE progress bar
##  - DONE verbosity



read.metharray2 <- function(basenames,
                            extended = FALSE,
                            verbose = FALSE,
                            dir = "my_h5_se",
                            nArrays = 1,
                            BPPARAM = bpparam("SerialParam"),
                            nArrays.chunk = 1) {
    ## This function assumes that all IDATs are identical.
    ## This is not a safe assumption; will be relaxed later.
    ## At least I check for it.
    
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

    ## First we check whether these are all the same array type and has the
    ## same number of probes

    nProbes <- sapply(G.files, readIDAT, what = "nSNPsRead") # FIXME: potentially parallize this
    arrayTypes <- cbind(do.call(rbind, lapply(nProbes, minfi:::.guessArrayTypes)),
                        size = nProbes)
    sameLength <- (length(unique(arrayTypes[, "size"])) == 1)
    sameArray <- (length(unique(arrayTypes[, "array"])) == 1)
    if (!sameArray) {
        cat("[read.metharray2] Trying to parse IDAT files from different arrays.\n")
        cat("  Inferred Array sizes and types:\n")
        print(arrayTypes[, c("array", "size")])
        stop("[read.metharray2] Trying to parse different IDAT files, of ",
             "different size and type." )
    }
    if (sameArray && !sameLength && !force) {
        stop("[read.metharray] Trying to parse IDAT files with different ",
             "array size but seemingly all of the same type.\n  You can force ",
             "this by 'force=TRUE', see the man page ?read.metharray")
    }
    if (sameArray && !sameLength && force) {
        ## get the rownames as bpiterate
        iterator <- function(files, nArrays) {
            grid <- RegularArrayGrid(length(files), nArrays)
            gridIdx <- lapply(seq_along(grid), function(xx) {
                ir <- ranges(grid[[xx]])
                seq(from = start(ir), to = end(ir))
            })
            b <- 0L
            function() {
                if(b == length(grid)) return(NULL)
                b <<- b + 1L
                gridIdx[[b]]
            }
        }
        computer <- function(files) {
            Reduce("intersect", lapply(files, readIDAT, what = "IlluminaIDs"))
        }
        reducer <- function(x, y) {
            intersect(x,y)
        }
        rownames.out <- bpiterate(iterator(G.files, nArrays = nArrays),
                                  FUN = computer,
                                  REDUCE = reducer,
                                  BPPARAM = BPPARAM)
    }
    if (sameArray && sameLength) {
        first.file <- readIDAT(G.files[1])
        rownames.out <- rownames(first.file$Quants)
    }
    ans_nrow <- length(rownames.out)
    grid <- RegularArrayGrid(c(ans_nrow, ans_ncol),
                             c(ans_nrow, nArrays)) 
    readAndCheck <- function(b, files, grid, sink, sink_lock, rownames.out) {
        ## Before calling this function, we ensure that rownames.out
        ## is a subset of the rownames in the file.
        cat(".")
        viewport <- grid[[b]]
        col_range <- ranges(viewport)[2L]
        col_indices <- seq(from = start(col_range), to = end(col_range))
        Quants <- lapply(files[col_indices], function(xx) {
            readIDAT(xx)[["Quants"]]
        })
        Mean <- do.call(cbind, lapply(Quants, function(xx) xx[rownames.out, "Mean", drop=FALSE]))
        write_block(sink, viewport = viewport, block = Mean)
    }
    Gmean_sink <- HDF5RealizationSink(
        dim = c(ans_nrow, ans_ncol),
        filepath = file.path(dir, "assays.h5"),
        dimnames = NULL, # NOTE: Never allow dimnames.
        type = "integer",
        name = "Green",
        chunkdim = c(ans_nrow, nArrays.chunk))
    on.exit(close(Gmean_sink), add = TRUE)
    sink_lock <- ipcid()
    on.exit(ipcremove(sink_lock), add = TRUE)
    if(verbose)
        message(sprintf("[read.metharray2]: reading %i Green files\n", length(G.files)))
    Gmean <- bplapply(seq_along(grid),
                      FUN = readAndCheck,
                      files = G.files,
                      grid = grid, 
                      sink = Gmean_sink,
                      sink_lock = sink_lock,
                      rownames.out = rownames.out,
                      BPPARAM = BPPARAM)
    ipcremove(sink_lock)
    Gmean <- as(Gmean_sink, "DelayedArray")
    Rmean_sink <- HDF5RealizationSink(
        dim = c(ans_nrow, ans_ncol),
        filepath = file.path(dir, "assays.h5"),
        dimnames = NULL, # NOTE: Never allow dimnames.
        type = "integer",
        name = "Red",
        chunkdim = c(ans_nrow, nArrays.chunk))
    on.exit(close(Rmean_sink), add = TRUE)
    if(verbose)
        message(sprintf("[read.metharray2]: reading %i Red files\n", length(R.files)))
    Rmean <- bplapply(seq_along(grid),
                      FUN = readAndCheck,
                      files = R.files,
                      grid = grid, 
                      sink = Rmean_sink,
                      sink_lock = sink_lock,
                      rownames.out = rownames.out,
                      BPPARAM = BPPARAM)
    Rmean <- as(Rmean_sink, "DelayedArray")
    out <- RGChannelSet(Red = Rmean, Green = Gmean)
    rownames(out) <- rownames.out
    colnames(out) <- names(G.files)
    colData(out)$filename
    ## TODO: set annotation slot
    out
}

