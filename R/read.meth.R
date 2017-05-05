.guessArrayTypes <- function(nProbes) {
    if(nProbes >= 622000 && nProbes <= 623000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation)
    } else if(nProbes >= 1050000 && nProbes <= 1053000) {
        ## "Current EPIC scan type"
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = .default.epic.annotation)
    } else if(nProbes >= 1032000 && nProbes <= 1033000) {
        ## "Old EPIC scan type"
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = .default.epic.annotation) 
    } else if(nProbes >= 54000 && nProbes <= 56000) {
        arrayAnnotation <- c(array = "IlluminaHumanMethylation27k", annotation = .default.27k.annotation)
    } else {
        arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
    }
    arrayAnnotation
}

read.metharray <- function(basenames, extended = FALSE, verbose = FALSE, force = FALSE) {
    basenames <- sub("_Grn\\.idat.*", "", basenames)
    basenames <- sub("_Red\\.idat.*", "", basenames)
    stopifnot(!anyDuplicated(basenames))
    
    G.files <- paste(basenames, "_Grn.idat", sep = "")
    names(G.files) <- basename(basenames)
    these.dont.exists <- !file.exists(G.files)
    if(any(these.dont.exists))
        G.files[these.dont.exists] <- paste0(G.files[these.dont.exists], ".gz")
    if(!all(file.exists(G.files))) {
        noexits <- sub("\\.gz", "", G.files[!file.exists(G.files)])
        stop("The following specified files do not exist:", paste(noexits, collapse = ", "))
    }
    R.files <- paste(basenames, "_Red.idat", sep = "")
    names(R.files) <- basename(basenames)
    these.dont.exists <- !file.exists(R.files)
    if(any(these.dont.exists))
        R.files[these.dont.exists] <- paste0(R.files[these.dont.exists], ".gz")
    if(!all(file.exists(R.files))) {
        noexits <- R.files[!file.exists(G.files)]
        stop("The following specified files do not exist:", paste(noexits, collapse = ", "))
    }
    stime <- system.time({
        G.idats <- lapply(G.files, function(xx) {
            if(verbose) message("[read.metharray] Reading", basename(xx))
            readIDAT(xx)
        })
        R.idats <- lapply(R.files, function(xx) {
            if(verbose) message("[read.metharray] Reading", basename(xx))
            readIDAT(xx)
        })
    })[3]
    if(verbose) message(sprintf("[read.metharray] Read idat files in %.1f seconds", stime))
    if(verbose) message("[read.metharray] Creating data matrices ... ", appendLF = FALSE)
    ptime1 <- proc.time()
    allNProbes <- sapply(G.idats, function(xx) nrow(xx$Quants))
    arrayTypes <- cbind(do.call(rbind, lapply(allNProbes, .guessArrayTypes)),
                        size = allNProbes)
    sameLength <- (length(unique(arrayTypes[, "size"])) == 1)
    sameArray <- (length(unique(arrayTypes[, "array"])) == 1)
    
    if(!sameLength && !sameArray) {
        cat("[read.metharray] Trying to parse IDAT files from different arrays.\n")
        cat("  Inferred Array sizes and types:\n")
        print(arrayTypes[, c("array", "size")])
        stop("[read.metharray] Trying to parse different IDAT files, of different size and type.")
    }
    if(!sameLength && sameArray && !force) {
        stop("[read.metharray] Trying to parse IDAT files with different array size but seemingly all of the same type.\n  You can force this by 'force=TRUE', see the man page ?read.metharray")
    }
    commonAddresses <- as.character(Reduce("intersect", lapply(G.idats, function(xx) rownames(xx$Quants))))
    GreenMean <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
    RedMean <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[commonAddresses, "Mean"]))
    if(extended) {
        GreenSD <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[commonAddresses, "SD"]))
        RedSD <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[commonAddresses, "SD"]))
        NBeads <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[commonAddresses, "NBeads"]))
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) message(sprintf("done in %.1f seconds", stime))
    if(verbose) message("[read.metharray] Instantiating final object ... ", appendLF = FALSE)
    ptime1 <- proc.time()
    if(extended) {
        out <- RGChannelSetExtended(Red = RedMean, Green = GreenMean,
                                    RedSD = RedSD, GreenSD = GreenSD, NBeads = NBeads)
    } else {
        out <- RGChannelSet(Red = RedMean, Green = GreenMean)
    }
    rownames(out) <- commonAddresses
    out@annotation <- c(array = arrayTypes[1,1], annotation = arrayTypes[1,2])
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) message(sprintf("done in %.1f seconds", stime))
    out
}

read.metharray.sheet <- function(base, pattern = "csv$", ignore.case = TRUE,
                                 recursive = TRUE, verbose = TRUE) {
    readSheet <- function(file) {
        dataheader <- grep("^\\[DATA\\]", readLines(file), ignore.case = TRUE)
        if(length(dataheader) == 0)
            dataheader <- 0
        df <- read.csv(file, stringsAsFactor = FALSE, skip = dataheader)
        if(length(nam <- grep("Sentrix_Position", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Array <- as.character(df[, nam])
            df[, nam] <- NULL
        }
        if(length(nam <- grep("Array[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Array <- as.character(df[, nam])
            df[, nam] <- NULL
        }
        if(! "Array" %in% names(df))
            warning(sprintf("Could not infer array name for file: %s", file))
        if(length(nam <- grep("Sentrix_ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Slide <- as.character(df[, nam])
            df[, nam] <- NULL
        }
        if(length(nam <- grep("Slide[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Slide <- as.character(df[, nam])
            df[, nam] <- NULL
        }
        if(! "Slide" %in% names(df))
            warning(sprintf("Could not infer slide name for file: %s", file))
        else
            df[, "Slide"] <- as.character(df[, "Slide"])
        if(length(nam <- grep("Plate[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Plate <- as.character(df[, nam])
            df[, nam] <- NULL
        }
        for(nam in c("Pool_ID", "Sample_Plate", "Sample_Well")) {
            if(nam %in% names(df)) {
                df[[nam]] <- as.character(df[[nam]])
            }
        }
            
        if(!is.null(df$Array)) {
            patterns <- sprintf("%s_%s_Grn.idat", df$Slide, df$Array)
            allfiles <- list.files(dirname(file), recursive = recursive, full.names = TRUE)
            basenames <- sapply(patterns, function(xx) grep(xx, allfiles, value = TRUE))
            names(basenames) <- NULL
            basenames <- sub("_Grn\\.idat.*", "", basenames, ignore.case = TRUE)
            df$Basename <- basenames
        }
        df
    }
    if(!all(file.exists(base)))
        stop("'base' does not exists")
    info <- file.info(base)
    if(!all(info$isdir) && !all(!info$isdir))
        stop("'base needs to be either directories or files")
    if(all(info$isdir)) {
        csvfiles <- list.files(base, recursive = recursive, pattern = pattern,
                               ignore.case = ignore.case, full.names = TRUE)
        if(verbose) {
            message("[read.metharray.sheet] Found the following CSV files:")
            print(csvfiles)
        }
    } else
        csvfiles <- list.files(base, full.names = TRUE)
    dfs <- lapply(csvfiles, readSheet)
    namesUnion <- Reduce(union, lapply(dfs, names))
    df <- do.call(rbind, lapply(dfs, function(df) {
        newnames <- setdiff(namesUnion, names(df))
        newdf <- matrix(NA, ncol = length(newnames), nrow = nrow(df), dimnames = list(NULL, newnames))
        cbind(df, as.data.frame(newdf))
    }))
    df
}
    

read.metharray.exp <- function(base = NULL, targets = NULL, extended = FALSE, 
                               recursive = FALSE, verbose = FALSE, force = FALSE) {
    if(!is.null(targets)) {
        if(! "Basename" %in% names(targets))
            stop("Need 'Basename' amongst the column names of 'targets'")
        if(!is.null(base)) {
            files <- file.path(base, basename(targets$Basename))
        } else {
            files <- targets$Basename
        }
        rgSet <- read.metharray(files, extended = extended, verbose = verbose, force = force)
        pD <- targets
        pD$filenames <- files
        rownames(pD) <- colnames(rgSet)
        colData(rgSet) <- as(pD, "DataFrame")
        return(rgSet)
    }
    ## Now we just read all files in the directory
    Grn.files <- list.files(base, pattern = "_Grn.idat$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    Red.files <- list.files(base, pattern = "_Red.idat$", recursive = recursive,
                            ignore.case = TRUE, full.names = TRUE)
    if(length(Grn.files) == 0 || length(Red.files) == 0)
        stop("No IDAT files were found")
    commonFiles <- intersect(sub("_Grn.idat$", "", Grn.files), sub("_Red.idat$", "", Red.files))
    if(length(commonFiles) == 0)
        stop("No IDAT files with both Red and Green channel were found")
    commonFiles.Grn <- paste(commonFiles, "_Grn.idat", sep = "")
    commonFiles.Red <- paste(commonFiles, "_Red.idat", sep = "")
    if(!setequal(commonFiles.Grn, Grn.files))
        warning(sprintf("the following files only exists for the green channel: %s",
                        paste(setdiff(Grn.files, commonFiles.Grn), collapse = ", ")))
    if(!setequal(commonFiles.Red, Red.files))
        warning(sprintf("the following files only exists for the red channel: %s",
                        paste(setdiff(Red.files, commonFiles.Red), collapse = ", ")))
    rgSet <- read.metharray(basenames = commonFiles, extended = extended, verbose = verbose, force = force)
    rgSet
}
