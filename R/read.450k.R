read.450k <- function(basenames, extended = FALSE, verbose = FALSE) {
    basenames <- sub("_Grn\\.idat$", "", basenames)
    basenames <- sub("_Red\\.idat$", "", basenames)
    G.files <- paste(basenames, "_Grn.idat", sep = "")
    names(G.files) <- basename(basenames)
    if(!all(file.exists(G.files))) {
        noexits <- G.files[!file.exists(G.files)]
        stop("The following specified files do not exist:", paste(noexits, collapse = ", "))
    }
    R.files <- paste(basenames, "_Red.idat", sep = "")
    names(R.files) <- basename(basenames)
    if(!all(file.exists(R.files))) {
        noexits <- R.files[!file.exists(G.files)]
        stop("The following specified files do not exist:", paste(noexits, collapse = ", "))
    }
    stime <- system.time({
        G.idats <- lapply(G.files, function(xx) {
            if(verbose) cat("[read.450k] Reading", basename(xx), "\n")
            readIDAT(xx)
        })
        R.idats <- lapply(R.files, function(xx) {
            if(verbose) cat("[read.450k] Reading", basename(xx), "\n")
            readIDAT(xx)
        })
    })[3]
    if(verbose) cat("[read.450k] Read idat files in ", stime, "seconds\n")
    if(verbose) cat("[read.450k] Creating data matrices ... ")
    ptime1 <- proc.time()
    GreenMean <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "Mean"]))
    RedMean <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[, "Mean"]))
    if(extended) {
        GreenSD <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "SD"]))
        RedSD <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[, "SD"]))
        NBeads <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "NBeads"]))
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat("done in", stime, "seconds\n")
    if(verbose) cat("[read.450k] Instantiating final object ... ")
    ptime1 <- proc.time()
    if(extended) {
        out <- new("RGChannelSetExtended", Red = RedMean, Green = GreenMean,
                   RedSD = RedSD, GreenSD = GreenSD, NBeads = NBeads)
    } else {
        out <- new("RGChannelSet", Red = RedMean, Green = GreenMean)
    }
    featureNames(out) <- rownames(G.idats[[1]]$Quants)
    annotation(out) <- c(array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat("done in", stime, "seconds\n")
    out
}

read.450k.sheet <- function(base, pattern = "csv$", ignore.case = TRUE,
                            recursive = TRUE, verbose = TRUE) {
    readSheet <- function(file) {
        dataheader <- grep("^\\[DATA\\]", readLines(file), ignore.case = TRUE)
        if(length(dataheader) == 0)
            dataheader <- 0
        df <- read.csv(file, stringsAsFactor = FALSE, skip = dataheader)
        if(length(nam <- grep("Sentrix_Position", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Array <- df[, nam]
            df[, nam] <- NULL
        }
        if(length(nam <- grep("Array[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Array <- df[, nam]
            df[, nam] <- NULL
        }
        if(! "Array" %in% names(df))
            warning(sprintf("Could not infer array name for file: %s", file))
        if(length(nam <- grep("Sentrix_ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Slide <- df[, nam]
            df[, nam] <- NULL
        }
        if(length(nam <- grep("Slide[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Slide <- df[, nam]
            df[, nam] <- NULL
        }
        if(! "Slide" %in% names(df))
            warning(sprintf("Could not infer slide name for file: %s", file))
        if(length(nam <- grep("Plate[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
            df$Plate <- df[, nam]
            df[, nam] <- NULL
        }
        if(!is.null(df$Array)) {
            patterns <- sprintf("%s_%s_Grn.idat", df$Slide, df$Array)
            allfiles <- list.files(dirname(file), recursive = recursive, full.names = TRUE)
            basenames <- sapply(patterns, function(xx) grep(xx, allfiles, value = TRUE))
            names(basenames) <- NULL
            basenames <- sub("_Grn\\.idat", "", basenames, ignore.case = TRUE)
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
            cat("[read.450k.sheet] Found the following CSV files:\n")
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
    

read.450k.exp <- function(base, targets = NULL, extended = FALSE, 
                          recursive = FALSE, verbose = FALSE) {
    if(!is.null(targets)) {
        if(! "Basename" %in% names(targets))
            stop("Need 'Basename' amongst the column names of 'targets'")
        files <- targets$Basename
        rgSet <- read.450k(files, extended = extended, verbose = verbose)
        pD <- targets
        pD$filenames <- files
        rownames(pD) <- sampleNames(rgSet)
        pData(rgSet) <- pD
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
    rgSet <- read.450k(basenames = commonFiles, extended = extended, verbose = verbose)
    rgSet
}

read.GenomeStudio <- function(filename) {
    colnames <- strsplit(readLines(filename, n = 9)[9], "\t")[[1]]
    colClasses <- rep("NULL", length(colnames))
    names(colClasses) <- colnames
    colClasses["TargetID"] <- "character"
    colClasses[grep("AVG_Beta", colnames)] <- "numeric"
    colClasses[grep("Signal_A", colnames)] <- "numeric"
    colClasses[grep("Signal_B", colnames)] <- "numeric"
    colClasses <- c(colClasses, dummy = "NULL") # File has trailing tab...
    mat <- read.table(filename, header = TRUE, sep = "\t",
                      comment.char = "", quote = "", 
                      skip = 8, colClasses = colClasses)
    sampleNames <- sub("\\.AVG_Beta", "", grep("AVG_Beta", colnames(mat), value = TRUE))
    beta <- mat[, grep("AVG_Beta", colnames(mat))]
    colnames(beta) <- sampleNames
    rownames(beta) <- mat$TargetID
    SignalA <- mat[, grep("Signal_A", colnames(mat))]
    colnames(SignalA) <- sampleNames
    rownames(SignalA) <- mat$TargetID
    SignalB <- mat[, grep("Signal_B", colnames(mat))]
    colnames(SignalB) <- sampleNames
    rownames(SignalB) <- mat$TargetID
    list(beta = beta, SignalA = SignalA, SignalB = SignalB)
}
