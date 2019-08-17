# Internal functions -----------------------------------------------------------

read.GenomeStudio <- function(filename) {
    colnames <- strsplit(readLines(filename, n = 9)[9], "\t")[[1]]
    colClasses <- rep("NULL", length(colnames))
    names(colClasses) <- colnames
    colClasses["TargetID"] <- "character"
    colClasses[grep("AVG_Beta", colnames)] <- "numeric"
    colClasses[grep("Signal_A", colnames)] <- "numeric"
    colClasses[grep("Signal_B", colnames)] <- "numeric"
    # NOTE: File has trailing tab...
    colClasses <- c(colClasses, dummy = "NULL")
    mat <- read.table(
        file = filename,
        header = TRUE,
        sep = "\t",
        comment.char = "",
        quote = "",
        skip = 8,
        colClasses = colClasses)
    sampleNames <- sub(
        "\\.AVG_Beta", "", grep("AVG_Beta", colnames(mat), value = TRUE))
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

# Exported functions -----------------------------------------------------------

makeGenomicRatioSetFromMatrix <- function(mat,rownames = NULL, pData = NULL,
                                          array = "IlluminaHumanMethylation450k",
                                          annotation = .default.450k.annotation,
                                          mergeManifest = FALSE,
                                          what = c("Beta","M")) {

    what <- match.arg(what)

    if (!is.matrix(mat)) {
        stop(sprintf("'mat' must be a matrix. It is a %s.", class(mat)))
    }

    if (is.null(rownames)) {
        rownames <- rownames(mat)
    } else {
        if (length(rownames) != nrow(mat)) {
            stop("Number of rows of mat and length of rownames must match.")
        }
        rownames(mat) <- rownames
    }

    if (is.data.frame(pData)) pData <- as(pData,"DataFrame")
    if (is.null(colnames(mat))) colnames(mat) <- seq_len(ncol(mat))
    if (is.null(pData)) {
        pData <- DataFrame(X1 = seq_len(ncol(mat)), row.names = colnames(mat))
    }
    if (!is(pData, "DataFrame")) {
        stop(sprintf(
            "'pData' must be DataFrame or data.frame. It is a %s.",
            class(pData)))
    }

    # Create GRanges
    ann <- .getAnnotationString(c(array = array, annotation = annotation))
    if (!require(ann, character.only = TRUE)) {
        stop(sprintf("cannot load annotation package %s", ann))
    }
    object <- get(ann)

    gr <- getLocations(object, mergeManifest = mergeManifest,
                       orderByLocation = TRUE)

    locusNames <- names(gr)

    # NOTE:Tthis might return NAs but it's ok
    # TODO: Fix this. return only what is sent
    common <- intersect(locusNames, rownames(mat))
    if (length(common) == 0) {
        stop("No rowname matches. 'rownames' need to match ",
             "IlluminaHumanMethylation450k probe names.")
    }

    # NOTE: we give no warning if some of the rownames have no match.
    ind1 <- match(common,rownames(mat))
    ind2 <- match(common,locusNames)

    preprocessing <- c(
        rg.norm = "Matrix converted with makeGenomicRatioSetFromMatrix")

    if (what == "Beta") {
        out <- GenomicRatioSet(
            gr = gr[ind2, ],
            Beta = mat[ind1, ,drop = FALSE],
            M = NULL,
            CN = NULL,
            colData = pData,
            annotation = c(array = array, annotation = annotation),
            preprocessMethod = preprocessing)
    } else {
        out <- GenomicRatioSet(
            gr = gr[ind2, ],
            Beta = NULL,
            M = mat[ind1, ,drop = FALSE],
            CN = NULL,
            colData = pData,
            annotation = c(array = array, annotation = annotation),
            preprocessMethod = preprocessing)
    }
    out
}


getGenomicRatioSetFromGEO <- function(GSE = NULL, path = NULL,
                                      array = "IlluminaHumanMethylation450k",
                                      annotation = .default.450k.annotation,
                                      what = c("Beta","M"),
                                      mergeManifest = FALSE,
                                      i = 1) {

    what <- match.arg(what)

    if (is.null(GSE) && is.null(path)) {
        stop("Either GSE or path must be supplied.")
    }

    # Read the GEO Main Files Information
    if (!is.null(GSE)) {
        gset <- getGEO(GSE)
    } else {
        gset <- getGEO(file.path(path, list.files(path, pattern = ".soft")))
    }

    if (length(gset) == 0) stop("Empty list retrieved from GEO.")
    if (length(gset) > 1) {
        warning(
            "More than one ExpressionSet found:\n",
            names(gset),
            "\nUsing entry ",
            i)
        gset <- gset[[i]]
    } else gset <- gset[[1]]
    platform <- annotation(gset)

    if (platform != "GPL13534") {
        warning(
            sprintf("%s is not the platform ID associated with IlluminaHumanMethylation450k. Should be GPL13534.",
                    platform))
    }
    if (what ==" Beta" && (min(exprs(gset)[, 1], na.rm = TRUE) < 0 ||
                           max(exprs(gset)[, 1], na.rm = TRUE) > 1 )) {
        warning("Values outside [0,1] detected. 'what' argument should not ",
                "be Beta.")
    }

    # GEO data is read, now ready to create a minfi object
    ann <- .getAnnotationString(c(array = array, annotation = annotation))
    if (!require(ann, character.only = TRUE)) {
        stop(sprintf("cannot load annotation package %s", ann))
    }
    object <- get(ann)

    gr <- getLocations(
        object = object,
        mergeManifest = mergeManifest,
        orderByLocation = TRUE)

    locusNames <- names(gr)

    sampleNames(gset) <- gset$title

    # TODO: We could call makeGenomicRatioSetFromMatrix() but rewrite to avoid
    # a copy of exprs(gset)
    common <- intersect(locusNames, fData(gset)$Name)
    if (length(common) == 0) {
        stop("No rowname matches. 'rownames' need to match ",
             "IlluminaHumanMethylation450k probe names.")
    }

    # NOTE: We give no warning if some rownames have no match
    ind1 <- match(common,fData(gset)$Name)
    ind2 <- match(common,locusNames)

    preprocessing <- c(rg.norm = paste0("See GEO ", GSE, " for details"))

    if (what == "Beta") {
        out <- GenomicRatioSet(
            gr = gr[ind2, ],
            Beta = exprs(gset)[ind1, ,drop = FALSE],
            M = NULL,
            CN = NULL,
            colData = as(pData(gset), "DataFrame"),
            annotation = c(array = array, annotation = annotation),
            preprocessMethod = preprocessing)
    } else {
        out <- GenomicRatioSet(
            gr = gr[ind2, ],
            Beta = NULL,
            M = exprs(gset)[ind1, ,drop = FALSE],
            CN = NULL,
            colData = as(pData(gset), "DataFrame"),
            annotation = c(array = array, annotation = annotation),
            preprocessMethod = preprocessing)
    }

    out
}

readTCGA <- function(filename, sep = "\t", keyName = "Composite Element REF",
                     Betaname = "Beta_value", pData = NULL,
                     array = "IlluminaHumanMethylation450k",
                     annotation = .default.450k.annotation,
                     mergeManifest = FALSE,
                     showProgress = TRUE) {
    # NOTE: We assume first column are sample names and second column are the
    #       value identifiers
    colnames <- strsplit(readLines(filename, n = 2), sep)

    select <- sort(
        c(grep(keyName, colnames[[2]]), grep(Betaname,colnames[[2]])))

    mat <- fread(
        input = filename,
        header = FALSE,
        sep = sep,
        select = select,
        showProgress = showProgress,
        skip = 2)
    rowNames <- as.matrix(mat[, 1, with = FALSE])
    mat <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat) <- rowNames
    colnames(mat) <- colnames[[1]][select][-1]
    # TODO: Shouldn't be necessary to rm() anything
    rm(rowNames,colnames)

    makeGenomicRatioSetFromMatrix(
        mat = mat,
        pData = pData,
        array = array,
        annotation = annotation,
        mergeManifest = mergeManifest, what = "Beta")
}


readGEORawFile <- function(filename, sep = ",", Uname = "Unmethylated signal",
                           Mname = "Methylated signal", row.names = 1,
                           pData = NULL, array = "IlluminaHumanMethylation450k",
                           annotation = .default.450k.annotation,
                           mergeManifest = FALSE, showProgress = TRUE, ...) {
    colnames <- strsplit(readLines(filename, n = 1), sep)[[1]]

    if (all(!grepl(Uname, colnames))) {
        stop("No columns contain Uname. Use readLines or look at file header ",
             "to see column names.")
    }

    if (all(!grepl(Mname, colnames))) {
        stop("No columns contain Mname. Use readLines or look at file header ",
             "to see column names.")
    }

    select <- sort(c(row.names, grep(Uname,colnames), grep(Mname,colnames)))

    mat <- fread(
        input = filename,
        header = TRUE,
        sep = sep,
        select = select,
        showProgress = showProgress,
        ...)

    rowNames <- as.matrix(mat[, 1, with = FALSE])
    mat <- as.matrix(mat[,-1, with = FALSE])
    rownames(mat) <- rowNames
    rm(rowNames)

    uindex <- grep(Uname,colnames(mat))
    mindex <- grep(Mname,colnames(mat))

    trim <- function(x){
        x <- gsub("^\\s+|\\s+$", "", x)
        x <- gsub("^\\.+|\\.+$", "", x)
        x <- gsub("^\\_+|\\_$", "", x)
        return(x)
    }

    UsampleNames <- trim(sub(Uname, "", colnames(mat)[uindex]))
    MsampleNames <- trim(sub(Mname, "", colnames(mat)[mindex]))

    index <- match(UsampleNames,MsampleNames)
    MsampleNames <- MsampleNames[index]
    mindex <- mindex[index]

    if (!identical(UsampleNames,MsampleNames)) {
        stop("Sample names do not match for Meth and Unmeth channels.")
    }

    if (is.data.frame(pData)) pData <- as(pData,"DataFrame")

    if (is.null(pData)) {
        pData <- DataFrame(
            X1 = seq_along(UsampleNames),
            row.names = UsampleNames)
    }

    ann <- .getAnnotationString(c(array = array, annotation = annotation))
    if (!require(ann, character.only = TRUE)) {
        stop(sprintf("cannot load annotation package %s", ann))
    }
    object <- get(ann)

    gr <- getLocations(
        object = object,
        mergeManifest = mergeManifest,
        orderByLocation = TRUE)

    locusNames <- names(gr)

    # NOTE: This might return NAs but it's ok
    # TODO: Fix this. return only what is sent
    common <- intersect(locusNames,rownames(mat))
    if (length(common) == 0) {
        stop("No rowname matches. 'rownames' need to match ",
             "IlluminaHumanMethylation450k probe names.")
    }

    # NOTE: We give no warning if some of the rownames have no match.
    ind1 <- match(common,rownames(mat))
    ind2 <- match(common,locusNames)

    preprocessing <- c(rg.norm = paste0("Data read from file ", filename, "."))
    colnames(mat) <- NULL
    GenomicMethylSet(
        gr =  gr[ind2, ],
        Meth = mat[ind1, mindex],
        Unmeth = mat[ind1, uindex],
        colData = pData,
        preprocessMethod = preprocessing,
        annotation = c(array = array, annotation = annotation))
}

# TODOs ------------------------------------------------------------------------

# TODO: Lots of duplicated code; DRY
