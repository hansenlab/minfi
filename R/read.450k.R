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
        else
            df[, "Slide"] <- as.character(df[, "Slide"])
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
    

read.450k.exp <- function(base = NULL, targets = NULL, extended = FALSE, 
                          recursive = FALSE, verbose = FALSE) {
    if(!is.null(targets)) {
        if(! "Basename" %in% names(targets))
            stop("Need 'Basename' amongst the column names of 'targets'")
        if(!is.null(base)) {
            files <- file.path(base, targets$Basename)
        } else {
            files <- targets$Basename
        }
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

makeGenomicRatioSetFromMatrix <- function(mat,rownames=NULL,
                                          pData=NULL,
                                          array = "IlluminaHumanMethylation450k",
                                          annotation=.default.450k.annotation,
                                          mergeManifest = FALSE,
                                          what=c("Beta","M")){

    what <- match.arg(what)

    if(!is.matrix(mat)) stop(sprintf("'mat' must be a matrix. It is a %s.",class(mat)))

    if(is.null(rownames)) rownmaes <- rownames(mat) else{
        if(length(rownames)!=nrow(mat))
            stop("Number of rows of mat and length of rownames must match.")
        rownames(mat) <- rownames
    }

    if(is.data.frame(pData)) pData <- as(pData,"DataFrame")

    if(is.null(colnames(mat))) colnames(mat) <- 1:ncol(mat)

    if(is.null(pData))  pData <- DataFrame( X1=1:ncol(mat), row.names=colnames(mat))

    if(class(pData)!="DataFrame")
        stop(sprintf("'pData' must be DataFrame or data.frame. It is a %s.",class(pData)))
    
    ##Create granges
    ann <- .getAnnotationString(c(array=array,annotation=annotation))
    if(!require(ann, character.only = TRUE))
        stop(sprintf("cannot load annotation package %s", ann))
    object <- get(ann)

    gr <- getLocations(object, mergeManifest = mergeManifest,
                                 orderByLocation = TRUE)

    locusNames <- names(gr)
     
    ##this might return NAs but it's ok
    ###fix this. return only what is sent
    common <- intersect(locusNames,rownames(mat))
    if(length(common)==0)
        stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")
    ##Note we give no warning if some of the rownmaes have no match.

    ind1 <- match(common,rownames(mat))
    ind2 <- match(common,locusNames)

    preprocessing <- c(rg.norm='Matrix converted with makeGenomicRatioSetFromMatrix')
    
    if(what=="Beta"){
        out <- GenomicRatioSet(gr=gr[ind2,],
                               Beta=mat[ind1,,drop=FALSE],
                               M =NULL,
                               CN=NULL,
                               pData=pData,
                               annotation=c(array=array,annotation=annotation),
                               preprocessMethod=preprocessing)
    } else {
        out <- GenomicRatioSet(gr=gr[ind2,],
                               Beta=NULL,
                               M=mat[ind1,,drop=FALSE],
                               CN=NULL,
                               pData=pData,
                               annotation=c(array=array,annotation=annotation),
                               preprocessMethod=preprocessing)
    }
    return(out)
}

    
getGenomicRatioSetFromGEO <- function(GSE=NULL,path=NULL,
                                      array = "IlluminaHumanMethylation450k",
                                      annotation=  .default.450k.annotation,
                                      what=c("Beta","M"),
                                      mergeManifest=FALSE,
                                      i=1) { 
    
    what <- match.arg(what)
    
    if(is.null(GSE) && is.null(path))
        stop("Either GSE or path must be supplied.")

    ##Reading the GEO Main Files Information
    if(!is.null(GSE)) gset <- GEOquery::getGEO(GSE) else  gset <- GEOquery::getGEO(filename = file.path(path, list.files(path, pattern = ".soft")))

    if(length(gset)==0) stop("Empty list retrieved from GEO.")
    if(length(gset)>1){
        warning("More than one ExpressionSet found:\n",names(gset),"\nUsing entry ",i)
        gset <- gset[[i]]
    } else gset <- gset[[1]]
    platform <- annotation(gset)

    if(platform!="GPL13534")
        warning(sprintf("%s is not the platform ID associated with IlluminaHumanMethylation450k. Should be GPL13534.",platform))
    if(what=="Beta" & (min(exprs(gset)[,1],na.rm=TRUE)<0 | max(exprs(gset)[,1],na.rm=TRUE)>1 ))
        warning("Values outside [0,1] detected. 'what' argument should not be Beta.")

    ##GEO data read, now ready to create a minfi object
    ann <- .getAnnotationString(c(array=array,annotation=annotation))
    if(!require(ann, character.only = TRUE))
        stop(sprintf("cannot load annotation package %s", ann))
    object <- get(ann)

    gr <- getLocations(object, mergeManifest = mergeManifest,
                                 orderByLocation = TRUE)
    
    locusNames <- names(gr)
     
    sampleNames(gset) <- gset$title
        
    ##we could call makeGenomicRatioSetFromMatrix but rewrite to
    ##avoide  a copy of exprs(gset)
    common <- intersect(locusNames,fData(gset)$Name)
    if(length(common)==0)
        stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")
    ##Note we give no warning if some rownames have no match
    
    ind1 <- match(common,fData(gset)$Name)
    ind2 <- match(common,locusNames)

    preprocessing <- c(rg.norm=paste0('See GEO ',GSE,' for details'))
    
    if(what=="Beta"){
        out <- GenomicRatioSet(gr=gr[ind2,],
                               Beta=exprs(gset)[ind1,,drop=FALSE],
                               M =NULL,
                               CN=NULL,
                               pData=pData(gset),
                               annotation=c(array=array,annotation=annotation),
                               preprocessMethod=preprocessing)
    } else {
        out <- GenomicRatioSet(gr=gr[ind2,],
                               Beta=NULL,
                               M=exprs(gset)[ind1,,drop=FALSE],
                               CN=NULL,
                               pData=pData(gset),
                               annotation=c(array=array,annotation=annotation),
                               preprocessMethod=preprocessing)

    }

    return(out)
}

readTCGA <- function(filename,sep="\t",
                      keyName="Composite Element REF",
                      Betaname="Beta_value",
                      pData=NULL,
                      array = "IlluminaHumanMethylation450k",
                      annotation=.default.450k.annotation,
                      mergeManifest = FALSE,
                      showProgress=TRUE){

    if (!requireNamespace("data.table", quietly = TRUE))
        stop("You need to install the data.table package from CRAN.")

    ##we assume first column are sample names
    ## and second column are the value identifiers
    colnames <- strsplit(readLines(filename, n = 2), sep)
    
    select <- sort(c(grep(keyName, colnames[[2]]), grep(Betaname,colnames[[2]])))
    
    mat <- data.table::fread(filename, header = FALSE, sep = sep, select=select,
                             showProgress=showProgress,skip=2)
    rowNames <- as.matrix(mat[, 1, with=FALSE])
    mat <- as.matrix(mat[, -1, with=FALSE])
    rownames(mat) <- rowNames
    colnames(mat)<-colnames[[1]][select][-1]
    rm(rowNames,colnames)

    return(makeGenomicRatioSetFromMatrix(mat, pData=pData, array=array,
                                         annotation=annotation,
                                         mergeManifest=mergeManifest,what="Beta"))
}


readGEORawFile <- function(filename,sep=",",
                            Uname="Unmethylated signal",
                            Mname="Methylated signal",
                            row.names=1,
                            pData=NULL,
                            array = "IlluminaHumanMethylation450k",
                            annotation=.default.450k.annotation,
                            mergeManifest = FALSE,
                            showProgress=TRUE){

    if (!requireNamespace("data.table", quietly = TRUE))
        stop("You need to install the data.table package from CRAN.")

    colnames <- strsplit(readLines(filename, n = 1), sep)[[1]]

    if(all(!grepl(Uname,colnames)))
        stop("No columns contain Uname. Use readLines or look at file header to see column names.")

    if(all(!grepl(Mname,colnames)))
        stop("No columns contain Mname. Use readLines or look at file header to see column names.")   

    select <- sort(c(row.names, grep(Uname,colnames),grep(Mname,colnames)))

    mat <- data.table::fread(filename, header = TRUE, sep = sep, select=select,
                             showProgress=showProgress)

    rowNames <- as.matrix(mat[,1,with=FALSE])
    mat <- as.matrix(mat[,-1,with=FALSE])
    rownames(mat) <- rowNames
    rm(rowNames)
    
    uindex <- grep(Uname,colnames(mat))
    mindex <- grep(Mname,colnames(mat))
  
    trim <- function (x){
        x<-gsub("^\\s+|\\s+$", "", x)
        x<-gsub("^\\.+|\\.+$", "", x)
        x<-gsub("^\\_+|\\_$", "", x)
        return(x)
    }
    
    UsampleNames <- trim(sub(Uname, "", colnames(mat)[uindex]))
    MsampleNames <- trim(sub(Mname, "", colnames(mat)[mindex]))

    index <- match(UsampleNames,MsampleNames)
    MsampleNames <- MsampleNames[index]
    mindex <- mindex[index]
    
    if(!identical(UsampleNames,MsampleNames))
        stop("Sample names do not match for Meth and Unmeth channels.")
    
    if(is.data.frame(pData)) pData <- as(pData,"DataFrame")

    if(is.null(pData))  pData <- DataFrame(
        X1=seq_along(UsampleNames),
        row.names=UsampleNames)

    ann <- .getAnnotationString(c(array=array,annotation=annotation))
    if(!require(ann, character.only = TRUE))
        stop(sprintf("cannot load annotation package %s", ann))
    object <- get(ann)

    gr <- getLocations(object, mergeManifest = mergeManifest,
                                 orderByLocation = TRUE)

    locusNames <- names(gr)
     
    ##this might return NAs but it's ok
    ###fix this. return only what is sent
    common <- intersect(locusNames,rownames(mat))
    if(length(common)==0)
        stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")
    ##Note we give no warning if some of the rownmaes have no match.

    ind1 <- match(common,rownames(mat))
    ind2 <- match(common,locusNames)

    preprocessing <- c(rg.norm=paste0("Data read from file ",filename,"."))
    
    return(GenomicMethylSet(gr =  gr[ind2,],
                            Meth = mat[ind1,mindex],
                            Unmeth = mat[ind1,uindex],
                            pData = pData,
                            preprocessMethod = preprocessing,
                            annotation = c(array=array,annotation=annotation)))
}




