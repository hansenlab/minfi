read.manifest <- function(file, returnAll = FALSE) {
    if(!require(Biostrings))
        stop("read.manifest requires the package Biostrings")
    ## As is, requires grep
    control.line <- system(sprintf("grep -n Controls %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    ## Column headers is in line 8, hardcoded
    colNames <- readLines(file, n = 8)[8]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) = colNames
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(file, header = TRUE,
                           sep = ",", comment.char = "", quote = "", 
                           skip = 7, colClasses = colClasses, nrows = control.line - 9)
    TypeI <- manifest[manifest$Infinium_Design_Type == "I",
                        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel",
                          "Next_Base", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2,3,4,5,6,7)] <- c("AddressA", "AddressB", "Color", "NextBase",
                                        "ProbeSeqA", "ProbeSeqB")
    TypeI$nCpG <- as.integer(oligonucleotideFrequency(DNAStringSet(TypeI$ProbeSeqB),
                                                      width = 2)[, "CG"] - 1)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeII <- manifest[manifest$Infinium_Design_Type == "II",
                        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("Address", "ProbeSeq")
    TypeII$nCpG <- as.integer(letterFrequency(BStringSet(TypeII$ProbeSeq), letters = "R"))
    controls <- read.table(file, skip = control.line,
                           sep = ",", comment.char = "", quote = "",
                           colClasses = c(rep("character", 5)))
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    if(returnAll)
        list(manifestList = list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl),
             manifest = manifest, controls = controls)
    else 
        list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl)
}


read.manifest.27k <- function(file) {
    if(!require(Biostrings))
        stop("read.manifest requires the package Biostrings")
    ## As is, requires grep
    control.line <- system(sprintf("grep -a -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    assay.line <- system(sprintf("grep -a -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))
    ## Column headers is in line 8, hardcoded
    colNames <- tail(readLines(file, n = assay.line + 1), n= 1)
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) = colNames
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(file, header = TRUE,
                           sep = ",", comment.char = "", quote = "", 
                           skip = assay.line, colClasses = colClasses,
                           nrows = control.line - (assay.line + 1), fill = TRUE)
    TypeI <- manifest[ c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel",
                         "Next_Base", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2,3,4,5,6,7)] <- c("AddressA", "AddressB", "Color", "NextBase",
                                      "ProbeSeqA", "ProbeSeqB")
    TypeI$nCpG <- as.integer(oligonucleotideFrequency(DNAStringSet(TypeI$ProbeSeqB),
                                                      width = 2)[, "CG"] - 1)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeII <- manifest[manifest$Infinium_Design_Type == "II",
                        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("Address", "ProbeSeq")
    TypeII$nCpG <- as.integer(letterFrequency(BStringSet(TypeII$ProbeSeq), letters = "R"))
    controls <- read.table(file, skip = control.line,
                           sep = ",", comment.char = "", quote = "",
                           colClasses = c(rep("character", 5)))
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    list(manifestList = list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl),
             manifest = manifest, controls = controls)
}
