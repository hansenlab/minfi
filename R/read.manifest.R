read.manifest <- function(file) {
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
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(oligonucleotideFrequency(TypeI$ProbeSeqB,
                                                      width = 2)[, "CG"] - 1)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name),]
    TypeI <- TypeI[-grep("^rs", TypeI$Name),]
    
    TypeII <- manifest[manifest$Infinium_Design_Type == "II",
                        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name),]
    TypeII <- TypeII[-grep("^rs", TypeII$Name),]
    
    controls <- read.table(file, skip = control.line,
                           sep = ",", comment.char = "", quote = "",
                           colClasses = c(rep("character", 5)))
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")
    
    list(manifestList = list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl,
         TypeSnpI = TypeSnpI, TypeSnpII = TypeSnpII),
         manifest = manifest, controls = controls)
}


read.manifest.27k <- function(file) {
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
    TypeI <- TypeI[TypeI$Name != "",]
    names(TypeI)[c(2,3,4,5,6,7)] <- c("AddressA", "AddressB", "Color", "NextBase",
                                      "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(oligonucleotideFrequency(TypeI$ProbeSeqB,
                                                      width = 2)[, "CG"] - 1)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L

    TypeII <- manifest[manifest$Infinium_Design_Type == "II",
                        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- BStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeq, letters = "R"))
    
    controls <- read.table(file, skip = control.line,
                           sep = ",", comment.char = "", quote = "",
                           colClasses = c(rep("character", 5)))
    
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")
    
    snps <- TypeControl[TypeControl$Type == "Genotyping",]
    TypeControl <- TypeControl[TypeControl$Type != "Genotyping",]
    rsname <- sub("_[AB]", "", snps$ExtendedType)
    snps.sp <- split(snps, rsname)
    snps.sp <- lapply(names(snps.sp), function(rs) {
        snp <- snps.sp[[rs]]
        DataFrame(Name = rs,
                  AddressA = snp[grep("_A", snp$ExtendedType), "Address"],
                  AddressB = snp[grep("_B", snp$ExtendedType), "Address"],
                  Color = "Unknown")
    })
    TypeSnpI <- do.call(rbind, snps.sp)
    TypeSnpII <- TypeSnpI[0,]
                  
    list(manifestList = list(TypeI = TypeI, TypeII = TypeII, TypeControl = TypeControl,
         TypeSnpI = TypeSnpI, TypeSnpII = TypeSnpII),
         manifest = manifest, controls = controls)
}
