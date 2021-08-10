# Internal functions -----------------------------------------------------------

read.manifest.Mammal <- function(file, type = 2) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 &&
                  is.integer(control.line) &&
                  !is.na(control.line))
    if(type == 1) {
        assay.line <- system(
            sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE)
        assay.line <- as.integer(sub(":.*", "", assay.line))
        stopifnot(length(assay.line) == 1 &&
                  is.integer(assay.line) &&
                  !is.na(assay.line))
    } else {
        assay.line  <- 0
    }
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"
    if(type == 1) {
        colClasses <- c(colClasses, drop = "logical")
    }
    manifest <- read.table(
        file = file,
        header = FALSE,
        col.names = names(colClasses),
        sep = ",",
        comment.char = "",
        quote = "",
        skip = assay.line + 1L,
        colClasses = colClasses,
        nrows = control.line - assay.line - 2L)
    manifest$drop <- NULL
    if(type == 2) {
        names(manifest)[c(1,2,30)] <- c("Name", "Internal_Name", "Internal_Name2")
    }
    
    manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)
    
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c(
        "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[-grep("^rs", TypeI$Name), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[-grep("^rs", TypeII$Name), ]

    if(type == 1) {
        controls <- read.table(
            file = file,
            skip = control.line,
            sep = ",",
            comment.char = "",
            quote = "",
            colClasses = c(rep("character", 5)))[, 1:5]
        TypeControl <- controls[, 1:4]
        names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
        TypeControl <- as(TypeControl, "DataFrame")
    } else {
        controls <- read.table(
            file = file,
            skip = control.line,
            sep = ",",
            comment.char = "",
            quote = "",
            colClasses = c(rep("character", 4)))[, 1:4]
        TypeControl <- controls[, 1:4]
        names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
        TypeControl$Type  <- toupper(TypeControl$Type)
        TypeControl <- as(TypeControl, "DataFrame")
    }
        
    
    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}

read.manifest.Allergy <- function(file) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 &&
                  is.integer(control.line) &&
                  !is.na(control.line))
    assay.line <- system(
        sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))
    stopifnot(length(assay.line) == 1 &&
                  is.integer(assay.line) &&
                  !is.na(assay.line))
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = FALSE,
        col.names = names(colClasses),
        sep = ",",
        comment.char = "",
        quote = "",
        skip = assay.line + 1L,
        colClasses = colClasses,
        nrows = control.line - assay.line - 2L)
    manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("IlmnID", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(1, 2, 3, 4, 5, 6, 7)] <- c("Name",
        "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[grep("^rs", TypeI$Name, invert=TRUE), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("IlmnID", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(1,2,3)] <- c("Name", "AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[grep("^rs", TypeII$Name, invert=TRUE), ]

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")

    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}


read.manifest.sesame.Allergy <- function(file) {
    ID_column <- "IlmnID"
    colNames <- readLines(file, n = 1)
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = TRUE,
        sep = ",",
        comment.char = "",
        colClasses = colClasses,
        quote = "")
    names(manifest)[names(manifest) == "Probe_ID"] <- "IlmnID"
    names(manifest)[names(manifest) == "U"] <- "AddressA_ID"
    names(manifest)[names(manifest) == "M"] <- "AddressB_ID"
    manifest$Infinium_Design_Type <- sub("1", "I", manifest$Infinium_Design_Type)
    manifest$Infinium_Design_Type <- sub("2", "II", manifest$Infinium_Design_Type)
    manifest$Color_Channel <- sub("Both", "", manifest$Color_Channel)
    controls.idx <- grep("^ctl_", manifest$IlmnID)
    nocontrols.idx <- grep("^ctl_", manifest$IlmnID, invert = TRUE)
    controls <- manifest[controls.idx,]
    manifest <- manifest[nocontrols.idx,]
    dupNames <- unique(manifest$Name[duplicated(manifest$Name)])
    manifest <- manifest[! manifest$Name %in% dupNames,]
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c(ID_column, "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI) <- c("Name", "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI$NextBase[is.na(TypeI$NextBase)] <- ""
    TypeI$Color[is.na(TypeI$Color)] <- ""
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[grep("^rs", TypeI$Name, invert = TRUE), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c(ID_column, "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII) <- c("Name", "AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[grep("^rs", TypeII$Name, invert=TRUE), ]

    TypeControl <- as(controls[, c("AddressA_ID", "Probe_Type", "IlmnID")], "DataFrame")
    names(TypeControl) <- c("Address", "Type", "ExtendedType")
    TypeControl$Type <- sub("BISULFITE_CONVERSION_", "BISULFITE CONVERSION ", TypeControl$Type)
    TypeControl$Type <- sub("SPECIFICITY_", "SPECIFICITY ", TypeControl$Type)
    TypeControl$Type <- sub("TARGET_REMOVAL", "TARGET REMOVAL", TypeControl$Type)
    TypeControl$ExtendedType <- sub("^ctl_", "", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Biotin_5K", "Biotin(5K)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Biotin_Bkg", "Biotin (Bkg)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Biotin_High", "Biotin (High)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("BS_Conversion_", "BS Conversion ", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("I_", "I-", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("DNP_20K", "DNP(20K)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("DNP_Bkg", "DNP (Bkg)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("DNP_High", "DNP (High)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Extension_A", "Extension (A)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Extension_C", "Extension (C)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Extension_G", "Extension (G)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Extension_T", "Extension (T)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("GT_Mismatch_", "GT Mismatch ", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("_MM", " (MM)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("_PM", " (PM)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Hyb_High", "Hyb (High)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Hyb_Low", "Hyb (Low)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Hyb_Medium", "Hyb (Medium)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_A", "NP (A)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_C", "NP (C)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G_1", "NP (G) 1", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G_2", "NP (G) 2", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G_3", "NP (G) 3", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G_4", "NP (G) 4", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G_5", "NP (G) 5", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_G", "NP (G)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("NP_T", "NP (T)", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Negative_", "Negative ", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Specificity_", "Specificity ", TypeControl$ExtendedType)
    TypeControl$ExtendedType <- sub("Target_Removal_", "Target Removal ", TypeControl$ExtendedType)
    
    return(list(manifestList = list(
                    TypeI = TypeI,
                    TypeII = TypeII,
                    TypeControl = TypeControl,
                    TypeSnpI = TypeSnpI,
                    TypeSnpII = TypeSnpII),
                manifest = manifest, controls = controls))
}


read.manifest.EPIC <- function(file) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 &&
                  is.integer(control.line) &&
                  !is.na(control.line))
    assay.line <- system(
        sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))
    stopifnot(length(assay.line) == 1 &&
                  is.integer(assay.line) &&
                  !is.na(assay.line))
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = FALSE,
        col.names = names(colClasses),
        sep = ",",
        comment.char = "",
        quote = "",
        skip = assay.line + 1L,
        colClasses = colClasses,
        nrows = control.line - assay.line - 2L)
    manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c(
        "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[-grep("^rs", TypeI$Name), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[-grep("^rs", TypeII$Name), ]

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")

    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}

read.manifest.450k <- function(file) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 &&
                  is.integer(control.line) &&
                  !is.na(control.line))
    assay.line <- system(
        sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))
    stopifnot(length(assay.line) == 1 &&
                  is.integer(assay.line) &&
                  !is.na(assay.line))

    # NOTE: Column headers is in line 8, hardcoded
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = TRUE,
        sep = ",",
        comment.char = "",
        quote = "",
        skip = 7,
        colClasses = colClasses,
        nrows = control.line - 9)
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2, 3, 4, 5, 6 , 7)] <-
        c("AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[-grep("^rs", TypeI$Name), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2, 3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[-grep("^rs", TypeII$Name), ]

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")

    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}


read.manifest.27k <- function(file) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -a -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    assay.line <- system(
        sprintf("grep -a -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))

    # NOTE: Column headers is in line 8, hardcoded
    colNames <- tail(readLines(file, n = assay.line + 1), n = 1)
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) = colNames
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = TRUE,
        sep = ",",
        comment.char = "",
        quote = "",
        skip = assay.line,
        colClasses = colClasses,
        nrows = control.line - (assay.line + 1),
        fill = TRUE)
    TypeI <- manifest[
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    TypeI <- TypeI[TypeI$Name != "", ]
    names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c(
        "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2, 3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- BStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeq, letters = "R"))

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
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
        DataFrame(
            Name = rs,
            AddressA = snp[grep("_A", snp$ExtendedType), "Address"],
            AddressB = snp[grep("_B", snp$ExtendedType), "Address"],
            Color = "Unknown")
    })
    TypeSnpI <- do.call(rbind, snps.sp)
    TypeSnpII <- TypeSnpI[0, ]

    list(manifestList =
             list(TypeI = TypeI,
                  TypeII = TypeII,
                  TypeControl = TypeControl,
                  TypeSnpI = TypeSnpI,
                  TypeSnpII = TypeSnpII),
         manifest = manifest,
         controls = controls)
}

# TODOs ------------------------------------------------------------------------

# TODO: Lots of duplicated code; DRY
