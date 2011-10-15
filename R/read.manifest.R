read.manifest <- function(file, returnAll = FALSE) {
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
    TypeII <- manifest[manifest$Infinium_Design_Type == "II",
                        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("Address", "ProbeSeq")
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
