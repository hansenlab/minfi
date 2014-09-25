######################################################
## Functional normalization of the 450k array
## Jean-Philippe Fortin 
## Sept 24 2013
#####################################################

## Return a normalized beta matrix
##
## need stuff from gmodels.  Where?
preprocessFunnorm <- function(rgSet, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE) {
    .isRG(rgSet)
    rgSet <- updateObject(rgSet) ## FIXM: might not KDH: technically, this should not be needed, but might be nice

    # Background correction and dye bias normalization:
    if (bgCorr){
        if(verbose && dyeCorr) {
            cat("[preprocessFunnorm] Background and dye bias correction with noob \n") 
        } else {
            cat("[preprocessFunnorm] Background correction with noob \n") 
        }
        gmSet <- preprocessNoob(rgSet, dyeCorr = dyeCorr)
        if(verbose) cat("[preprocessFunnorm] Mapping to genome\n")
        gmSet <- mapToGenome(gmSet)
    } else {
        if(verbose) cat("[preprocessFunnorm] Mapping to genome\n")
        gmSet <- mapToGenome(rgSet)
    }
  
    subverbose <- max(as.integer(verbose) - 1L, 0)
    
    if(verbose) cat("[preprocessFunnorm] Quantile extraction\n")
    extractedData <- .extractFromRGSet450k(rgSet)
    rm(rgSet)

    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    if(verbose) cat("[preprocessFunnorm] Normalization\n")
    CN <- getCN(gmSet)
    gmSet <- .normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                   sex = sex, nPCs = nPCs, verbose = subverbose)
    grSet <- ratioConvert(gmSet, type = "Illumina")
    assay(grSet, "CN") <- CN
    grSet@preprocessMethod <- c(preprocessMethod(gmSet),
                                mu.norm = sprintf("Funnorm, nPCs=%s", nPCs))
    return(grSet)
 }	

 .getFunnormIndices <- function(object) {
     .isGenomic(object)
     probeType <- getProbeType(object, withColor = TRUE)
     autosomal <- (seqnames(object) %in% paste0("chr", 1:22))
     indices <- list(IGrn = which(probeType == "IGrn" & autosomal),
                     IRed = which(probeType == "IRed" & autosomal),
                     II = which(probeType == "II" & autosomal),
                     X = which(seqnames(object) == "chrX"),
                     Y = which(seqnames(object) == "chrY"))
     indices
 }

 .normalizeFunnorm450k <- function(object, extractedData, nPCs, sex, verbose = TRUE) {
     normalizeQuantiles <- function(matrix, indices, sex = NULL) {
         matrix <- matrix[indices,,drop=FALSE]
         ## uses probs, model.matrix, nPCS, through scoping)
         oldQuantiles <- t(colQuantiles(matrix, probs = probs))
         if(is.null(sex)) {
             newQuantiles <- .returnFit(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs)
         } else {
             newQuantiles <- .returnFitBySex(controlMatrix = model.matrix, quantiles = oldQuantiles, nPCs = nPCs, sex = sex)
         }
         .normalizeMatrix(matrix, newQuantiles)
     }

    indicesList <- .getFunnormIndices(object) 
    model.matrix <- .buildControlMatrix450k(extractedData)
    probs <- seq(from = 0, to = 1, length.out = 500)
    Meth <- getMeth(object)
    Unmeth <- getUnmeth(object)
    for (type in c("IGrn", "IRed", "II")) {
        indices <- indicesList[[type]]
        if(length(indices) > 0) {
            if(verbose) cat(sprintf("[normalizeFunnorm450k] Normalization of the %s probes\n", type))
            Unmeth[indices,] <- normalizeQuantiles(Unmeth, indices = indices, sex = NULL)
            Meth[indices,] <- normalizeQuantiles(Meth, indices = indices, sex = NULL)
        }
    }

    indices <- indicesList[["X"]]
    if(length(indices) > 0) {
        if(verbose) cat("[normalizeFunnorm450k] Normalization of the X-chromosome")
        Unmeth[indices,] <- normalizeQuantiles(Unmeth, indices = indices, sex = sex)
        Meth[indices,] <- normalizeQuantiles(Meth, indices = indices, sex = sex)
    }

    indices <- indicesList[["Y"]]
    if(length(indices) > 0) {
        if(verbose) cat("[normalizeFunnorm450k] Normalization of the Y-chromosome")
        sex <- as.character(sex)
        levels <- unique(sex)
        nSexes <- length(levels)
        if (nSexes == 2) {
            level1 <- levels[1]
            level2 <- levels[2]
        }
        if (nSexes == 2) {
            if (sum(sex == level1)>1) {
                Meth[indices, sex==level1]   <- preprocessCore::normalize.quantiles(Meth[indices, sex == level1, drop=FALSE])
                Unmeth[indices, sex==level1] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level1,drop=FALSE])
            }
            if (sum(sex == level2)>1) {
                Meth[indices, sex==level2]   <- preprocessCore::normalize.quantiles(Meth[indices, sex == level2,drop=FALSE])
                Unmeth[indices, sex==level2] <- preprocessCore::normalize.quantiles(Unmeth[indices, sex == level2,drop=FALSE])
            }
        } else {
            Meth[indices,] <- preprocessCore::normalize.quantiles(Meth[indices,])
            Unmeth[indices,] <- preprocessCore::normalize.quantiles(Unmeth[indices,])
        }
    }
    assay(object, "Meth") <- Meth
    assay(object, "Unmeth") <- Unmeth
    return(object)
}



### To extract quantiles and control probes from rgSet
.extractFromRGSet450k <- function(rgSet) {
    rgSet <- updateObject(rgSet)
    controlType <- c("BISULFITE CONVERSION I",
                     "BISULFITE CONVERSION II",
                     "EXTENSION",
                     "HYBRIDIZATION",
                     "NEGATIVE",
                     "NON-POLYMORPHIC",
                     "NORM_A",
                     "NORM_C",
                     "NORM_G",
                     "NORM_T",
                     "SPECIFICITY I",
                     "SPECIFICITY II",
                     "TARGET REMOVAL",
                     "STAINING")

    controlAddr <- getControlAddress(rgSet, controlType = controlType, asList = TRUE)
    controlAddr <- lapply(controlAddr, function(addr) {
        addr[addr %in% featureNames(rgSet)]
    })

    ## New code
    ctrls <- getProbeInfo(rgSet, type = "Control")
    ## FIXME: should be fixed in the manifest object
    ctrls <- ctrls[ctrls$Address %in% featureNames(rgSet),]
    ctrlsList <- split(ctrls, ctrls$Type)[controlType]
    ## End new code
    
    redControls <- getRed(rgSet)[ctrls$Address,]
    redControls <- lapply(ctrlsList, function(ctl) redControls[ctl$Address,])
    greenControls <- getGreen(rgSet)[ctrls$Address,]
    greenControls <- lapply(ctrlsList, function(ctl) greenControls[ctl$Address,])
    
    ## Extraction of the undefined negative control probes
    oobRaw <- getOOB(rgSet)
    probs <- c(0.01, 0.50, 0.99)
    greenOOB <- t(colQuantiles(oobRaw$Grn, na.rm = TRUE, probs = probs))
    redOOB   <- t(colQuantiles(oobRaw$Red, na.rm=TRUE,  probs = probs))
    oob      <- list(greenOOB = greenOOB, redOOB = redOOB)                       
    
    return(list(
        greenControls = greenControls,
        redControls = redControls,
        oob = oob, ctrlsList = ctrlsList))
}


## Extraction of the Control matrix
.buildControlMatrix450k <- function(extractedData) {
    getCtrlsAddr <- function(exType, index) {
        ctrls <- ctrlsList[[index]]
        addr <- ctrls$Address
        names(addr) <- ctrls$ExtendedType
        addr[exType]
    }
    
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)
    ctrlsList <- extractedData$ctrlsList
        
    ## Bisulfite conversion extraction for probe type II:
    index <- match("BISULFITE CONVERSION II", controlNames)
    redControls.current <- redControls[[ index ]]
    bisulfite2 <- colMeans(redControls.current, na.rm = TRUE)
    
    ## Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3), index = index)
    greenControls.current <- redControls[[ index ]][addr,]
    addr <- getCtrlsAddr(exType = sprintf("BS Conversion I-C%s", 4:6), index = index)
    redControls.current <- redControls[[ index ]][addr,]
    bisulfite1 <- colMeans(redControls.current + greenControls.current, na.rm = TRUE)
    
    ## Staining
    index <- match("STAINING", controlNames)
    addr <- getCtrlsAddr(exType = "Biotin (High)", index = index)
    stain.green <- greenControls[[ index ]][addr,]
    addr <- getCtrlsAddr(exType = "DNP (High)", index = index)
    stain.red <- redControls[[ index ]][addr, ] 
    
    ## Extension
    index <-    match("EXTENSION", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("A", "T")), index = index)
    extension.red <- t(redControls[[index]][addr,])
    colnames(extension.red) <- paste0("extRed", 1:ncol(extension.red))
    addr <- getCtrlsAddr(exType = sprintf("Extension (%s)", c("C", "G")), index = index)
    extension.green <- t(greenControls[[index]][addr,])
    colnames(extension.green) <- paste0("extGrn", 1:ncol(extension.green))
    
    ## Hybridization should be monitored only in the green channel
    index <- match("HYBRIDIZATION", controlNames)
    hybe <- t(greenControls[[index]])
    colnames(hybe) <- paste0("hybe", 1:ncol(hybe))
    
    ## Target removal should be low compared to hybridization probes
    index <- match("TARGET REMOVAL", controlNames)
    targetrem <- t(greenControls[[index]])
    colnames(targetrem) <- paste0("targetrem", 1:ncol(targetrem))
    
    ## Non-polymorphic probes
    index <- match("NON-POLYMORPHIC", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("A", "T")), index = index)
    nonpoly.red <- t(redControls[[index]][addr, ])
    colnames(nonpoly.red) <- paste0("nonpolyRed", 1:ncol(nonpoly.red))
    addr <- getCtrlsAddr(exType = sprintf("NP (%s)", c("C", "G")), index = index)
    nonpoly.green <- t(greenControls[[index]][addr, ])
    colnames(nonpoly.green) <- paste0("nonpolyGrn", 1:ncol(nonpoly.green))
    
    ## Specificity II
    index <- match("SPECIFICITY II", controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    spec2.green <- t(greenControls.current)
    colnames(spec2.green) <- paste0("spec2Grn", 1:ncol(spec2.green))
    spec2.red <- t(redControls.current)
    colnames(spec2.red) <- paste0("spec2Red", 1:ncol(spec2.red))
    spec2.ratio <- colMeans(greenControls.current, na.rm = TRUE) /
        colMeans(redControls.current, na.rm = TRUE)
    
    ## Specificity I
    index <- match("SPECIFICITY I", controlNames)
    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 1:3), index = index)
    greenControls.current <- greenControls[[index]][addr,]
    redControls.current <- redControls[[index]][addr,]
    spec1.green <- t(greenControls.current)
    colnames(spec1.green) <- paste0("spec1Grn", 1:ncol(spec1.green))
    spec1.ratio1 <- colMeans(redControls.current, na.rm = TRUE) /
        colMeans(greenControls.current, na.rm = TRUE)
    

    addr <- getCtrlsAddr(exType = sprintf("GT Mismatch %s (PM)", 4:6), index = index)
    greenControls.current <- greenControls[[index]][addr,]
    redControls.current <- redControls[[index]][addr,]
    spec1.red <- t(redControls.current)
    colnames(spec1.red) <- paste0("spec1Red", 1:ncol(spec1.red))
    spec1.ratio2 <- colMeans(greenControls.current, na.rm = TRUE) / 
        colMeans(redControls.current, na.rm = TRUE)
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2) / 2
    
    ## Normalization probes:
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans(greenControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans(greenControls[[index]], na.rm = TRUE)
    
    dyebias <- (normC + normG)/(normA + normT)

    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    oob.ratio <- oobG[2,]/oobR[2,]
    oobG <- t(oobG)
    colnames(oobG) <- paste0("oob", c(1,50,99))
    
    model.matrix <- cbind(
        bisulfite1, bisulfite2, extension.green, extension.red, hybe,
        stain.green, stain.red, nonpoly.green, nonpoly.red,
        targetrem, spec1.green, spec1.red, spec2.green, spec2.red, spec1.ratio1,
        spec1.ratio, spec2.ratio, spec1.ratio2, normA, normC, normT, normG, dyebias,
        oobG, oob.ratio)
    
    ## Imputation
    for (colindex in 1:ncol(model.matrix)) {
        if(any(is.na(model.matrix[,colindex]))) {
            column <- model.matrix[,colindex]
            column[is.na(column)]    <- mean(column, na.rm = TRUE)
            model.matrix[ , colindex] <- column
        }
    }
    
    ## Scaling   
    model.matrix <- scale(model.matrix)   
    
    ## Fixing outliers 
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3   
    
    ## Rescaling
    model.matrix <- scale(model.matrix) 
    
    return(model.matrix)
}


### Return the normalized quantile functions
.returnFit <- function(controlMatrix, quantiles, nPCs) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    ## Fixing potential problems with extreme quantiles
    quantiles[1,] <- 0
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles) - 1,] + 1000
    meanFunction <- rowMeans(quantiles)
    res <- quantiles - meanFunction
    controlPCs <- prcomp(controlMatrix)$x[,1:nPCs,drop=FALSE]
    design <- model.matrix(~controlPCs - 1)
    fits <- lm.fit(x = design, y = t(res))
    newQuantiles <- meanFunction + t(fits$residuals)
    return(newQuantiles)
}

.returnFitBySex <- function(controlMatrix, quantiles, nPCs, sex) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(controlMatrix))
    stopifnot(ncol(quantiles) == nrow(controlMatrix))
    sex    <- as.character(sex)
    levels <- unique(sex)
    nSexes <- length(levels)
    if (nSexes == 2) {
        sex1 <- sum(sex == levels[1])
        sex2 <- sum(sex == levels[2])
        
    } else {
        sex1 <- sum(sex == levels[1])
        sex2 <- 0
    }
    
    ## When normalization should not be performed by sex separately:
    if ((sex1 <= 10) | (sex2 <= 10)) {
        newQuantiles <- .returnFit(controlMatrix = controlMatrix,
                                  quantiles = quantiles, 
                                  nPCs = nPCs)
    } else {
        quantiles1 <- quantiles[, sex == levels[1]]
        controlMatrix1 <- controlMatrix[sex == levels[1], ]
        
        newQuantiles1 <- .returnFit(controlMatrix = controlMatrix1,
                                   quantiles = quantiles1, 
                                   nPCs = nPCs)
        
        quantiles2 <- quantiles[, sex == levels[2]]
        controlMatrix2 <- controlMatrix[sex == levels[2], ]
        
        newQuantiles2 <- .returnFit(controlMatrix = controlMatrix2,
                                   quantiles = quantiles2, 
                                   nPCs = nPCs)
        
        newQuantiles <- quantiles
        newQuantiles[, sex == levels[1]] <- newQuantiles1
        newQuantiles[, sex == levels[2]] <- newQuantiles2
    }
    
    return(newQuantiles)
}

### Normalize a matrix of intensities  
.normalizeMatrix <- function(intMatrix, newQuantiles) {
    ## normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix)) 
    n <- nrow(newQuantiles)
    normMatrix <- sapply(1:ncol(intMatrix), function(i) {
        crtColumn <- intMatrix[ , i]
        crtColumn.reduced <- crtColumn[!is.na(crtColumn)]
        ## Generation of the corrected intensities:
        target <- sapply(1:(n-1), function(j) {
            start <- newQuantiles[j,i]
            end <- newQuantiles[j+1,i]
            sequence <- seq(start, end,( end-start)/n)[-n]
            return(sequence)
        })
        target <- as.vector(target)
        result <- preprocessCore::normalize.quantiles.use.target(matrix(crtColumn.reduced), target)
        return(result)
    })
    return(normMatrix)
}

