# A helper function to check whether two SummarizedExperiment objects, one
# matrix-backed and one DelayedMatrix-backed, are 'equal'
checkEqualSummarizedExperiments <- function(SE, SE_with_DA) {
    stopifnot(is(SE, "SummarizedExperiment"),
              is(SE_with_DA, "SummarizedExperiment"),
              is(assay(SE), "matrix"),
              is(assay(SE_with_DA), "DelayedMatrix"))
    assays(SE_with_DA) <- endoapply(assays(SE_with_DA), as.matrix)
    all.equal(SE, SE_with_DA)
}

test_bumphunter_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    gmSet <- preprocessQuantile(MsetEx)
    design <- model.matrix(~ gmSet$status)
    bumps <- bumphunter(gmSet, design = design, B = 0,
                        type = "Beta", cutoff = 0.25)
    bumps0 <- bumps

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    gmSet <- preprocessQuantile(MsetEx)
    design <- model.matrix(~ gmSet$status)
    bumps <- bumphunter(gmSet, design = design, B = 0,
                        type = "Beta", cutoff = 0.25)

    # Check equivalency
    checkEquals(bumps0, bumps)
}

test_combineArrays_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(minfiDataEPIC))
    stopifnot(require(DelayedArray))

    # Original example
    data(RGsetEx.sub)
    data(RGsetEPIC)
    rgSet <- combineArrays(RGsetEPIC, RGsetEx.sub)
    rgSet0 <- rgSet

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    assays(RGsetEPIC) <- endoapply(assays(RGsetEPIC), DelayedArray)
    rgSet <- combineArrays(RGsetEPIC, RGsetEx.sub)

    # Check equivalency
    checkEqualSummarizedExperiments(rgSet0, rgSet)
}

test_compartments_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    GMset <- mapToGenome(MsetEx)
    comps <- compartments(GMset, res = 10^6)
    comps0 <- comps

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    GMset <- mapToGenome(MsetEx)
    comps <- compartments(GMset, res = 10^6)

    # Check equivalency
    checkEquals(comps0, comps)
}

test_controlStripPlot_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stop("Currently no way to test that the 2 plots are identical",
         call. = TRUE)
    # TODO: Conceptually, the below is what is needed

    # Original example
    names <- pData(RGsetEx)$Sample_Name
    plot0 <- controlStripPlot(RGsetEx, controls=c("BISULFITE CONVERSION I"),
                              sampNames=names)


    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    names <- pData(RGsetEx)$Sample_Name
    plot <- controlStripPlot(RGsetEx, controls=c("BISULFITE CONVERSION I"),
                             sampNames=names)

    # Check equivalency
    checkIdentical(plot0, plot)
}

test_convertArray_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    data(RGsetEx.sub)
    rgSet <- convertArray(RGsetEx.sub, outType = "IlluminaHumanMethylationEPIC")
    rgSet0 <- rgSet

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    rgSet <- convertArray(RGsetEx.sub, outType = "IlluminaHumanMethylationEPIC")

    # Check equivalency
    checkEqualSummarizedExperiments(rgSet0, rgSet)
}

test_densityBeanPlot_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stop("Currently no way to test that the 2 plots are identical",
         call. = TRUE)
    # TODO: Conceptually, the below is what is needed

    # Original example
    names <- pData(RGsetEx)$Sample_Name
    groups <- pData(RGsetEx)$Sample_Group
    par(mar=c(5,6,4,2))
    plot0 <- densityBeanPlot(RGsetEx, sampNames=names, sampGroups=groups)

    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    names <- pData(RGsetEx)$Sample_Name
    groups <- pData(RGsetEx)$Sample_Group
    par(mar=c(5,6,4,2))
    plot <- densityBeanPlot(RGsetEx, sampNames=names, sampGroups=groups)

    # Check equivalency
    checkIdentical(plot0, plot)
}

test_detectionP_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    detP <- detectionP(RGsetEx.sub)
    failed <- detP>0.01
    colMeans(failed) # Fraction of failed positions per sample
    sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
    failed0 <- failed

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    detP <- detectionP(RGsetEx.sub)
    failed <- detP>0.01
    colMeans(failed) # Fraction of failed positions per sample
    sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?

    # Check equivalency
    checkIdentical(failed0, failed)
}

test_dmpFinder_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    grp <- pData(MsetEx)$Sample_Group
    MsetExSmall <- MsetEx[1:1e4,] # To speed up the example
    M <- getM(MsetExSmall, type = "beta", betaThreshold = 0.001)
    dmp <- dmpFinder(M, pheno=grp, type="categorical")
    sum(dmp$qval < 0.05, na.rm=TRUE)
    head(dmp)
    dmp0 <- dmp

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    grp <- pData(MsetEx)$Sample_Group
    MsetExSmall <- MsetEx[1:1e4,] # To speed up the example
    M <- getM(MsetExSmall, type = "beta", betaThreshold = 0.001)
    dmp <- dmpFinder(M, pheno=grp, type="categorical")
    sum(dmp$qval < 0.05, na.rm=TRUE)
    head(dmp)

    # Check equivalency
    checkIdentical(dmp0, dmp)
}

test_estimateCellCounts_examples <- function() {
    stopifnot(require(FlowSorted.Blood.450k))
    stopifnot(require(DelayedArray))

    # Original example
    wh.WBC <- which(FlowSorted.Blood.450k$CellType == "WBC")
    wh.PBMC <- which(FlowSorted.Blood.450k$CellType == "PBMC")
    RGset <- FlowSorted.Blood.450k[, c(wh.WBC, wh.PBMC)]
    ## The following line is purely to work around an issue with repeated
    ## sampleNames and Biobase::combine()
    sampleNames(RGset) <- paste(RGset$CellType,
                                c(seq(along = wh.WBC), seq(along = wh.PBMC)), sep = "_")
    counts <- estimateCellCounts(RGset, meanPlot = FALSE)
    round(counts, 2)
    counts0 <- counts

    # DelayedArray-backed version
    assays(FlowSorted.Blood.450k) <- endoapply(assays(FlowSorted.Blood.450k), DelayedArray)
    wh.WBC <- which(FlowSorted.Blood.450k$CellType == "WBC")
    wh.PBMC <- which(FlowSorted.Blood.450k$CellType == "PBMC")
    RGset <- FlowSorted.Blood.450k[, c(wh.WBC, wh.PBMC)]
    ## The following line is purely to work around an issue with repeated
    ## sampleNames and Biobase::combine()
    sampleNames(RGset) <- paste(RGset$CellType,
                                c(seq(along = wh.WBC), seq(along = wh.PBMC)), sep = "_")
    counts <- estimateCellCounts(RGset, meanPlot = FALSE)
    round(counts, 2)

    # Check equivalency
    checkIdentical(counts0, counts)
}

test_fixMethOutliers_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    MsetEx <- fixMethOutliers(MsetEx)
    MsetEx0 <- MsetEx

    # DelayedArray-backed version
    MsetEx <- minfiData::MsetEx
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    MsetEx <- fixMethOutliers(MsetEx)

    # Check equivalency
    checkEqualSummarizedExperiments(MsetEx0, MsetEx)
}

test_gaphunter_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    gapres <- gaphunter(MsetEx.sub, threshold=0.3, keepOutliers=TRUE)
    gapres0 <- gapres

    # DelayedArray-backed version
    assays(MsetEx.sub) <- endoapply(assays(MsetEx.sub, DelayedArray))
    gapres <- gaphunter(MsetEx.sub, threshold=0.3, keepOutliers=TRUE)

    # Check equivalency
    checkEqualSummarizedExperiments(gapres0, gapres)
}

test_getAnnotation_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    table(getIslandStatus(MsetEx))
    anno0 <- getAnnotation(MsetEx, what = "Manifest")

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    table(getIslandStatus(MsetEx))
    anno <- getAnnotation(MsetEx, what = "Manifest")

    # Check equivalency
    checkIdentical(anno0, anno)
}

test_getQC_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    qc <- getQC(MsetEx)
    MsetEx <- addQC(MsetEx, qc = qc)
    ## plotQC(qc)
    MsetEx0 <- MsetEx

    # DelayedArray-backed version
    MsetEx <- minfiData::MsetEx
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    qc <- getQC(MsetEx)
    MsetEx <- addQC(MsetEx, qc = qc)
    ## plotQC(qc)

    # Check equivalency
    checkEqualSummarizedExperiments(MsetEx0, MsetEx)
}

test_getSex_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    GMsetEx <- mapToGenome(MsetEx)
    estSex <- getSex(GMsetEx)
    GMsetEx <- addSex(GMsetEx, sex = estSex)
    GMsetEx0 <- GMsetEx0


    # DelayedArray-backed version
    MsetEx <- minfiData::MsetEx
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    GMsetEx <- mapToGenome(MsetEx)
    estSex <- getSex(GMsetEx)
    GMsetEx <- addSex(GMsetEx, sex = estSex)

    # Check equivalency
    checkEqualSummarizedExperiments(GMsetEx0, GMsetEx)
}

test_logit2_examples <- function() {
    stopifnot(require(DelayedArray))

    # Original example
    x <- matrix(c(0.25, 0.5, 0.75))
    val0 <- logit2(x)

    # DelayedArray-backed version
    x <- DelayedArray(x)
    val <- logit2(x)

    # Check equivalency
    checkEquals(as.matrix(val0), as.matrix(val))
}

test_makeGenomicRatioSetFromMatrix_examples <- function() {
    stopifnot(require(DelayedArray))

    # Original example
    mat <- matrix(10,5,2)
    rownames(mat) <- c( "cg13869341", "cg14008030","cg12045430", "cg20826792","cg00381604")
    grset <- makeGenomicRatioSetFromMatrix(mat)
    grset0 <- grset

    # DelayedArray-backed version
    mat <- DelayedArray(matrix(10,5,2))
    rownames(mat) <- c( "cg13869341", "cg14008030","cg12045430", "cg20826792","cg00381604")
    grset <- makeGenomicRatioSetFromMatrix(mat)

    # Check equivalency
    checkEqualSummarizedExperiments(grset0, grset)
}

test_mapToGenome_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    GMsetEx.sub <- mapToGenome(MsetEx.sub)
    GMsetEx.sub0 <- GMsetEx.sub

    # DelayedArray-backed version
    assays(MsetEx.sub) <- endoapply(assays(MsetEx.sub), DelayedArray)
    GMsetEx.sub <- mapToGenome(MsetEx.sub)

    checkEqualSummarizedExperiments(GMsetEx.sub0, GMsetEx.sub)
}

test_mdsPlot_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stop("Currently no way to test that the 2 plots are identical",
         call. = TRUE)
    # TODO: Conceptually, the below is what is needed

    # Original example
    names <- pData(MsetEx)$Sample_Name
    groups <- pData(MsetEx)$Sample_Group
    plot0 <- mdsPlot(MsetEx, sampNames=names, sampGroups=groups)

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    names <- pData(MsetEx)$Sample_Name
    groups <- pData(MsetEx)$Sample_Group
    plot <- mdsPlot(MsetEx, sampNames=names, sampGroups=groups)

    # Check equivalency
    checkIdentical(plot0, plot)
}

test_minfiQC_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    out <- minfiQC(MsetEx)
    ## plotQC(out$qc)
    ## plotSex(out$sex)
    out0 <- out

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    out <- minfiQC(MsetEx)
    ## plotQC(out$qc)
    ## plotSex(out$sex)

    # Check equivalency
    checkEquals(out0, out)
}

test_plotBetasByType_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stop("Currently no way to test that the 2 plots are identical",
         call. = TRUE)
    # TODO: Conceptually, the below is what is needed

    # Original example
    Mset.swan <- preprocessSWAN(RGsetEx, MsetEx)
    par(mfrow=c(1,2))
    plot_raw0 <- plotBetasByType(MsetEx[,1], main="Raw")
    plot_swan0 <- plotBetasByType(Mset.swan[,1], main="SWAN")


    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    Mset.swan <- preprocessSWAN(RGsetEx, MsetEx)
    par(mfrow=c(1,2))
    plot_raw <- plotBetasByType(MsetEx[,1], main="Raw")
    plot_swan <- plotBetasByType(Mset.swan[,1], main="SWAN")

    # Check equivalency
    checkEquals(plot_raw0, plot_raw)
    checkEquals(plot_swan0, plot_swan)
}

test_plotCpg_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stop("Currently no way to test that the 2 plots are identical",
         call. = TRUE)
    # TODO: Conceptually, the below is what is needed

    # Original example
    grp <- pData(MsetEx)$Sample_Group
    cpgs <- c("cg00050873", "cg00212031", "cg26684946", "cg00128718")
    par(mfrow=c(2,2))
    plot0 <- plotCpg(MsetEx, cpg=cpgs, pheno=grp, type="categorical")
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)

    # DelayedArray-backed version
    grp <- pData(MsetEx)$Sample_Group
    cpgs <- c("cg00050873", "cg00212031", "cg26684946", "cg00128718")
    par(mfrow=c(2,2))
    plot <- plotCpg(MsetEx, cpg=cpgs, pheno=grp, type="categorical")

    # Check equivalency
    checkEquals(plot0, plot)
}

test_preprocessFunnorm_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    Mset.sub.funnorm <- preprocessFunnorm(RGsetEx.sub)
    Mset.sub.funnorm0 <- Mset.sub.funnorm

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- assays(RGsetEx.sub, DelayedArray)
    Mset.sub.funnorm <- preprocessFunnorm(RGsetEx.sub)

    # Check equivalency
    checkEqualSummarizedExperiments(Mset.sub.funnorm0, Mset.sub.funnorm)
}

test_preprocessIllumina_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    dat <- preprocessIllumina(RGsetEx, bg.correct=FALSE, normalize="controls")
    slot(name="preprocessMethod", dat)[1]
    dat0 <- dat

    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    dat <- preprocessIllumina(RGsetEx, bg.correct=FALSE, normalize="controls")
    slot(name="preprocessMethod", dat)[1]

    # Check equality
    checkEqualSummarizedExperiments(dat0, dat)
}

test_preprocessNoob_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    MsetEx.sub.noob <- preprocessNoob(RGsetEx.sub)
    dyeMethods <- c(ssNoob="single", refNoob="reference")
    GRsets <- lapply(dyeMethods,
                     function(m) preprocessNoob(RGsetEx, dyeMethod=m))
    # NOTE: Was all.equal() in example, here changed to checkEquals()
    checkEquals(getBeta(GRsets$refNoob), getBeta(GRsets$ssNoob))
    GRsets0 <- GRsets

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    MsetEx.sub.noob <- preprocessNoob(RGsetEx.sub)
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    dyeMethods <- c(ssNoob="single", refNoob="reference")
    GRsets <- lapply(dyeMethods,
                     function(m) preprocessNoob(RGsetEx, dyeMethod=m))
    # NOTE: Was all.equal() in example, here changed to checkEquals()
    checkEquals(getBeta(GRsets$refNoob), getBeta(GRsets$ssNoob))

    # Check equivalency
    checkEqualSummarizedExperiments(GRsets0, GRsets)
}

test_preprocessQuantile_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    GMset.sub.quantile <- preprocessSWAN(RGsetEx.sub)
    GMset.sub.quantile0 <- GMset.sub.quantile
    GMset <- preprocessQuantile(RGsetEx)
    GMset0 <- GMset

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    GMset.sub.quantile <- preprocessSWAN(RGsetEx.sub)
    assays(RGset.sub) <- endoapply(assays(RGset.sub), DelayedArray)
    GMset <- preprocessQuantile(RGsetEx)

    # Check equivalency
    checkEqualSummarizedExperiments(GMset.sub.quantile0, GMset.sub.quantile)
    checkEqualSummarizedExperiments(GMset0, GMset)
}

test_preprocessRaw_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    dat <- preprocessRaw(RGsetEx)
    dat0 <- dat

    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    dat <- preprocessRaw(RGsetEx)

    # Check equivalency
    checkEqualSummarizedExperiments(dat0, dat)
}

test_preprocessSwan_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original example
    MsetEx.sub.swan <- preprocessSWAN(RGsetEx.sub)
    dat <- preprocessRaw(RGsetEx)
    preprocessMethod(dat)
    datSwan <- preprocessSWAN(RGsetEx, mSet = dat)
    datIlmn <- preprocessIllumina(RGsetEx)
    preprocessMethod(datIlmn)
    datIlmnSwan <- preprocessSWAN(RGsetEx, mSet = datIlmn)
    dat0 <- dat
    datSwan0 <- datSwan
    datIlmn0 <- datIlmn
    datIlmnSwan0 <- datIlmnSwan

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    MsetEx.sub.swan <- preprocessSWAN(RGsetEx.sub)
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    dat <- preprocessRaw(RGsetEx)
    preprocessMethod(dat)
    datSwan <- preprocessSWAN(RGsetEx, mSet = dat)
    datIlmn <- preprocessIllumina(RGsetEx)
    preprocessMethod(datIlmn)
    datIlmnSwan <- preprocessSWAN(RGsetEx, mSet = datIlmn)

    # Check equivalency
    checkEqualSummarizedExperiments(dat0, dat)
    checkEqualSummarizedExperiments(datSwan0, datSwan)
    checkEqualSummarizedExperiments(datIlmn0, datIlmn)
    checkEqualSummarizedExperiments(datIlmnSwan0, datIlmnSwan)
}

test_qcReport_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))
    stopifnot(require(tools))

    tmpdir <- tempdir()

    # Original examples
    names <- pData(RGsetEx)$Sample_Name
    groups <- pData(RGsetEx)$Sample_Group
    qcReport(RGsetEx, sampNames=names, sampGroups=groups,
             pdf=file.path(tmpdir, "qcReport0.pdf"))

    # DelayedArray-backed version
    assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
    names <- pData(RGsetEx)$Sample_Name
    groups <- pData(RGsetEx)$Sample_Group
    qcReport(RGsetEx, sampNames=names, sampGroups=groups,
             pdf=file.path(tmpdir, "qcReport.pdf"))

    # Check equivalency
    checkIdentical(md5sum(file.path(tmpdir, "qcReport0.pdf")),
                   md5sum(file.path(tmpdir, "qcReport.pdf")))
}

test_ratioConvert_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original examples
    RsetEx.sub <- ratioConvert(MsetEx.sub, keepCN = TRUE)
    RsetEx.sub0 <- RsetEx.sub

    # DelayedArray-backed version
    assays(MsetEx.sub) <- endoapply(assays(MsetEx.sub), DelayedArray)
    RsetEx.sub <- ratioConvert(MsetEx.sub, keepCN = TRUE)

    # Check equivalency
    checkEqualSummarizedExperiments(RGsetEx.sub0, RGsetEx.sub)
}

test_subsetByLoci_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original examples
    loci <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611")
    a0 <- subsetByLoci(RGsetEx.sub, includeLoci = loci)
    b0 <- subsetByLoci(RGsetEx.sub, excludeLoci = loci)

    # DelayedArray-backed version
    assays(RGsetEx.sub) <- endoapply(assays(RGsetEx.sub), DelayedArray)
    a <- subsetByLoci(RGsetEx.sub, includeLoci = loci)
    b <- subsetByLoci(RGsetEx.sub, excludeLoci = loci)

    # Check equivalency
    checkEqualSummarizedExperiments(a0, a)
    checkEqualSummarizedExperiments(b0, b)
}

test_utils_examples <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(DelayedArray))

    # Original examples
    a0 <- head(getMethSignal(MsetEx, what = "Beta"))

    # DelayedArray-backed version
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    a <- head(getMethSignal(MsetEx, what = "Beta"))

    # Check equivalency
    checkEquals(a0, a)
}
