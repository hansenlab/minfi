## # A helper function to check whether two SummarizedExperiment objects are
## # equivalent
## checkEquivalentSummarizedExperiments <- function(SE1, SE2) {
##     stopifnot(require(SummarizedExperiment))
##     stopifnot(is(SE1, "SummarizedExperiment"),
##               is(SE2, "SummarizedExperiment"))
##     assays(SE1) <- endoapply(assays(SE1), as.matrix)
##     assays(SE2) <- endoapply(assays(SE2), as.matrix)
##     if (isTRUE(all.equal(SE1, SE2))) {
##         return(all.equal(assays(SE1), assays(SE2)))
##     }
##     FALSE
## }

## # TODO: Uncomment once bumphunter() supports DelayedArray-backed minfi
## #       objects
## test_bumphunter_examples <- function() {

##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))
##     # Original example
##     # gmSet <- preprocessQuantile(MsetEx)
##     # design <- model.matrix(~ gmSet$status)
##     # bumps <- bumphunter(
##     #     object = gmSet,
##     #     design = design,
##     #     B = 0,
##     #     type = "Beta",
##     #     cutoff = 0.25)
##     # bumps0 <- bumps

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     checkException(gmSet <- preprocessQuantile(MsetEx),
##                    silent = TRUE)
##     # design <- model.matrix(~ gmSet$status)
##     # checkException(bumps <- bumphunter(
##     #     object = gmSet,
##     #     design = design,
##     #     B = 0,
##     #     type = "Beta",
##     #     cutoff = 0.25))

##     # Check equivalency
##     # checkIdentical(bumps0, bumps)
## }

## test_combineArrays_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(minfiDataEPIC))
##     stopifnot(require(DelayedArray))

##     # Original example
##     data(RGsetEx.sub)
##     data(RGsetEPIC)
##     rgSet0 <- combineArrays(RGsetEPIC, RGsetEx.sub)

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     RGsetEPIC <- realize(RGsetEPIC)
##     rgSet <- combineArrays(RGsetEPIC, RGsetEx.sub)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(rgSet0, rgSet)
## }

## # TODO: Uncomment once compartments() supports DelayedArray-backed minfi objects
## test_compartments_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     # GMset <- mapToGenome(MsetEx)
##     # comps <- compartments(GMset, res = 10^6)
##     # comps0 <- comps

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     GMset <- mapToGenome(MsetEx)
##     checkException(comps <- compartments(GMset, res = 10^6), silent = TRUE)

##     # Check equivalency
##     # checkEquals(comps0, comps)
## }

## test_convertArray_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     data(RGsetEx.sub)
##     rgSet <- convertArray(RGsetEx.sub, outType = "IlluminaHumanMethylationEPIC")
##     rgSet0 <- rgSet

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     rgSet <- convertArray(RGsetEx.sub, outType = "IlluminaHumanMethylationEPIC")

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(rgSet0, rgSet)
## }

## test_detectionP_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     detP <- detectionP(RGsetEx.sub)
##     failed0 <- detP > 0.01

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     detP <- detectionP(RGsetEx.sub)
##     failed <- detP > 0.01

##     # Check equivalency
##     # NOTE: detectionP() returns a DelayedMatrix when applied to a
##     #       DelayedArray-backed object, hence need for as.matrix()
##     checkIdentical(failed0, as.matrix(failed))
## }

## # TODO: Uncomment once dmpFinder() supports DelayedArray-backed minfi objects
## test_dmpFinder_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     MsetExSmall <- MsetEx[seq_len(1e4), ]
##     # grp <- pData(MsetEx)$Sample_Group
##     # M <- getM(MsetExSmall, type = "beta", betaThreshold = 0.001)
##     # dmp0 <- dmpFinder(M, pheno = grp, type = "categorical")

##     # DelayedArray-backed version
##     MsetExSmall <- realize(MsetExSmall)
##     M <- getM(MsetExSmall, type = "beta", betaThreshold = 0.001)
##     checkException(dmp <- dmpFinder(M, pheno = grp, type = "categorical"),
##                    silent = TRUE)

##     # Check equivalency
##     # checkIdentical(dmp0, dmp)
## }

## # TODO: Uncomment once estimateCellCounts() supports DelayedArray-backed minfi
## #       objects
## test_estimateCellCounts_examples <- function() {
##     stopifnot(require(FlowSorted.Blood.450k))
##     stopifnot(require(DelayedArray))

##     # Original example
##     wh.WBC <- which(FlowSorted.Blood.450k$CellType == "WBC")
##     wh.PBMC <- which(FlowSorted.Blood.450k$CellType == "PBMC")
##     RGset <- FlowSorted.Blood.450k[, c(wh.WBC, wh.PBMC)]
##     # The following line is purely to work around an issue with repeated
##     # sampleNames and Biobase::combine()
##     sampleNames(RGset) <- paste(
##         RGset$CellType,
##         c(seq(along = wh.WBC), seq(along = wh.PBMC)),
##         sep = "_")
##     # counts0 <- estimateCellCounts(RGset, meanPlot = FALSE)

##     # DelayedArray-backed version
##     RGset <- realize(RGset)
##     checkException(counts <- estimateCellCounts(RGset, meanPlot = FALSE),
##                    silent = TRUE)

##     # Check equivalency
##     # checkIdentical(counts0, counts)
## }

## test_fixMethOutliers_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     MsetEx0 <- fixMethOutliers(MsetEx)

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     MsetEx <- fixMethOutliers(MsetEx)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(MsetEx0, MsetEx)
## }

## # TODO: Uncomment once gaphunter() supports DelayedArray-backed minfi objects
## test_gaphunter_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     # gapres0 <- gaphunter(MsetEx.sub, threshold = 0.3, keepOutliers = TRUE)

##     # DelayedArray-backed version
##     MsetEx.sub <- realize(MsetEx.sub)
##     checkException(
##         gapres <- gaphunter(MsetEx.sub, threshold = 0.3, keepOutliers = TRUE),
##         silent = TRUE)

##     # Check equivalency
##     # checkEquivalentSummarizedExperiments(gapres0, gapres)
## }

## test_getAnnotation_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     anno0 <- getAnnotation(MsetEx, what = "Manifest")

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     anno <- getAnnotation(MsetEx, what = "Manifest")

##     # Check equivalency
##     checkIdentical(anno0, anno)
## }

## test_getQC_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     qc <- getQC(MsetEx)
##     MsetEx0 <- addQC(MsetEx, qc = qc)

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     qc <- getQC(MsetEx)
##     MsetEx <- addQC(MsetEx, qc = qc)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(MsetEx0, MsetEx)
## }

## test_getSex_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     GMsetEx <- mapToGenome(MsetEx)
##     estSex <- getSex(GMsetEx)
##     GMsetEx0 <- addSex(GMsetEx, sex = estSex)

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     GMsetEx <- mapToGenome(MsetEx)
##     estSex <- getSex(GMsetEx)
##     GMsetEx <- addSex(GMsetEx, sex = estSex)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(GMsetEx0, GMsetEx)
## }

## test_logit2_examples <- function() {
##     stopifnot(require(DelayedArray))

##     # Original example
##     x <- matrix(c(0.25, 0.5, 0.75))
##     val0 <- logit2(x)

##     # DelayedArray-backed version
##     x <- DelayedArray(x)
##     val <- logit2(x)

##     # Check equivalency
##     # NOTE: logi2() returns a DelayedMatrix when applied to a DelayedArra,
##     #       hence need for as.matrix()
##     checkIdentical(val0, as.matrix(val))
## }

## # TODO: Uncomment once makeGenomicRatioSetFromMatrix() supports DelayedArray
## test_makeGenomicRatioSetFromMatrix_examples <- function() {
##     stopifnot(require(DelayedArray))

##     # Original example
##     mat <- matrix(10, 5, 2)
##     rownames(mat) <- c(
##         "cg13869341", "cg14008030","cg12045430", "cg20826792","cg00381604")
##     # grset0 <- makeGenomicRatioSetFromMatrix(mat)

##     # DelayedArray-backed version
##     mat <- realize(mat)
##     checkException(grset <- makeGenomicRatioSetFromMatrix(mat), silent = TRUE)

##     # Check equivalency
##     # checkEquivalentSummarizedExperiments(grset0, grset)
## }

## test_mapToGenome_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     GMsetEx.sub0 <- mapToGenome(MsetEx.sub)

##     # DelayedArray-backed version
##     MsetEx.sub <- realize(MsetEx.sub)
##     GMsetEx.sub <- mapToGenome(MsetEx.sub)

##     checkEquivalentSummarizedExperiments(GMsetEx.sub0, GMsetEx.sub)
## }

## test_minfiQC_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     out0 <- minfiQC(MsetEx)

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     out <- minfiQC(MsetEx)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(out0$object, out$object)
##     checkIdentical(out0$qc, out$qc)
## }

## # TODO: Uncomment once preprocessFunnorm() supports DelayedArray-backed minfi
## #       objects
## test_preprocessFunnorm_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     # Mset.sub.funnorm0 <- preprocessFunnorm(RGsetEx.sub)

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     checkException(Mset.sub.funnorm <- preprocessFunnorm(RGsetEx.sub),
##                    silent = TRUE)

##     # Check equivalency
##     # checkEquivalentSummarizedExperiments(Mset.sub.funnorm0, Mset.sub.funnorm)
## }

## test_preprocessIllumina_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     dat0 <- preprocessIllumina(
##         rgSet = RGsetEx,
##         bg.correct = FALSE,
##         normalize = "controls")

##     # DelayedArray-backed version
##     RGsetEx <- realize(RGsetEx)
##     dat <- preprocessIllumina(
##         rgSet = RGsetEx,
##         bg.correct = FALSE,
##         normalize = "controls")

##     # Check equality
##     checkEquivalentSummarizedExperiments(dat0, dat)
## }

## test_preprocessNoob_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     MsetEx.sub.noob <- preprocessNoob(RGsetEx.sub)
##     dyeMethods <- c(ssNoob = "single", refNoob = "reference")
##     GRsets0 <- lapply(
##         X = dyeMethods,
##         FUN = function(m) preprocessNoob(RGsetEx, dyeMethod = m))

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     MsetEx.sub.noob <- preprocessNoob(RGsetEx.sub)
##     GRsets <- lapply(
##         X = dyeMethods,
##         FUN = function(m) preprocessNoob(RGsetEx, dyeMethod = m))

##     # Check equivalency
##     Map(checkEquivalentSummarizedExperiments, GRsets0, GRsets)
## }

## # TODO: Uncomment once preprocessQuantile() supports DelayedArray-backed minfi
## #       objects
## test_preprocessQuantile_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     # GMset.sub.quantile0 <- preprocessQuantile(RGsetEx.sub)
##     # GMset0 <- preprocessQuantile(RGsetEx)

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     checkException(GMset.sub.quantile <- preprocessQuantile(RGsetEx.sub),
##                    silent = TRUE)
##     RGsetEx <- realize(RGsetEx)
##     checkException(GMset <- preprocessQuantile(RGsetEx),
##                    silent = TRUE)

##     # Check equivalency
##     # checkEquivalentSummarizedExperiments(
##     #     GMset.sub.quantile0,
##     #     GMset.sub.quantile)
##     # checkEquivalentSummarizedExperiments(GMset0, GMset)
## }

## test_preprocessRaw_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     dat0 <- preprocessRaw(RGsetEx)

##     # DelayedArray-backed version
##     RGsetEx <- realize(RGsetEx)
##     dat <- preprocessRaw(RGsetEx)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(dat0, dat)
## }

## test_preprocessSwan_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original example
##     set.seed(666)
##     MsetEx.sub.swan0 <- preprocessSWAN(RGsetEx.sub)
##     dat0 <- preprocessRaw(RGsetEx)
##     set.seed(777)
##     datSwan0 <- preprocessSWAN(RGsetEx, mSet = dat0)
##     datIlmn0 <- preprocessIllumina(RGsetEx)
##     set.seed(888)
##     datIlmnSwan0 <- preprocessSWAN(RGsetEx, mSet = datIlmn0)

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     RGsetEx <- realize(RGsetEx)
##     set.seed(666)
##     MsetEx.sub.swan <- preprocessSWAN(RGsetEx.sub)
##     dat <- preprocessRaw(RGsetEx)
##     set.seed(777)
##     datSwan <- preprocessSWAN(RGsetEx, mSet = dat)
##     datIlmn <- preprocessIllumina(RGsetEx)
##     set.seed(888)
##     datIlmnSwan <- preprocessSWAN(RGsetEx, mSet = datIlmn)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(dat0, dat)
##     checkEquivalentSummarizedExperiments(datSwan0, datSwan)
##     checkEquivalentSummarizedExperiments(datIlmn0, datIlmn)
##     checkEquivalentSummarizedExperiments(datIlmnSwan0, datIlmnSwan)
## }

## test_ratioConvert_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original examples
##     RsetEx.sub0 <- ratioConvert(MsetEx.sub, keepCN = TRUE)

##     # DelayedArray-backed version
##     MsetEx.sub <- realize(MsetEx.sub)
##     RsetEx.sub <- ratioConvert(MsetEx.sub, keepCN = TRUE)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(RsetEx.sub0, RsetEx.sub)
## }

## test_subsetByLoci_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original examples
##     loci <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611")
##     a0 <- subsetByLoci(RGsetEx.sub, includeLoci = loci)
##     b0 <- subsetByLoci(RGsetEx.sub, excludeLoci = loci)

##     # DelayedArray-backed version
##     RGsetEx.sub <- realize(RGsetEx.sub)
##     a <- subsetByLoci(RGsetEx.sub, includeLoci = loci)
##     b <- subsetByLoci(RGsetEx.sub, excludeLoci = loci)

##     # Check equivalency
##     checkEquivalentSummarizedExperiments(a0, a)
##     checkEquivalentSummarizedExperiments(b0, b)
## }

## test_utils_examples <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))

##     # Original examples
##     a0 <- head(getMethSignal(MsetEx, what = "Beta"))

##     # DelayedArray-backed version
##     MsetEx <- realize(MsetEx)
##     a <- head(getMethSignal(MsetEx, what = "Beta"))

##     # Check equivalency
##     # NOTE: getMethSignal() returns a DelayedMatrix when applied to a
##     #       DelayedArray-backed object, hence need for as.matrix()
##     checkEquals(a0, as.matrix(a))
## }
