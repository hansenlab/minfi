# a "best practices" baseline for most studies
noobPipeline <- function(rgSet, pCutoff=0.01) {
  pval <- detectionP(rgSet)
  message("Running TCGA-style noob pipeline on ", ncol(rgSet), " samples...")
  grSet <- mapToGenome(ratioConvert(preprocessNoob(rgSet)))
  message("Masking probes with detection p-value > ", pCutoff, "...")
  is.na(assays(grSet)$Beta) <- (pval[rownames(grSet),] >= pCutoff)
  message("Dropping redundant \"M\" assay matrix...")
  assays(grSet) <- assays(grSet)[c("Beta","CN")]
  message("Placing SNP probe betas in metadata(grSet)$SNPs...")
  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  message("Done.")
  return(grSet)
}
