# a "best practices" baseline for most studies
tcgaPipeline <- function(rgSet, pCutoff=0.05) {
  pval <- detectionP(rgSet)
  message("Running TCGA-style (noob) pipeline on ", ncol(rgSet), " samples...")
  grSet <- ratioConvert(mapToGenome(preprocessNoob(rgSet)))
  message("Masking probes with detection p-value > ", pCutoff, "...")
  is.na(assays(grSet)$Beta) <- (pval[rownames(grSet),] >= pCutoff)
  message("Placing SNP probe betas in metadata(grSet)$SNPs...")
  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  message("Done.")
  return(grSet)
}
