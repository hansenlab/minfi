setMethod("bumphunter", signature(object = "GenomicRatioSet"),
          function(object, design, cluster=NULL,
                   coef = 2, cutoff = NULL, cutoffQ = 0.99,
                   maxGap = 500, smooth = FALSE,
                   smoothFunction = loessByCluster,
                   useWeights = FALSE,
                   B = 1000, verbose = TRUE,
                   type = c("M", "Beta"), ...){

              type <- match.arg(type)
              bumphunterEngine(getMethSignal(object, type),
                               design = design,
                               chr = as.factor(seqnames(object)),
                               pos = start(object),
                               cluster = cluster,
                               coef = coef,
                               cutoff = cutoff,
                               cutoffQ = cutoffQ,
                               maxGap = maxGap,
                               smooth = smooth,
                               smoothFunction = smoothFunction,
                               useWeights = useWeights,
                               B = B, verbose = verbose, ...)
          })
