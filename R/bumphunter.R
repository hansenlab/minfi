setMethod("bumphunter", signature(object = "GenomicRatioSet"),
          function(object, design, cluster=NULL,
                   coef = 2,
                   cutoff=NULL, pickCutoff=FALSE, pickCutoffQ=0.99, 
                   maxGap = 500, smooth = FALSE,
                   smoothFunction = locfitByCluster,
                   useWeights = FALSE,
                   B=ncol(permutations), permutations=NULL,
                   verbose = TRUE,
                   type = c("Beta","M"), ...){
              
              type <- match.arg(type)
              bumphunterEngine(getMethSignal(object, type),
                               design = design,
                               chr = as.factor(seqnames(object)),
                               pos = start(object),
                               cluster = cluster,
                               coef = coef,
                               cutoff=cutoff,
                               pickCutoff=pickCutoff,
                               pickCutoffQ=pickCutoffQ, 
                               maxGap = maxGap,
                               smooth = smooth,
                               smoothFunction = smoothFunction,
                               useWeights = useWeights,
                               B=B,
                               permutations=permutations,
                               verbose = verbose, ...)
          })
