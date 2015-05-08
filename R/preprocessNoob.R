preprocessNoob <- function(rgSet, offset=15, dyeCorr=TRUE, verbose = TRUE) { 
    .isRG(rgSet)
    subverbose <- max(as.integer(verbose) - 1L, 0)

    ## Extraction of the out-of-band controls
    controls <- getOOB(rgSet)
    names(controls) <- c("Cy3", "Cy5")
    mset <- preprocessRaw(rgSet)
    meth <- getMeth(mset)
    unmeth <- getUnmeth(mset)

    if (any(meth<=0)){
        meth[which(meth<=0)] <- 1
    }
    if (any(unmeth<=0)){
        unmeth[which(unmeth<=0)] <- 1
    }

    probe.type <- getProbeType(mset, withColor=TRUE)
    cy3.probes <- which(probe.type=="IGrn")
    cy5.probes <- which(probe.type=="IRed")
    d2.probes <- which(probe.type=="II")

    dat <- list(Cy3 = list(M =  as.matrix(meth[cy3.probes,]), 
                           U =  as.matrix(unmeth[cy3.probes,]),
                           D2 = as.matrix(meth[d2.probes,])), 
                Cy5 = list(M =  as.matrix(meth[cy5.probes,]), 
                           U =  as.matrix(unmeth[cy5.probes,]),
                           D2 = as.matrix(unmeth[d2.probes,])))

    rows <- lapply(dat, function(ch) {
        sapply(names(ch), function(nch) {
            nrow(ch[[nch]])
        })
    })
    last <- lapply(rows, cumsum)
    first <- lapply(names(last), function(nch) {
        last[[nch]] - rows[[nch]] + 1
    })
    names(first) <- names(last)

    estimates <- lapply(names(dat), function(nch) { 
        xf <- rbind(dat[[nch]][['M']], dat[[nch]][['U']], dat[[nch]][['D2']])
        xs <- normexp.get.xs(xf = xf, controls = controls[[nch]], offset=offset, verbose = subverbose)
        names(xs[['params']]) <- paste(names(xs[['params']]), nch, sep='.')
        names(xs[['meta']]) <- paste(names(xs[['meta']]), nch, sep='.')
        xs
    })
    names(estimates) <- names(dat)

    cy3.M <- first[['Cy3']][['M']]:last[['Cy3']][['M']]
    meth[cy3.probes, ] <- estimates[['Cy3']][['xs']][cy3.M,]

    cy3.U <- first[['Cy3']][['U']]:last[['Cy3']][['U']]
    unmeth[cy3.probes,] <- estimates[['Cy3']][['xs']][cy3.U,]

    cy5.M <- first[['Cy5']][['M']]:last[['Cy5']][['M']]
    meth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.M,]

    cy5.U <- first[['Cy5']][['U']]:last[['Cy5']][['U']]
    unmeth[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.U,]

    d2.M <- first[['Cy3']][['D2']]:last[['Cy3']][['D2']]
    d2.U <- first[['Cy5']][['D2']]:last[['Cy5']][['D2']]

    meth[d2.probes,] <- estimates[['Cy3']][['xs']][d2.M,]
    unmeth[d2.probes,] <- estimates[['Cy5']][['xs']][d2.U,]

    ## This next code block does nothing because the rgSet is not returned
    ## and pData(rgSet) is not referenced below
    for(ch in names(estimates)) { 
        chnames <- names(estimates[[ch]][['params']])
        for(nm in chnames)
            pData(rgSet)[,nm] <- estimates[[ch]][['params']][[nm]]
        varMetadata(rgSet)[chnames,] <- paste(ch, estimates[[ch]][['meta']]) 
    } 

    ## Performing dye bias normalization
    if (dyeCorr){
        ## Correction of the Illumina control probes with the background correction:
        ctrls <- getProbeInfo(rgSet, type = "Control")
        ctrls <- ctrls[ctrls$Address %in% featureNames(rgSet),]
        redControls <- getRed(rgSet)[ctrls$Address,]
        greenControls <- getGreen(rgSet)[ctrls$Address,]
        rownames(redControls) <- rownames(greenControls) <- ctrls$Type
        internal.controls <- list(Cy3 = greenControls, Cy5 = redControls)
        xcs <- lapply(names(internal.controls), function(nch) {
          xcf <- as.matrix(internal.controls[[nch]])
          normexp.get.xcs(xcf = xcf, params=estimates[[nch]][['params']])
        })
        names(xcs) <- names(dat)
        internal.controls[['Cy3']] <- xcs[["Cy3"]]
        internal.controls[['Cy5']] <- xcs[["Cy5"]]

        CG.controls <- rownames(internal.controls[[1]]) %in% c("NORM_C", "NORM_G")
        AT.controls <- rownames(internal.controls[[1]]) %in% c("NORM_A", "NORM_T")

        ## Dye bias normalizastion with the corrected Illumina control probes:
        Green.avg <- colMeans(internal.controls[["Cy3"]][CG.controls,])
        Red.avg <- colMeans(internal.controls[["Cy5"]][AT.controls,])
        R.G.ratio <- Red.avg/Green.avg

        reference <- which.min(abs(R.G.ratio-1) )
        if(verbose)
            cat('[preprocessNoob] Using sample number', reference, 'as reference level...\n')
        ref <- (Green.avg + Red.avg)[reference]/2
        if(is.na(ref))
            stop("'reference' refers to an array that is not present")
        Grn.factor <- ref/Green.avg
        Red.factor <- ref/Red.avg

        Grn <- list(M = as.matrix(meth[cy3.probes,]), 
                    U = as.matrix(unmeth[cy3.probes,]),
                    D2 = as.matrix(meth[d2.probes,]))
        Red <- list(M = as.matrix(meth[cy5.probes,]), 
                    U = as.matrix(unmeth[cy5.probes,]),
                    D2 = as.matrix(unmeth[d2.probes,]))
        Grn <- lapply(Grn, function(y) sweep(y, 2, FUN="*", Grn.factor))
        Red <- lapply(Red, function(y) sweep(y, 2, FUN="*", Red.factor))

        meth[cy3.probes,] <- Grn$M
        unmeth[cy3.probes,] <- Grn$U
        meth[d2.probes,] <- Grn$D2

        meth[cy5.probes,] <- Red$M
        unmeth[cy5.probes,] <- Red$U
        unmeth[d2.probes,] <- Red$D2
    }

    assayDataElement(mset, "Meth") <- meth
    assayDataElement(mset, "Unmeth") <- unmeth

    mset@preprocessMethod <- c( mu.norm = sprintf("Noob, dyeCorr=%s", dyeCorr))
    return(mset)
}

normexp.get.xs <- function(xf, controls, offset=50, verbose = FALSE){
    if(verbose)
        cat("[normexp.get.xs] Background mean & SD estimated from", nrow(controls), "probes\n")
    mu <- sigma <- alpha <- rep(NA, ncol(xf))
    for( i in 1:ncol(xf) ) {
        ests <- huber(controls[, i]) # from MASS
        mu[i] <- ests$mu
        sigma[i] <- ests$s
        alpha[i] <- max(huber(xf[, i])$mu - mu[i], 10)
    }
    pars <- data.frame(mu=mu, lsigma=log(sigma), lalpha=log(alpha))
    for(i in 1:ncol(xf))
        xf[,i] <- normexp.signal(as.numeric(pars[i,]), xf[,i]) # from limma
    return(list(xs=xf+offset, 
                params=data.frame(mu=mu, sigma=sigma, alpha=alpha, offset=offset),
                meta=c('background mean','background SD','signal mean','offset')))
}

normexp.get.xcs <- function(xcf, params){
    stopifnot(any(grepl("mu", names(params))))
    stopifnot(any(grepl("sigma", names(params))))
    stopifnot(any(grepl("alpha", names(params))))
    stopifnot(any(grepl("offset", names(params))))
    pars <- data.frame(mu=params[[grep("mu", names(params), value=T)]],
                       sigma = log(params[[grep("sigma", names(params), value=T)]]),
                       alpha = log(params[[grep("alpha", names(params), value=T)]]))
    for(i in 1:ncol(xcf))
        xcf[,i] <- normexp.signal(as.numeric((pars[i,])), xcf[,i] ) # from limma
    return( xcf + params[[grep('offset', names(params), value=T)]][1] )
}

