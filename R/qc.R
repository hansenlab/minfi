qcReport <- function(rgSet, sampNames=NULL, sampGroups=NULL, pdf="qcReport.pdf", maxSamplesPerPage=24, 
        controls=c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "EXTENSION", "HYBRIDIZATION", 
        "NON-POLYMORPHIC", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL")) {
    .isRG(rgSet)
    if (is.null(sampNames)) sampNames <- sampleNames(rgSet)
    n <- ncol(rgSet)
    o <- rev(order(sampNames))
    rgSet <- rgSet[,o]
    sampNames <- sampNames[o]
    sampGroups <- sampGroups[o]
    sampGroups <- as.factor(sampGroups)
    if (is.null(sampGroups)) sampGroups <- rep(1, n)
    numPages <- ceiling(n/maxSamplesPerPage)
    samplesPerPage <- ceiling(n/numPages)
    sampleIdxs <- suppressWarnings(split(1:n, rep(1:numPages, each=samplesPerPage)))
    pdf(file=pdf, width=8, height=11)
    par(mfrow=c(2,1))
    densityPlot(rgSet, sampGroups=sampGroups, main="Beta", xlab="Beta")
    plot.new()
    par(mfrow=c(1,1), oma=c(2,10,1,1))
    for (sampleIdx in sampleIdxs) {
        densityBeanPlot(rgSet[,sampleIdx], sampGroups=sampGroups[sampleIdx], sampNames=sampNames[sampleIdx])
    }    
    for (controlType in controls) {
        for (sampleIdx in sampleIdxs) {
            controlStripPlot(rgSet[,sampleIdx], sampNames=sampNames[sampleIdx], controls=controlType)
        }
    }
    dev.off()
}


controlStripPlot <- function(rgSet, controls=c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II"), 
    sampNames=NULL, xlim=c(5, 17)) {
    .isRG(rgSet)
    
    r <- getRed(rgSet)
    g <- getGreen(rgSet)
    
    for (controlType in controls) {
        ctrlAddress <- getControlAddress(rgSet, controlType = controlType)
        
        ## Red channel
        ctlWide <- log2(r[ctrlAddress,])
        if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
        ctlR <- melt(ctlWide, varnames=c("address", "sample"))
        
        ## Green channel
        ctlWide <- log2(g[ctrlAddress,])
        if (!is.null(sampNames)) colnames(ctlWide) <- sampNames
        ctlG<- melt(ctlWide, varnames=c("address", "sample"))
        
        ## Plot
        ctl <- rbind(cbind(channel="Red", ctlR), cbind(channel="Green", ctlG))
        if (any((ctl$value<xlim[1]) | (ctl$value>xlim[2]))) message("Warning: ", controlType, " probes outside plot range")
    fig <- xyplot(sample ~ value | channel, groups=channel, horizontal=TRUE, pch=19, 
                  col=c("darkred", "darkgreen"),
                  xlab="Log2 Intensity", xlim=xlim,
                  main=paste("Control:", controlType), layout=c(2,1), data=ctl,
                  panel = function(x,y,...) {
                      panel.stripplot(x,y,...)
                      panel.abline(h=(as.numeric(y)-0.5), lty=3, col="grey70")
                  })
        print(fig) 
    } 
}


densityBeanPlot <- function(dat, sampGroups=NULL, sampNames=NULL, main=NULL, pal=brewer.pal(8, "Dark2"), numPositions=10000) {
    if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
        b <- getBeta(dat)
    } else if (is(dat, "matrix")) {
        b <- dat
    } else {
        stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet'  or matrix.")
    }
    n <- ncol(b)
    if (!is.null(sampNames)) colnames(b) <- sampNames
    if (is.null(main)) main <- "Beta"
    if (is.null(sampGroups))  sampGroups <- rep(1, n)
    sampGroups <- as.factor(sampGroups)
    col <- lapply(sampGroups, function(x) rep(pal[x],4))
    if(is.null(numPositions))
        idx <- 1:dim(dat)[1]
    else
        idx <- sample(nrow(b), numPositions)
    x <- melt(b[idx, ], varnames=c("cpg", "sample"))    
    o <- order(colnames(b))
    beanplot(value ~ sample, horizontal=TRUE, what=c(0,1,1,0), log="", las=1, ylim=c(0,1), 
        xlab="Beta", main=main, col=col[o], data=x, cex.lab=0.9, beanlinewd=1, border=NA)
    abline(h=1:(n+1)-0.5, lty=3, col="grey70")
}


densityPlot <- function(dat, sampGroups=NULL, main="", xlab="Beta", pal=brewer.pal(8, "Dark2"), xlim, ylim, add=TRUE, legend=TRUE, ...) {
    if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
        b <- getBeta(dat)
    } else if (is(dat, "matrix")) {
        b <- dat
    } else {
        stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet'  or matrix.")
    }
    d <- apply(b, 2, density, na.rm=TRUE)
    if (missing(ylim)) ylim <- range(sapply(d, function(i) range(i$y)))
    if (missing(xlim)) xlim <- range(sapply(d, function(i) range(i$x)))
    if (is.null(sampGroups)) {
        sampGroups <- rep(1, ncol(b))
    } else if (length(sampGroups)==1) {
        sampGroups <- rep(sampGroups, ncol(b))
    }
    sampGroups <- as.factor(sampGroups)
    if (add) {
        plot(0, type="n", ylim=ylim, xlim=xlim, ylab="Density", xlab=xlab, main=main, ...)
        abline(h=0, col="grey80")
    }
    for (i in 1:length(d)) {
        lines(d[[i]], col=pal[sampGroups[i]])
    }
    if (legend & length(levels(sampGroups))>1) legend("topright", legend=levels(sampGroups), text.col=pal)
}
