mdsPlot <- function(dat, numPositions=1000, sampNames=NULL, sampGroups=NULL, xlim, ylim, pch=1, pal=brewer.pal(8, "Dark2"), legendPos="bottomleft", legendNCol, main=NULL) {
    if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
        b <- getBeta(dat)
        if(is.null(main)) main <- sprintf("Beta MDS\n%d most variable positions", numPositions)
    } else if (is(dat, "matrix")) {
        b <- dat
        if(is.null(main)) main <- sprintf("MDS\n%d most variable positions", numPositions)
    } else {
        stop("dat must be an 'MethylSet' or 'matrix'.")
    }
    o <- order(-rowVars(b))[1:numPositions]
    d <- dist(t(b[o,]))
    fit <- cmdscale(d)
    if (missing(xlim)) xlim <- range(fit[,1]) * 1.2
    if (missing(ylim)) ylim <- range(fit[,2]) * 1.2
    if (is.null(sampGroups)) sampGroups <- rep(1, numPositions)
    sampGroups <- as.factor(sampGroups)
    col <- pal[sampGroups]
    if (is.null(sampNames)) {
        plot(fit[,1], fit[,2], col=col, pch=pch, xlim=xlim, ylim=ylim, xlab="", ylab="", main=main)
    } else {
        plot(0,0, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", main=main)
        text(fit[,1], fit[,2], sampNames, col=col)
    }
    numGroups <- length(levels(sampGroups))
    if (missing(legendNCol)) legendNCol <- numGroups
    if (numGroups>1) {
        legend(legendPos, legend=levels(sampGroups), ncol=legendNCol, text.col=pal[1:numGroups])
    }
}

plotCpg <- function(dat, cpg, pheno, type=c("categorical", "continuous"), measure=c("beta", "M"), ylim=NULL, ylab = NULL, xlab="", fitLine=TRUE, mainPrefix=NULL, mainSuffix=NULL) {
    if (is.numeric(cpg)) cpg <- rownames(dat)[cpg]
    type <- match.arg(type)
    measure <- match.arg(measure)
    if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
        if (measure=="beta") {
            x <- getBeta(dat[cpg,])
            if(is.null(ylab)) ylab <- "Beta"
            if (is.null(ylim)) ylim <- c(0,1)
        } else if (measure=="M") {
            x <- getM(dat[cpg,])
            if(is.null(ylab)) ylab <- "M"
            if (is.null(ylim)) ylim <- range(x)
        }
    } else {
        if(is.null(ylab)) ylab <- "unknown"
        x <- dat
        if (is.vector(x)) 
            x <- matrix(x, ncol = 1)
    }
    main <- paste(mainPrefix, cpg, mainSuffix)
    names(main) <- cpg
    for (i in cpg) {
        if (type=="categorical") {
            stripchart(x[i,] ~ pheno, vertical=TRUE, method="jitter", jitter=0.15, ylab=ylab, ylim=ylim, main=main[i])
        } else if (type=="continuous") {
            plot(pheno, x[i,], ylab=ylab, xlab=xlab, ylim=ylim, main=main[i])
            if (fitLine) abline(lm(x[i,]~pheno), col="blue")
        }
    }
}



