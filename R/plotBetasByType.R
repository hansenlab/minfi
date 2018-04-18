plotBetasByType <- function(data, probeTypes = NULL, legendPos = "top",
                            colors = c("black", "red", "blue"), main = "",
                            lwd = 3, cex.legend = 1) {
    # Check inputs
    if (!(is(data, "MethylSet") || is.matrix(data) || is.vector(data) ||
          is(data, "DelayedMatrix"))) {
        stop("'data' needs to be a 'MethylSet' or a matrix, 'DelayedMatrix', ",
             "or vector of beta values")
    }
    if (!is.vector(data)) {
        if (ncol(data) > 1) stop("'data' must only contain one sample")
    }
    if (!is(data, "MethylSet")) {
        r <- range(data, na.rm = TRUE)
        if (!(r[1] >= 0 && r[2] <= 1)) {
            stop("'data' should be beta values in the range [0, 1]")
        }
    }
    if (is(data, "MethylSet")) {
        if (is.null(probeTypes)) {
            typeI <- getProbeInfo(data, type = "I")[, c("Name", "nCpG")]
            typeII <- getProbeInfo(data, type = "II")[, c("Name", "nCpG")]
            probeTypes <- rbind(typeI, typeII)
            probeTypes$Type <- rep(
                x = c("I", "II"),
                times = c(nrow(typeI), nrow(typeII)))
        }
    }
    if (!all(c("Name", "Type") %in% colnames(probeTypes))) {
        stop("'probeTypes' must be a data.frame with a column 'Name' of probe ",
             "IDs and a column 'Type' indicating their design type")
    }

    # Construct 1-column matrix of beta values
    if (is(data, "MethylSet")) {
        betas <- as.matrix(getBeta(data))
    } else if (is.vector(data)) {
        betas <- matrix(
            data = data,
            ncol = 1,
            dimnames = list(names(data), NULL))
    } else {
        betas <- as.matrix(data)
    }

    # Compute densities
    betas <- betas[!is.na(betas), , drop = FALSE]
    n_probes <- nrow(betas)
    n_type1 <- sum(probeTypes$Type == "I" &
                       probeTypes$Name %in% rownames(betas))
    n_type2 <- sum(probeTypes$Type == "II" &
                       probeTypes$Name %in% rownames(betas))
    betas_density <- density(betas)
    typeI_density <- suppressWarnings(
        density(
            x = betas[rownames(betas) %in%
                          probeTypes$Name[probeTypes$Type == "I"]],
            weights = rep(1 / n_probes, n_type1)))
    typeII_density <- suppressWarnings(
        density(
            x = betas[rownames(betas) %in%
                          probeTypes$Name[probeTypes$Type == "II"]],
            weights = rep(1 / n_probes, n_type2)))

    # Plot densities
    plot(
        x = betas_density,
        main = main,
        xlab = "Beta values",
        ylim = c(0, max(betas_density$y)),
        lwd = lwd,
        col = colors[1])
    lines(
        x = typeI_density,
        col = colors[2],
        lty = 2,
        lwd = lwd)
    lines(
        x = typeII_density,
        col = colors[3],
        lty = 2,
        lwd = lwd)
    legend(
        x = legendPos,
        y = c("All probes", "Infinium I", "Infinium II"),
        col = colors,
        lwd = lwd,
        bg = "white",
        cex = cex.legend,
        lty = c(1, 2, 2))
}
