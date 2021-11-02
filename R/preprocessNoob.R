# Global variables -------------------------------------------------------------

utils::globalVariables("oob")

# Internal functions -----------------------------------------------------------

# TODO: Profile
normexp.get.xs <- function(xf, controls, offset = 50, verbose = FALSE) {
    if (verbose) {
        message(
            "[normexp.get.xs] Background mean & SD estimated from ",
            nrow(controls),
            " probes")
    }
    mu <- sigma <- alpha <- rep(NA_real_, ncol(xf))
    for (i in seq_len(ncol(xf))) {
        ests <- huber(controls[, i])
        mu[i] <- ests$mu
        sigma[i] <- ests$s
        alpha[i] <- max(huber(xf[, i])$mu - mu[i], 10)
    }
    pars <- data.frame(mu = mu, lsigma = log(sigma), lalpha = log(alpha))
    for (i in seq_len(ncol(xf))) {
        xf[, i] <- normexp.signal(as.numeric(pars[i, ]), xf[, i])
    }
    list(
        xs = xf + offset,
        params = data.frame(
            mu = mu, sigma = sigma, alpha = alpha, offset = offset),
        meta = c("background mean", "background SD", "signal mean", "offset"))
}

# TODO: Profile
normexp.get.xcs <- function(xcf, params) {
    stopifnot(any(grepl("mu", names(params))),
              any(grepl("sigma", names(params))),
              any(grepl("alpha", names(params))),
              any(grepl("offset", names(params))))
    pars <- data.frame(
        mu = params[[grep("mu", names(params), value = TRUE)]],
        sigma = log(params[[grep("sigma", names(params), value = TRUE)]]),
        alpha = log(params[[grep("alpha", names(params), value = TRUE)]]))
    for (i in seq_len(ncol(xcf))) {
        xcf[, i] <- normexp.signal(as.numeric(pars[i, ]), xcf[, i])
    }

    xcf + params[[grep("offset", names(params), value = TRUE)]][1]
}

# From TJT (2016-06-16)
# Performing dye bias normalization
#
# "single" = just reciprocate out the dye bias, don't use a reference. (similar
#            to, but implemented differently from, unmaintained "asmn" package
#            by Decker et al., doi:10.4161/epi.26037)
# "reference" = use the least-worst sample in the batch (previous default)
#
# "single" is now the default: it provides single-sample preprocessing and
# betas/M-values produced by this method are identical to those from the
# "reference" version used in (e.g.) the TCGA data processing pipeline.
dyeCorrection <- function(Meth, Unmeth, Red, Green, control_probes,
                          Green_probes, Red_probes, d2.probes, estimates,
                          array_type, dyeMethod, verbose) {

    # Background correct the Illumina normalization controls
    redControls <- Red[control_probes$Address, , drop = FALSE]
    greenControls <- Green[control_probes$Address, ,drop = FALSE]
    rownames(redControls) <- rownames(greenControls) <-
        control_probes$Type
    internal.controls <- list(Green = greenControls, Red = redControls)
    xcs <- lapply(names(internal.controls), function(nch) {
        xcf <- as.matrix(internal.controls[[nch]])
        normexp.get.xcs(xcf = xcf, params = estimates[[nch]][["params"]])
    })
    names(xcs) <- names(internal.controls)
    internal.controls[["Green"]] <- xcs[["Green"]]
    internal.controls[["Red"]] <- xcs[["Red"]]

    if (array_type %in% c("IlluminaHumanMethylation27k")) {
        CG.controls <- which(
            rownames(internal.controls[[1]]) %in% c("Normalization-Green"))
        AT.controls <- which(
            rownames(internal.controls[[1]]) %in% c("Normalization-Red"))
    } else {
        CG.controls <- which(
            rownames(internal.controls[[1]]) %in% c("NORM_C", "NORM_G"))
        AT.controls <- which(
            rownames(internal.controls[[1]]) %in% c("NORM_A", "NORM_T"))
    }

    # Dye bias normalization with the corrected Illumina control probes
    Green.avg <- colMeans2(x = internal.controls[["Green"]], rows = CG.controls)
    Red.avg <- colMeans2(x = internal.controls[["Red"]], rows = AT.controls)
    R.G.ratio <- Red.avg / Green.avg

    if (dyeMethod == "single") {
        if (verbose) {
            message(
                "[dyeCorrection] Applying R/G ratio flip to fix dye bias")
        }
        Red.factor <- 1 / R.G.ratio
        Green.factor <- 1
    } else if (dyeMethod == "reference") {
        reference <- which.min(abs(R.G.ratio - 1))
        if (verbose) {
            message(
                "[dyeCorrection] Using sample number ",
                reference,
                " as reference level")
        }
        ref <- (Green.avg + Red.avg)[reference] / 2
        if (is.na(ref)) {
            stop("'reference' refers to an array that is not present")
        }
        Green.factor <- ref / Green.avg
        Red.factor <- ref / Red.avg
    }

    Green <- list(
        M = Meth[Green_probes,, drop=FALSE],
        U = Unmeth[Green_probes,, drop=FALSE],
        D2 = Meth[d2.probes,, drop=FALSE])
    Red <- list(
        M = Meth[Red_probes,, drop=FALSE],
        U = Unmeth[Red_probes,, drop=FALSE],
        D2 = Unmeth[d2.probes,, drop=FALSE ])

    # NOTE: Adjust Red regardless of reference or equalization approach
    #       but only adjust Green if using older reference method
    Red <- lapply(
        X = Red,
        FUN = function(x) {
            sweep(x, MARGIN = 2, STATS = Red.factor, FUN = "*")
        })
    Meth[Red_probes, ] <- Red$M
    Unmeth[Red_probes, ] <- Red$U
    Unmeth[d2.probes, ] <- Red$D2
    if (dyeMethod == "reference") {
        Green <- lapply(
            X = Green,
            FUN = function(x) {
                sweep(x, MARGIN = 2, STATS = Green.factor, FUN = "*")
            })
        Meth[Green_probes, ] <- Green$M
        Unmeth[Green_probes, ] <- Green$U
        Meth[d2.probes, ] <- Green$D2
    }

    list(Meth = Meth, Unmeth = Unmeth)
}

# Internal generics ------------------------------------------------------------

# `...` are additional arguments passed to methods.
setGeneric(
    ".preprocessNoob",
    function(Meth, Unmeth, GreenOOB, RedOOB, Green_probes, Red_probes,
             d2.probes, offset, dyeCorr, Red, Green, control_probes, array_type,
             dyeMethod, verbose, ...)
        standardGeneric(".preprocessNoob"),
    signature = c("Meth", "Unmeth")
)

# Internal methods -------------------------------------------------------------

setMethod(
    ".preprocessNoob",
    c("matrix", "matrix"),
    function(Meth, Unmeth, GreenOOB, RedOOB, Green_probes, Red_probes,
             d2.probes, offset, dyeCorr, Red, Green, control_probes, array_type,
             dyeMethod, verbose) {
        subverbose <- max(as.integer(verbose) - 1L, 0)

        # Threshold Meth and Unmeth to be positive
        Meth[Meth <= 0] <- 1L
        Unmeth[Unmeth <= 0] <- 1L
        
        # NormExp estimates for Green and Red
        dat <- list(
            Green = list(
                M =  Meth[Green_probes, , drop = FALSE],
                U =  Unmeth[Green_probes, , drop = FALSE],
                D2 = Meth[d2.probes, , drop = FALSE]),
            Red = list(
                M =  Meth[Red_probes, , drop = FALSE],
                U =  Unmeth[Red_probes, , drop = FALSE],
                D2 = Unmeth[d2.probes, , drop = FALSE]))
        estimates <- lapply(names(dat), function(nch, oob) {
            xf <- rbind(
                dat[[nch]][["M"]],
                dat[[nch]][["U"]],
                dat[[nch]][["D2"]])
            xs <- normexp.get.xs(
                xf = xf,
                controls = oob[[nch]],
                offset = offset,
                verbose = subverbose)
            names(xs[["params"]]) <- paste(
                names(xs[["params"]]), nch, sep = ".")
            names(xs[["meta"]]) <- paste(names(xs[["meta"]]), nch, sep = ".")
            xs
        }, oob = list(Green = GreenOOB, Red = RedOOB))
        names(estimates) <- names(dat)

        # Correct for Green and Red in Meth and Unmeth
        rows <- lapply(dat, function(x) vapply(x, nrow, integer(1L)))
        last <- lapply(rows, cumsum)
        first <- Map(function(last, rows) last - rows + 1, last, rows)
        if (length(Green_probes) > 0) {
            Green.M <- seq(first[["Green"]][["M"]], last[["Green"]][["M"]])
            Meth[Green_probes, ] <- estimates[["Green"]][["xs"]][Green.M, ]
            Green.U <- seq(first[["Green"]][["U"]], last[["Green"]][["U"]])
            Unmeth[Green_probes, ] <- estimates[["Green"]][["xs"]][Green.U, ]
        }
        if (length(Red_probes) > 0) {
            Red.M <- seq(first[["Red"]][["M"]], last[["Red"]][["M"]])
            Meth[Red_probes, ] <- estimates[["Red"]][["xs"]][Red.M, ]
            Red.U <- seq(first[["Red"]][["U"]], last[["Red"]][["U"]])
            Unmeth[Red_probes, ] <- estimates[["Red"]][["xs"]][Red.U, ]
        }
        if (length(d2.probes) > 0) {
            d2.M <- seq(first[["Green"]][["D2"]], last[["Green"]][["D2"]])
            d2.U <- seq(first[["Red"]][["D2"]], last[["Red"]][["D2"]])
            Meth[d2.probes, ] <- estimates[["Green"]][["xs"]][d2.M, ]
            Unmeth[d2.probes, ] <- estimates[["Red"]][["xs"]][d2.U, ]
        }

        if (!dyeCorr) {
            return(list(Meth = Meth, Unmeth = Unmeth))
        }
        dyeCorrection(
            Meth = Meth,
            Unmeth = Unmeth,
            Red = Red,
            Green = Green,
            control_probes = control_probes,
            Green_probes = Green_probes,
            Red_probes = Red_probes,
            d2.probes = d2.probes,
            estimates = estimates,
            array_type = array_type,
            dyeMethod = dyeMethod,
            verbose = verbose)
    }
)

setMethod(
    ".preprocessNoob",
    c("DelayedMatrix", "DelayedMatrix"),
    function(Meth, Unmeth, GreenOOB, RedOOB, Green_probes, Red_probes,
             d2.probes, offset, dyeCorr, Red, Green, control_probes,
             array_type, dyeMethod, verbose, BPREDO = list(),
             BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink objects of appropriate dimensions
        # and type
        # NOTE: These are ultimately coerced to the output DelayedMatrix
        #       objects, `Meth` and `Unmeth`
        # NOTE: Noob can return non-integer values, so write to a numeric sink
        # NOTE: Don't do `Unmeth_sink <- Meth__sink` or else these will
        #       reference the same object and clobber each other when written
        #       to!
        ans_type <- "double"
        M_sink <- DelayedArray::AutoRealizationSink(
            dim = dim(Meth),
            dimnames = dimnames(Meth),
            type = ans_type)
        on.exit(close(M_sink))
        U_sink <- DelayedArray::AutoRealizationSink(
            dim = dim(Unmeth),
            dimnames = dimnames(Unmeth),
            type = ans_type)
        on.exit(close(U_sink), add = TRUE)

        # Set up ArrayGrid instances over `Meth`, `Unmeth`, `M_sink`, and
        # `U_sink`. Also set up "parallel" ArrayGrid instances over `GreenOOB`
        # and `RedOOB`, as well as `Red` and `Green` if performing dye
        # correction.
        Meth_grid <- colAutoGrid(Meth)
        Unmeth_grid <- colAutoGrid(Unmeth)
        M_sink_grid <- RegularArrayGrid(
            refdim = dim(M_sink),
            spacings = c(nrow(M_sink), ncol(M_sink) / length(Meth_grid)))
        U_sink_grid <- M_sink_grid
        GreenOOB_grid <- RegularArrayGrid(
            refdim = dim(GreenOOB),
            spacings = c(nrow(GreenOOB), ncol(GreenOOB) / length(Meth_grid)))
        RedOOB_grid <- RegularArrayGrid(
            refdim = dim(RedOOB),
            spacings = c(nrow(RedOOB), ncol(RedOOB) / length(Meth_grid)))
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(identical(dim(Meth_grid), dim(Unmeth_grid)),
                  identical(dim(Meth_grid), dim(M_sink_grid)),
                  identical(dim(Meth_grid), dim(GreenOOB_grid)),
                  identical(dim(Meth_grid), dim(RedOOB_grid)))

        # Loop over blocks of `Meth`, `Unmeth`, `GreenOOB`, and `RedOOB`, as
        # well as `Red` and `Green` if doing dye correction, and write to
        # `M_sink` and `U_sink`.
        if (dyeCorr) {
            Red_grid <- RegularArrayGrid(
                refdim = dim(Red),
                spacings = c(nrow(Red), ncol(Red) / length(Meth_grid)))
            Green_grid <- RegularArrayGrid(
                refdim = dim(Green),
                spacings = c(nrow(Green), ncol(Green) / length(Meth_grid)))
            # Sanity check ArrayGrid objects have the same dim
            stopifnot(dim(Red_grid) == dim(Green_grid),
                      dim(Red_grid) == dim(Meth_grid))
            blockMapplyWithRealization(
                FUN = .preprocessNoob,
                Meth = Meth,
                Unmeth = Unmeth,
                GreenOOB = GreenOOB,
                RedOOB = RedOOB,
                Red = Red,
                Green = Green,
                MoreArgs = list(
                    Green_probes = Green_probes,
                    Red_probes = Red_probes,
                    d2.probes = d2.probes,
                    offset = offset,
                    dyeCorr = dyeCorr,
                    control_probes = control_probes,
                    array_type = array_type,
                    dyeMethod = dyeMethod,
                    verbose = verbose),
                sinks = list(M_sink, U_sink),
                dots_grids = list(Meth_grid, Unmeth_grid,
                                  GreenOOB_grid, RedOOB_grid,
                                  Red_grid, Green_grid),
                sinks_grids = list(M_sink_grid, U_sink_grid),
                BPREDO = BPREDO,
                BPPARAM = BPPARAM)
        } else {
            blockMapplyWithRealization(
                FUN = .preprocessNoob,
                Meth = Meth,
                Unmeth = Unmeth,
                GreenOOB = GreenOOB,
                RedOOB = RedOOB,
                MoreArgs = list(
                    oob = oob,
                    Green_probes = Green_probes,
                    Red_probes = Red_probes,
                    d2.probes = d2.probes,
                    offset = offset,
                    dyeCorr = dyeCorr,
                    control_probes = control_probes,
                    array_type = array_type,
                    dyeMethod = dyeMethod,
                    verbose = verbose),
                sinks = list(M_sink, U_sink),
                dots_grids = list(Meth_grid, Unmeth_grid,
                                  GreenOOB_grid, RedOOB_grid),
                sinks_grids = list(M_sink_grid, U_sink_grid),
                BPREDO = BPREDO,
                BPPARAM = BPPARAM)
        }

        # Return as DelayedMatrix objects
        M <- as(M_sink, "DelayedArray")
        U <- as(U_sink, "DelayedArray")

        list(Meth = M, Unmeth = U)
    }
)

# Exported functions -----------------------------------------------------------

# TODO: Document: because we simultaneously walk over column-blocks of `Meth`,
#       and `Unmeth`, as well as `Red` and `Green` if doing dye correction, the
#       number of elements loaded into memory is doubled or quadrupled. If
#       running into memory issues, try halving/quartering
#       getOption("DelayedArray.block.size")
preprocessNoob <- function(rgSet, offset = 15, dyeCorr = TRUE, verbose = FALSE,
                           dyeMethod = c("single", "reference")) {

    # Check inputs
    .isRGOrStop(rgSet)
    dyeMethod <- match.arg(dyeMethod)

    # Extract data to pass to low-level functions that construct `M` and `U`
    oob <- getOOB(rgSet)
    GreenOOB <- oob[["Grn"]]
    RedOOB <- oob[["Red"]]
    MSet <- preprocessRaw(rgSet)
    probe.type <- getProbeType(MSet, withColor = TRUE)
    Green_probes <- which(probe.type == "IGrn")
    Red_probes <- which(probe.type == "IRed")
    d2.probes <- which(probe.type == "II")
    Meth <- getMeth(MSet)
    Unmeth <- getUnmeth(MSet)

    if (dyeCorr) {
        control_probes <- getProbeInfo(rgSet, type = "Control")
        control_probes <- control_probes[
            control_probes$Address %in% rownames(rgSet), ]
        Red <- getRed(rgSet)
        Green <- getGreen(rgSet)
        array_type <- rgSet@annotation[["array"]]
    } else {
        control_probes <- NULL
        Red <- NULL
        Green <- NULL
        array_type <- NULL
    }
    M_and_U <- .preprocessNoob(
        Meth = Meth,
        Unmeth = Unmeth,
        GreenOOB = GreenOOB,
        RedOOB = RedOOB,
        Green_probes = Green_probes,
        Red_probes = Red_probes,
        d2.probes = d2.probes,
        offset = offset,
        dyeCorr = dyeCorr,
        Red = Red,
        Green = Green,
        control_probes = control_probes,
        array_type = array_type,
        dyeMethod = dyeMethod,
        verbose = verbose)

    # Construct MethylSet
    assay(MSet, "Meth") <- M_and_U[["Meth"]]
    assay(MSet, "Unmeth") <- M_and_U[["Unmeth"]]
    # TODO: Is this the correct way to add a preprocessing method? It is
    #       shown as NA by show,MethylSet-method
    MSet@preprocessMethod <- c(
        mu.norm = sprintf("Noob, dyeCorr=%s, dyeMethod=%s", dyeCorr, dyeMethod))
    MSet
}

# TODO: Switch from `Meth` to `M` and `Unmeth` to `U`
# TODO: Rationalise argument names and order across functions, set defaults
#       assuming dyeCorr is FALSE
