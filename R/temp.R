#' @export
getopqedata <- function(path,
                        index = F,
                        xsmethod = "centWave",
                        peakwidth = c(14, 25),
                        ppm = 2.5,
                        noise = 0,
                        snthresh = 10,
                        mzdiff = -0.00395,
                        prefilter = c(3, 100),
                        mzCenterFun = "wMean",
                        integrate = 1,
                        fitgauss = FALSE,
                        verbose.columns = FALSE,
                        BPPARAM = BiocParallel::SnowParam(workers = 12),
                        rmethod = "obiwarp",
                        plottype = "none",
                        distFunc = "cor_opt",
                        profStep = 1,
                        center = 2,
                        response = 1,
                        gapInit = 0.6176,
                        gapExtend = 2.4,
                        factorDiag = 2,
                        factorGap = 1,
                        localAlignment = 0,
                        gmethod = "density",
                        bw = 0.25,
                        mzwid = 0.0021748,
                        minfrac = 1,
                        minsamp = 1,
                        gmax = 50,
                        ...) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)

        if (index) {
                cdffiles <- cdffiles[index]
        }
        xset <- xcms::xcmsSet(
                cdffiles,
                method = xsmethod,
                snthresh = snthresh,
                mzdiff = mzdiff,
                BPPARAM = BPPARAM,
                peakwidth = peakwidth,
                ppm = ppm,
                noise = noise,
                prefilter = prefilter,
                mzCenterFun = mzCenterFun,
                integrate = integrate,
                fitgauss = fitgauss,
                verbose.columns = verbose.columns,
                ...
        )
        if (index & length(index) == 1) {
                xset3 <- xset
        } else{
                xset <- xcms::group(
                        xset,
                        method = gmethod,
                        bw = bw,
                        mzwid = mzwid,
                        minfrac = minfrac,
                        minsamp = minsamp,
                        max = gmax
                )
                xset2 <- xcms::retcor(
                        xset,
                        method = rmethod,
                        plottype = plottype,
                        distFunc = distFunc,
                        profStep = profStep,
                        center = center,
                        response = response,
                        gapInit = gapInit,
                        gapExtend = gapExtend,
                        factorDiag = factorDiag,
                        factorGap = factorGap,
                        localAlignment = localAlignment
                )
                # you need group the peaks again for this corrected data
                xset2 <- xcms::group(
                        xset2,
                        method = gmethod,
                        bw = bw,
                        mzwid = mzwid,
                        minfrac = minfrac,
                        minsamp = minsamp,
                        max = gmax
                )
                xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        }
        return(xset3)
}
