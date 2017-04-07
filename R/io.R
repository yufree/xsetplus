#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the references
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @export
getdata <-
        function(path,
                 index = F,
                 BPPARAM = BiocParallel::SnowParam(workers = 12),
                 pmethod = 'hplcorbitrap',
                 ...) {
                cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
                if (index) {
                        cdffiles <- cdffiles[index]
                }
                if (pmethod == 'hplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(10, 60),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcorbitrap') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 2.5,
                                        peakwidth = c(5, 20),
                                        prefilter = c(3, 5000),
                                        ...
                                )
                        xset <-
                                xcms::group(xset, bw = 2, mzwid = 0.015)
                        xset2 <- xcms::retcor(xset)
                        # you need group the peaks again for this corrected data
                        xset2 <-
                                xcms::group(xset2, bw = 2, mzwid = 0.015)
                        xset3 <-
                                xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                } else if (pmethod == 'hplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'hplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(10, 60),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 5,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplcqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 30,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.025)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else if (pmethod == 'uplchqtof') {
                        xset <-
                                xcms::xcmsSet(
                                        cdffiles,
                                        BPPARAM = BPPARAM,
                                        method = "centWave",
                                        ppm = 15,
                                        peakwidth = c(5, 20),
                                        prefilter = c(0, 0),
                                        ...
                                )
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <-
                                        xcms::group(xset2,
                                                    bw = 2,
                                                    mzwid = 0.015)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                } else{
                        xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM, ...)
                        if (index & length(index) == 1) {
                                xset3 <- xset
                        } else{
                                xset <- xcms::group(xset)
                                xset2 <- xcms::retcor(xset)
                                # you need group the peaks again for this corrected data
                                xset2 <- xcms::group(xset2)
                                xset3 <-
                                        xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
                        }
                }
                return(xset3)
        }
#' Get the csv files to be submitted to Metaboanalyst
#' @param xset the xcmsset object which you want to submitted to Metaboanalyst
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param name file name
#' @return dataframe with data needed for Metaboanalyst if your want to perform local analysis.
#' @export
getupload <-
        function(xset,
                 method = "medret",
                 intensity = 'inio',
                 name = 'Peaklist') {
                peakIntensities <- xcms::groupval(xset, method, intensity)
                if (intensity == "intb") {
                        peakIntensities[is.na(peakIntensities)] = 0
                }
                data <-
                        rbind(group = as.character(xcms::phenoData(xset)$class), peakIntensities)
                filename <- paste0(name, '.csv')
                utils::write.csv(data, file = filename)
                return(data)
        }
#' Get the report for technique replicates.
#' @param xset the xcmsset object which for all of your technique replicates for one sample
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data
#' @export
gettechrep <- function(xset,
                       method =  'medret',
                       intensity = 'into') {
        data <- t(xcms::groupval(xset, method, intensity))
        lv <- xset@phenoData[, 1]
        mean <- stats::aggregate(data, list(lv), mean)
        sd <- stats::aggregate(data, list(lv), sd)
        suppressWarnings(rsd <- sd / mean * 100)
        result <-
                data.frame(cbind(t(mean[, -1]), t(sd[, -1]), t(rsd[, -1])))
        name <- unique(lv)
        colnames(result) <-
                c(paste0(name, 'mean'),
                  paste0(name, 'sd'),
                  paste0(name, 'rsd%'))
        datap <- xcms::groups(xset)
        report <- cbind.data.frame(datap, result)
        return(report)
}
#' Get the report for samples with technique replicates
#' @param xset the xcmsset object all of samples with technique replicates
#' @param anno logical if set as True, it will return the table for further annotation, default false
#' @param peaklist logical if set as True, it will return csv files for metaboanalyst, default false
#' @param file file name for the peaklist
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @return dataframe with mean, standard deviation and RSD for those technique replicates combined with raw data for all of the samples if anno and peaklist are defaults false.
#' @export
gettechbiorep <-
        function(xset,
                 anno = F,
                 peaklist = F,
                 file = NULL,
                 method =  'medret',
                 intensity = 'into') {
                data <- t(xcms::groupval(xset, method, intensity))
                lv <- xset@phenoData[, 1]
                lv2 <- xset@phenoData[, 2]
                mean <- stats::aggregate(data, list(lv, lv2), mean)
                sd <- stats::aggregate(data, list(lv, lv2), sd)
                suppressWarnings(rsd <- sd / mean * 100)
                result <-
                        data.frame(cbind(t(mean[, -c(1:2)]), t(sd[, -c(1:2)]), t(rsd[, -c(1:2)])))
                name <- unique(c(paste0(lv, lv2)))
                colnames(result) <-
                        c(paste0(name, 'mean'),
                          paste0(name, 'sd'),
                          paste0(name, 'rsd%'))
                datap <- xcms::groups(xset)
                report <- cbind.data.frame(datap, result)
                if (anno) {
                        anno <-
                                as.data.frame(cbind(xset@groups[, 1], xset@groups[, 4], cbind(t(
                                        mean[, -c(1:2)]
                                ))))
                        colnames(anno) <- c('mz', 'time', name)
                        return(anno)
                } else if (peaklist) {
                        result <- data.frame(t(mean[, -c(1:2)]))
                        data <- rbind(group = name, result)
                        utils::write.csv(data, file = file)
                } else{
                        return(report)
                }
        }
#' get the data of QC compound for a group of data
#' @param path data path for your QC samples
#' @param mzrange mass of the QC compound
#' @param rtrange retention time of the QC compound
#' @param index index of the files contained QC compounds, default is all of the compounds
#' @return number vector, each number indicate the peak area of that mass and retention time range
#' @export
getQCraw <- function(path, mzrange, rtrange, index = NULL) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        if (index) {
                cdffiles <- cdffiles[index]
        }
        nsamples <- length(cdffiles)
        area <- numeric()
        for (i in 1:nsamples) {
                RAW <- xcms::xcmsRaw(cdffiles[i])
                peak <- xcms::rawEIC(RAW, mzrange, rtrange)
                area[i] <- sum(peak$intensity)
        }
        return(area)
}
