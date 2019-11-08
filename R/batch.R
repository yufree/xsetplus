#' LIMMA to test batch effects by p-value and q-value
#' @param data data as mzrt profile
#' @param lv vector for the group information
#' @param batch vector for the batch information
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no batch index, corresponding part would miss.
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{famtcor}}
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- limmafit(list$data,list$group$class)
#' }
#' @export
limmafit <- function(data, lv, batch = NULL){
        mod <- stats::model.matrix(~lv)
        mod0 <- as.matrix(c(rep(1, ncol(data))))
        datacor <- signal <- error <- pValues <-  qValues <- NULL
        if(is.null(batch)){
                batch <- NULL
                # limma fit
                lmfit <- limma::lmFit(data, mod)
                signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(mod[, 1:nlevels(lv)])
                error <- data - signal
                rownames(signal) <- rownames(error) <- rownames(data)
                colnames(signal) <- colnames(error) <- colnames(data)
                # find the peaks with significant differences by F test
                # with BH correction for fdr control without correction
                pValues = sva::f.pvalue(data, mod, mod0)
                qValues = stats::p.adjust(pValues, method = "BH")
        }else{
                modcor <- cbind(mod,batch)
                modcor0 <- cbind(mod0,batch)
                lmfit <- limma::lmFit(data, modcor)
                # data decomposition with batch
                batch <- lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + NCOL(batch))] %*% t(modcor[, (nlevels(lv) + 1):(nlevels(lv) + NCOL(batch))])
                signal <- lmfit$coef[, 1:nlevels(lv)] %*% t(modcor[, 1:nlevels(lv)])
                error <- data - signal - batch
                datacor <- signal + error
                rownames(datacor) <- rownames(batch) <- rownames(signal) <- rownames(error) <- rownames(data)
                colnames(datacor) <- colnames(batch) <- colnames(signal) <- colnames(error) <- colnames(data)

                # find the peaks with significant differences by F test
                # with BH correction for fdr control
                pValues = sva::f.pvalue(data, modcor, modcor0)
                qValues = stats::p.adjust(pValues, method = "BH")
        }
        # get the results as list
        li <- list(data, datacor, signal, batch, error, pValues, qValues)
        names(li) <- c("data","dataCorrected","signal","batch", "error", "p-values", "q-values")
        return(li)
}

#' Use Surrogate Variable Analysis(SVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Surrogate Variable Analysis(SVA) to correct the unknown batch effects
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no surrogate variables were found, corresponding part would miss.
#' @seealso \code{\link{isvacor}}, \code{\link{pcacor}},\code{\link{limmafit}},\code{\link{famtcor}}
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' }
#' @export
svacor <- function(data, lv) {
        mod <- stats::model.matrix(~lv)
        svafit <- sva::sva(data, mod)
        if (svafit$n.sv == 0) {
                message("No surrogate variable found")
                li <- limmafit(data,lv)
        } else {
                message("Data is correcting ...")
                batch <- svafit$sv
                li <- limmafit(data,lv,batch)
                message("Done!")
        }
        return(li)
}

#' Use Independent Surrogate Variable Analysis(ISVA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Independent Surrogate Variable Analysis(ISVA) to correct the unknown batch effects
#' @return list object with raw data, corrected data, signal part, batch part, random errors part,  p-values, q-values. If no independent surrogate variable found, corresponding part would miss.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' list <- isvacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}},\code{\link{famtcor}}
#' @export
isvacor <- function(data, lv) {
        isvafit <- isva::DoISVA(data, lv, factor.log = T)
        if (isvafit$nsv == 0) {
                message("No surrogate variable found")
                li <- limmafit(data,lv)
        } else {
                message("Data is correcting ...")
                batch <- isvafit$isv
                li <- limmafit(data,lv,batch)
                message("Done!")
        }
        return(li)
}

#' Use Principal component analysis(PCA) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Principal component analysis(PCA) to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- pcacor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{famtcor}}
#' @export
pcacor <- function(data, lv) {
        batch <- svd(data - rowMeans(data))$v[,1]
        message("Data is correcting ...")
        li <- limmafit(data,lv,batch)
        message("Done!")
        return(li)
}

#' Use Factor Analysis for Multiple Testing(FAMT) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Factor Analysis for Multiple Testing(FAMT) to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, p-values, q-values.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- famtcor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{pcacor}}
#' @export
famtcor <- function(data, lv){                                           covFAMT  <- data.frame(id    = colnames(data),
                                                                                                trmt  = as.factor(lv))

dataFAMT <- FAMT::as.FAMTdata(expression = data, covariates = covFAMT, idcovar    = 1)
fitFAMT  <- FAMT::modelFAMT(dataFAMT,                            x = 2, test = 2)
data <- data
datacor <- fitFAMT$adjdata$expression
pValues <- fitFAMT$pval
qValues <- stats::p.adjust(pValues, method = "BH")
# get the results as list
li <- list(data, datacor, pValues, qValues)
names(li) <- c("data","dataCorrected","p-values", "q-values")
return(li)
}

#' Use Random main effect and Random compound-specific error variance with a mixture structure(RRmix) to correct the unknown batch effects
#' @param data data as mzrt profile
#' @param lv factor vector for the group infomation
#' @details Use Random main effect and Random compound-specific error variance with a mixture structure(RRmix) to correct the unknown batch effects
#' @return list object with various components such raw data, posterior probability.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- rrmixcor(list$data,list$group$class)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{limmafit}},\code{\link{pcacor}}
#' @export
rrmixcor <- function(data, lv){
        n <- RRmix::nfactors(data,plot = F)
        re <- RRmix::runRRmix(t(data),lv,q.in = ifelse(n>10,1,n))
        posterior <- re[['b_g']]
        # get the results as list
        li <- list(data, posterior)
        names(li) <- c("data","posterior")
        return(li)
}
#' Plot ROC curve for simulation data
#' @param sim simulation results from `mzrtsim` or `simmzrt`
#' @param pvalue p value results from `limmafit` function or related batch correction method
#' @param points numbers of points for ROC curve
#' @export
limmaroc <- function(sim, pvalue, points = 100) {
        indexc <- sim$conp
        TPR <- FPR <- FNR <- c(0:points)
        for (i in 0:points) {
                indexi <- which(pvalue < i / points, arr.ind = T)
                TPR[i + 1] <-
                        length(intersect(indexi, indexc)) / length(indexc)
                FPR[i + 1] <-
                        (length(indexi) - length(intersect(indexi, indexc))) / (nrow(sim$data) - length(indexc))
                FNR[i + 1] <-
                        1 - length(intersect(indexi, indexc)) / length(indexc)
        }
        graphics::plot(
                TPR ~ FPR,
                col = 'red',
                pch = 19,
                type = 'l',
                xlim = c(0, 1),
                ylim = c(0, 1)
        )
        graphics::lines(FNR ~ FPR,
                        col = 'blue',
                        pch = 19,
                        type = 'l')

        height = (TPR[-1] + TPR[-length(TPR)]) / 2
        width = diff(FPR)
        auc <- sum(height * width)
        return(auc)
}

#' Principal component analysis(PCA) for limma result with batch
#' @param list results from `limmafit` function or related batch correction method
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv factor vector for the group infomation
#' @return plot
#' @examples
#' \dontrun{
#' library(faahKO)
#' library(enviGCMS)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' li <- svacor(list$data,list$group$class)
#' limmapca(li)
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}}
#' @export
limmapca <- function(list,
                     center = T,
                     scale = T,
                     lv = NULL) {
        data <- list$data
        Signal <- list$signal
        Batch <- list$batch
        error <- list$error
        datacor <- list$dataCorrected
        if (is.null(lv)) {
                pch = colnames(data)
        } else {
                pch = as.character(lv)
        }

        graphics::par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))

        pcao <-
                stats::prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)),
                          3) * 100
        graphics::plot(pcao, type = "l", main = "PCA")

        pca <- stats::prcomp(t(Signal), center = TRUE, scale = TRUE)
        pcaVars = signif(((pca$sdev) ^ 2) / (sum((pca$sdev) ^ 2)),
                         3) * 100
        graphics::plot(pca, type = "l", main = "PCA-signal")

        pcab <-
                stats::prcomp(t(Batch), center = center, scale = scale)
        pcabVars = signif(((pcab$sdev) ^ 2) / (sum((pcab$sdev) ^ 2)),
                          3) * 100
        graphics::plot(pcab, type = "l", main = "PCA-batch")

        pcae <-
                stats::prcomp(t(error), center = center, scale = scale)
        pcaeVars = signif(((pcae$sdev) ^ 2) / (sum((pcae$sdev) ^ 2)),
                          3) * 100
        graphics::plot(pcae, type = "l", main = "PCA-error")

        pcac <- stats::prcomp(t(datacor), center = center,
                              scale = scale)
        pcacVars = signif(((pcac$sdev) ^ 2) / (sum((pcac$sdev) ^ 2)),
                          3) * 100
        graphics::plot(pcac, type = "l", main = "PCA-corrected")

        graphics::plot(
                pcao$x[, 1],
                pcao$x[, 2],
                xlab = paste("PC1:",
                             pcaoVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcaoVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA"
        )

        graphics::plot(
                pca$x[, 1],
                pca$x[, 2],
                xlab = paste("PC1:",
                             pcaVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcaVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-signal"
        )

        graphics::plot(
                pcab$x[, 1],
                pcab$x[, 2],
                xlab = paste("PC1:",
                             pcabVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcabVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-batch"
        )

        graphics::plot(
                pcae$x[, 1],
                pcae$x[, 2],
                xlab = paste("PC1:",
                             pcaeVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcaeVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-error"
        )

        graphics::plot(
                pcac$x[, 1],
                pcac$x[, 2],
                xlab = paste("PC1:",
                             pcacVars[1], "% of Variance Explained"),
                ylab = paste("PC2:",
                             pcacVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-corrected"
        )
}

#' Filter the data with p value and q value and show them as heatmap
#' @param list results from `limmafit` function or related batch correction method
#' @param lv factor vector for the group infomation
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @param index index for selected peaks
#' @return heatmap for the data
#' @examples
#' \dontrun{
#' sim <- mzrtsim()
#' li <- svacor(log(sim$data), as.factor(sim$con))
#' limmaplot(li,as.factor(sim$con))
#' }
#' @seealso \code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}}
#' @export
limmaplot <- function(list,
                      lv,
                      pt = 0.05,
                      qt = 0.05,
                      index = NULL) {
        data <- list$data[, order(lv)]
        signal <- list$signal[, order(lv)]
        batch <- list$batch[, order(lv)]
        error <- list$error[, order(lv)]
        datacor <- list$dataCorrected[, order(lv)]
        pValues <- list$"p-values"
        qValues <- list$"q-values"
        if (!is.null(index)) {
                data <- data[index, order(lv)]
                signal <- signal[index, order(lv)]
                batch <- batch[index, order(lv)]
                error <- error[index, order(lv)]
                datacor <- datacor[index, order(lv)]
                pValues <- pValues[index]
                qValues <- qValues[index]
        }
        # line position
        pos <- cumsum(as.numeric(table(lv) / sum(table(lv)))) -
                as.numeric(table(lv) / sum(table(lv))) / 2
        posv <-
                cumsum(as.numeric(table(lv) / sum(table(lv))))[1:(nlevels(lv) -
                                                                          1)]
        icolors <-
                (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(100)
        # plot the scale
        plotchange <- function(zlim) {
                breaks <- seq(zlim[1], zlim[2], round((zlim[2] -
                                                               zlim[1]) /
                                                              10))
                poly <- vector(mode = "list", length(icolors))
                graphics::plot(
                        1,
                        1,
                        t = "n",
                        xlim = c(0, 1),
                        ylim = zlim,
                        xaxt = "n",
                        yaxt = "n",
                        xaxs = "i",
                        yaxs = "i",
                        ylab = "",
                        xlab = "",
                        frame.plot = F
                )
                graphics::axis(
                        4,
                        at = breaks,
                        labels = round(breaks),
                        las = 1,
                        pos = 0.4,
                        cex.axis = 0.8
                )
                p <- graphics::par("usr")
                graphics::text(
                        p[2] + 1,
                        mean(p[3:4]),
                        labels = "intensity",
                        xpd = NA,
                        srt = -90
                )
                bks <-
                        seq(zlim[1], zlim[2], length.out = (length(icolors) +
                                                                    1))
                for (i in seq(poly)) {
                        graphics::polygon(
                                c(0.1, 0.1, 0.3, 0.3),
                                c(bks[i],
                                  bks[i + 1], bks[i + 1], bks[i]),
                                col = icolors[i],
                                border = NA
                        )
                }
        }
        # plot without batch
        plotimage1 <- function(data, signal, error, zlim) {
                graphics::image(
                        t(data),
                        col = icolors,
                        xlab = "samples",
                        main = "peaks",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(
                        1,
                        at = pos,
                        labels = levels(lv),
                        cex.axis = 0.8
                )
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                data
                        ) -
                                1)),
                        labels = rownames(data),
                        cex.axis = 1,
                        las = 2
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(signal),
                        col = icolors,
                        xlab = "samples",
                        main = "signal",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(
                        1,
                        at = pos,
                        labels = levels(lv),
                        cex.axis = 0.8
                )
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                signal
                        ) -
                                1)),
                        labels = rownames(signal),
                        cex.axis = 1,
                        las = 2
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(error),
                        col = icolors,
                        xlab = "samples",
                        main = "error",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(
                        1,
                        at = pos,
                        labels = levels(lv),
                        cex.axis = 0.8
                )
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                error
                        ) -
                                1)),
                        labels = rownames(error),
                        cex.axis = 1,
                        las = 2
                )
                graphics::abline(v = posv)
        }
        # plot with batch
        plotimage2 <- function(data,
                               signal,
                               batch,
                               error,
                               datacor,
                               zlim) {
                graphics::image(
                        t(data),
                        col = icolors,
                        xlab = "samples",
                        main = "peaks",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                data
                        ) -
                                1)),
                        labels = rownames(data),
                        las = 1
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(signal),
                        col = icolors,
                        xlab = "samples",
                        main = "signal",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(
                        1,
                        at = pos,
                        labels = levels(lv),
                        cex.axis = 1
                )
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                signal
                        ) -
                                1)),
                        labels = rownames(signal),
                        las = 1
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(batch),
                        col = icolors,
                        xlab = "samples",
                        main = "batch",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                batch
                        ) -
                                1)),
                        labels = rownames(batch),
                        las = 1
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(error),
                        col = icolors,
                        xlab = "samples",
                        main = "error",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                error
                        ) -
                                1)),
                        labels = rownames(error),
                        las = 1
                )
                graphics::abline(v = posv)

                graphics::image(
                        t(datacor),
                        col = icolors,
                        xlab = "samples",
                        main = "peaks-corrected",
                        xaxt = "n",
                        yaxt = "n",
                        zlim = zlim
                )
                graphics::axis(1, at = pos, labels = levels(lv))
                graphics::axis(
                        2,
                        at = seq(0, 1, 1 / (nrow(
                                datacor
                        ) -
                                1)),
                        labels = rownames(datacor),
                        las = 1
                )
                graphics::abline(v = posv)

        }

        # plot heatmap
        if (is.null(batch)) {
                if (sum(pValues < pt & qValues <
                        qt) != 0) {
                        message("No batch while p-values and q-values have results")
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 1, 2, 2, 2, 3, 3,
                                  3, 4, 4), 12
                        ), 12, 12, byrow = TRUE))
                        data <- data[pValues < pt & qValues < qt,]
                        signal <-
                                signal[pValues < pt & qValues < qt,]
                        error <- error[pValues < pt & qValues < qt,]

                        zlim <- range(c(data, signal, error))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage1(data, signal, error, zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)

                        li <-
                                list(data, pValues < pt & qValues < qt)
                        names(li) <- c("data", "pqvalues")
                        return(li)
                } else {
                        message("No batch while p-values and q-values have no results")
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 1, 2, 2, 2, 3, 3,
                                  3, 4, 4), 12
                        ), 12, 12,  byrow = TRUE))
                        zlim <- range(c(data, signal, error))
                        graphics::par(mar = c(3, 5, 2, 1))
                        plotimage1(data, signal, error, zlim)
                        graphics::par(mar = c(3, 1, 2, 5))
                        plotchange(zlim)
                }
        } else {
                if (sum(pValues < pt & qValues <
                        qt) != 0) {
                        message("p-values and q-values have results")
                        zlim <- range(c(data, signal, batch, error,
                                        datacor))
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 2, 2,
                                  3, 3, 4, 4, 5, 5, 5, 6, 6), 14
                        ), 14, 14, byrow = TRUE))
                        data <- data[pValues < pt & qValues < qt,]
                        signal <- signal[pValues < pt & qValues <
                                                 qt,]
                        batch <- batch[pValues < pt & qValues <
                                               qt,]
                        error <- error[pValues < pt & qValues <
                                               qt,]
                        datacor <- datacor[pValues < pt & qValues <
                                                   qt,]
                        zlim <- range(c(data, signal, batch, error,
                                        datacor))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage2(data,
                                   signal,
                                   batch,
                                   error,
                                   datacor,
                                   zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)
                        li <- list(datacor, data, pValues < pt &
                                           qValues < qt)
                        names(li) <-
                                c("dataCorrected", "data", "pqvalues")
                        return(li)
                } else {
                        message("p-values and q-values have no results")
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 2, 2,
                                  3, 3, 4, 4, 5, 5, 5, 6, 6), 14
                        ), 14, 14, byrow = TRUE))
                        zlim <- range(c(signal, data, batch, error,
                                        datacor))
                        graphics::par(mar = c(3, 4, 2, 1))
                        plotimage2(data,
                                   signal,
                                   batch,
                                   error,
                                   datacor,
                                   zlim)
                        graphics::par(mar = c(3, 1, 2, 4))
                        plotchange(zlim)
                }
        }
}
