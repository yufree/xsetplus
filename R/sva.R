#' Surrogate variable analysis(SVA) to correct the unknown batch effects
#' @param xset xcmsset object
#' @param lv group information
#' @param annotation logical annotation with CAMERA package
#' @param polarity 'positive' or 'negative' for annotation
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @param nSlaves parameter for CAMARA, although xcms depute it
#' @details this is used for reviesed version of SVA to correct the unknown batch effects
#' @return list object with various components such raw data, corrected data, signal part, random errors part, batch part, p-values, q-values, mass, rt, Posterior Probabilities of Surrogate variables and Posterior Probabilities of Mod. If no surrogate variable found, corresponding part would miss.
#' @export
svacor <-
        function(xset,
                 lv = NULL,
                 annotation = F,
                 polarity = "positive",
                 method  = "medret",
                 intensity = 'into',
                 nSlaves = 12) {
                data <- xcms::groupval(xset, method, intensity)
                if (intensity == "intb") {
                        data[is.na(data)] = 0
                }
                if (is.null(lv)) {
                        lv <- xset@phenoData[, 1]
                }
                mz <- xset@groups[, 1]
                rt <- xset@groups[, 4]
                mod <- stats::model.matrix(~ lv)
                mod0 <- as.matrix(c(rep(1, ncol(data))))
                svafit <- sva::sva(data, mod)
                if (svafit$n.sv == 0) {
                        svaX <- stats::model.matrix(~ lv)
                        lmfit <- limma::lmFit(data, svaX)
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal
                        rownames(signal) <-
                                rownames(error) <- rownames(data)
                        colnames(signal) <-
                                colnames(error) <- colnames(data)
                        pValues = sva::f.pvalue(data, mod, mod0)
                        qValues = qvalue::qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        CAMERA::annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))),]
                                li <-
                                        list(data,
                                             signal,
                                             error,
                                             pValues,
                                             qValues,
                                             dreport,
                                             mz,
                                             rt)
                                names(li) <-
                                        c(
                                                'data',
                                                'signal',
                                                'error',
                                                'p-values',
                                                'q-values',
                                                'diffreport',
                                                'mz',
                                                'rt'
                                        )
                        } else{
                                li <- list(data,
                                           signal,
                                           error,
                                           pValues,
                                           qValues,
                                           mz,
                                           rt)
                                names(li) <-
                                        c(
                                                'data',
                                                'signal',
                                                'error',
                                                'p-values',
                                                'q-values',
                                                'mz',
                                                'rt'
                                        )
                        }
                }
                else{
                        message('Data is correcting ...')
                        svaX <- stats::model.matrix(~ lv + svafit$sv)
                        lmfit <- limma::lmFit(data, svaX)
                        batch <-
                                lmfit$coef[, (nlevels(lv) + 1):(nlevels(lv) + svafit$n.sv)] %*% t(svaX[, (nlevels(lv) +
                                                                                                                  1):(nlevels(lv) + svafit$n.sv)])
                        signal <-
                                lmfit$coef[, 1:nlevels(lv)] %*% t(svaX[, 1:nlevels(lv)])
                        error <- data - signal - batch
                        datacor <- signal + error
                        svaX2 <- stats::model.matrix(~ lv)
                        lmfit2 <- limma::lmFit(data, svaX2)
                        signal2 <-
                                lmfit2$coef[, 1:nlevels(lv)] %*% t(svaX2[, 1:nlevels(lv)])
                        error2 <- data - signal2
                        rownames(signal2) <-
                                rownames(error2) <-
                                rownames(datacor) <-
                                rownames(signal) <-
                                rownames(batch) <-
                                rownames(error) <- rownames(data)
                        colnames(signal2) <-
                                colnames(error2) <-
                                colnames(datacor) <-
                                colnames(signal) <-
                                colnames(batch) <-
                                colnames(error) <- colnames(data)

                        modSv = cbind(mod, svafit$sv)
                        mod0Sv = cbind(mod0, svafit$sv)
                        pValuesSv = sva::f.pvalue(data, modSv, mod0Sv)
                        qValuesSv = qvalue::qvalue(pValuesSv)
                        qValuesSv = qValuesSv$qvalues

                        pValues = sva::f.pvalue(data, mod, mod0)
                        qValues = qvalue::qvalue(pValues)
                        qValues = qValues$qvalues
                        if (annotation) {
                                dreport <-
                                        CAMERA::annotateDiffreport(
                                                xset,
                                                metlin = T,
                                                polarity = polarity,
                                                nSlaves = nSlaves
                                        )
                                dreport <-
                                        dreport[order(as.numeric(rownames(dreport))),]
                                li <-
                                        list(
                                                data,
                                                datacor,
                                                signal,
                                                batch,
                                                error,
                                                signal2,
                                                error2,
                                                pValues,
                                                qValues,
                                                pValuesSv,
                                                qValuesSv,
                                                dreport,
                                                svafit$pprob.gam,
                                                svafit$pprob.b,
                                                mz,
                                                rt
                                        )
                                names(li) <-
                                        c(
                                                'data',
                                                'dataCorrected',
                                                'signal',
                                                'batch',
                                                'error',
                                                'signal2',
                                                'error2',
                                                'p-values',
                                                'q-values',
                                                'p-valuesCorrected',
                                                'q-valuesCorrected',
                                                'diffreport',
                                                'PosteriorProbabilitiesSurrogate',
                                                'PosteriorProbabilitiesMod',
                                                'mz',
                                                'rt'
                                        )
                        }
                        else{
                                li <-
                                        list(
                                                data,
                                                datacor,
                                                signal,
                                                batch,
                                                error,
                                                signal2,
                                                error2,
                                                pValues,
                                                qValues,
                                                pValuesSv,
                                                qValuesSv,
                                                svafit$pprob.gam,
                                                svafit$pprob.b,
                                                mz,
                                                rt
                                        )
                                names(li) <-
                                        c(
                                                'data',
                                                'dataCorrected',
                                                'signal',
                                                'batch',
                                                'error',
                                                'signal2',
                                                'error2',
                                                'p-values',
                                                'q-values',
                                                'p-valuesCorrected',
                                                'q-valuesCorrected',
                                                'PosteriorProbabilitiesSurrogate',
                                                'PosteriorProbabilitiesMod',
                                                'mz',
                                                'rt'
                                        )
                        }
                        message('Done!')
                }
                return(li)
        }
#' Principal component analysis(PCA) for SVA corrected data and raw data
#' @param list results from svacor function
#' @param center parameters for PCA
#' @param scale parameters for scale
#' @param lv group information
#' @return plot
#' @export
svapca <- function(list,
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
        } else{
                pch = lv
        }

        graphics::par(mfrow = c(2, 5), mar = c(4, 4, 2.6, 1))

        pcao <- stats::prcomp(t(data), center = center, scale = scale)
        pcaoVars = signif(((pcao$sdev) ^ 2) / (sum((pcao$sdev) ^ 2)), 3) *
                100
        graphics::plot(pcao, type = "l", main = "PCA")

        pca <- stats::prcomp(t(Signal), center = TRUE, scale = TRUE)
        pcaVars = signif(((pca$sdev) ^ 2) / (sum((pca$sdev) ^ 2)), 3) *
                100
        graphics::plot(pca, type = "l", main = "PCA-signal")

        pcab <- stats::prcomp(t(Batch), center = center, scale = scale)
        pcabVars = signif(((pcab$sdev) ^ 2) / (sum((pcab$sdev) ^ 2)), 3) *
                100
        graphics::plot(pcab, type = "l", main = "PCA-batch")

        pcae <- stats::prcomp(t(error), center = center, scale = scale)
        pcaeVars = signif(((pcae$sdev) ^ 2) / (sum((pcae$sdev) ^ 2)), 3) *
                100
        graphics::plot(pcae, type = "l", main = "PCA-error")

        pcac <- stats::prcomp(t(datacor), center = center, scale = scale)
        pcacVars = signif(((pcac$sdev) ^ 2) / (sum((pcac$sdev) ^ 2)), 3) *
                100
        graphics::plot(pcac, type = "l", main = "PCA-corrected")

        graphics::plot(
                pcao$x[, 1],
                pcao$x[, 2],
                xlab = paste("PC1:", pcaoVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaoVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA"
        )

        graphics::plot(
                pca$x[, 1],
                pca$x[, 2],
                xlab = paste("PC1:", pcaVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-signal"
        )

        graphics::plot(
                pcab$x[, 1],
                pcab$x[, 2],
                xlab = paste("PC1:", pcabVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcabVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-batch"
        )

        graphics::plot(
                pcae$x[, 1],
                pcae$x[, 2],
                xlab = paste("PC1:", pcaeVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcaeVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-error"
        )

        graphics::plot(
                pcac$x[, 1],
                pcac$x[, 2],
                xlab = paste("PC1:", pcacVars[1], "% of Variance Explained"),
                ylab = paste("PC2:", pcacVars[2], "% of Variance Explained"),
                pch = pch,
                cex = 2,
                main = "PCA-corrected"
        )
}
#' Filter the data with p value and q value
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @return data, corrected data, mz and retention for fileted data
#' @export
svadata <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05) {
        data <- list$data
        signal2 <- list$signal2
        datacor <- list$dataCorrected
        pValues <- list$'p-values'
        qValues <- list$'q-values'
        pValuesSv <- list$'p-valuesCorrected'
        qValuesSv <- list$'q-valuesCorrected'
        mz <- list$mz
        rt <- list$rt
        if (is.null(signal2)) {
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('No SV while p-values and q-values have results')
                        data <- data[pValues < pt & qValues < qt,]
                        mz <- mz[pValues < pt &
                                         qValues < qt]
                        rt <- rt[pValues < pt &
                                         qValues < qt]
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt, mz, rt)
                        names(li) <-
                                c('data', 'pqvalues', 'mz', 'rt')
                        return(li)
                }
                else{
                        message('No SV while p-values and q-values have no results')
                }
        } else{
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('Have SVs while p-values and q-values have results')
                        data <- data[pValues < pt & qValues < qt,]
                        mz <- mz[pValues < pt &
                                         qValues < qt]
                        rt <- rt[pValues < pt &
                                         qValues < qt]
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt, mz, rt)
                        names(li) <- c('data', 'pqvalues', mz, rt)
                        return(li)
                }
                else if (pqvalues == "anova") {
                        message('Have SVs while p-values and q-values have no results')
                }
                else if (pqvalues == "sv" &
                         sum(pValuesSv < pt &
                             qValuesSv < qt) != 0) {
                        message('SVs corrected while p-values and q-values have results')
                        data <-
                                data[pValuesSv < pt &
                                             qValuesSv < qt,]
                        datacor <-
                                datacor[pValuesSv < pt &
                                                qValuesSv < qt,]
                        mz <- mz[pValuesSv < pt &
                                         qValuesSv < qt]
                        rt <- rt[pValuesSv < pt &
                                         qValuesSv < qt]
                        li <-
                                list(datacor,
                                     data,
                                     pValuesSv < pt &
                                             qValuesSv < qt,
                                     mz,
                                     rt)
                        names(li) <-
                                c('dataCorrected',
                                  'data',
                                  'pqvalues',
                                  'mz',
                                  'rt')
                        return(li)
                }
                else{
                        message('SVs corrected while p-values and q-values have no results')
                }
        }
}
#' Filter the data with p value and q value and show them
#' @param list results from svacor function
#' @param pqvalues method for ANOVA or SVA
#' @param pt threshold for p value, default is 0.05
#' @param qt threshold for q value, default is 0.05
#' @return heatmap for the data
#' @export
svaplot <- function(list,
                    pqvalues = "sv",
                    pt = 0.05,
                    qt = 0.05) {
        data <- list$data
        signal <- list$signal
        signal2 <- list$signal2
        batch <- list$batch
        error <- list$error
        error2 <- list$error2
        datacor <- list$dataCorrected
        pValues <- list$'p-values'
        qValues <- list$'q-values'
        pValuesSv <- list$'p-valuesCorrected'
        qValuesSv <- list$'q-valuesCorrected'

        icolors <-
                grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100)

        if (is.null(signal2)) {
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('No SV while p-values and q-values have results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        graphics::par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt,]
                        signal <-
                                signal[pValues < pt &
                                               qValues < qt,]
                        error <-
                                error[pValues < pt & qValues < qt,]
                        zlim <- range(c(data, signal, error))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )

                        graphics::image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )

                        graphics::image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else{
                        message('No SV while p-values and q-values have no results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        graphics::par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal, error))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::par(mar = c(3, 0, 2, 3))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        } else{
                if (pqvalues == "anova" & sum(pValues < pt & qValues < qt) != 0) {
                        message('Have SVs while p-values and q-values have results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5), 9
                        ), 9, 9, byrow = TRUE))
                        graphics::par(mar = c(3, 5, 1, 1))
                        data <- data[pValues < pt & qValues < qt,]
                        signal <-
                                signal2[pValues < pt &
                                                qValues < qt,]
                        error <-
                                error2[pValues < pt &
                                               qValues < qt,]
                        zlim <- range(c(data, signal, error))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )

                        graphics::image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )

                        graphics::image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 2
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(data, pValues < pt &
                                             qValues < qt)
                        names(li) <- c('data', 'pqvalues')
                        return(li)
                }
                else if (pqvalues == "anova") {
                        message('Have SVs while p-values and q-values have no results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4), 8
                        ), 8, 8, byrow = TRUE))
                        graphics::par(mar = c(3, 6, 2, 3))
                        zlim <- range(c(data, signal2, error2))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(signal2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(signal2) - 1
                                )),
                                labels = colnames(signal2),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(error2),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error2
                                ) - 1)),
                                labels = colnames(error2),
                                cex.axis = 0.618,
                                las = 2
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::par(mar = c(3, 0, 2, 3))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
                else if (pqvalues == "sv" &
                         sum(pValuesSv < pt &
                             qValuesSv < qt) != 0) {
                        message('SVs corrected while p-values and q-values have results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6), 11
                        ), 11, 11, byrow = TRUE))
                        graphics::par(mar = c(3, 4, 2, 1))
                        data <-
                                data[pValuesSv < pt &
                                             qValuesSv < qt,]
                        signal <- signal[pValuesSv < pt &
                                                 qValuesSv < qt,]
                        batch <-
                                batch[pValuesSv < pt &
                                              qValuesSv < qt,]
                        error <-
                                error[pValuesSv < pt &
                                              qValuesSv < qt,]
                        datacor <-
                                datacor[pValuesSv < pt &
                                                qValuesSv < qt,]
                        zlim <-
                                range(c(data, signal, batch, error, datacor))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 1
                        )

                        graphics::image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        signal
                                ) - 1)),
                                labels = rownames(signal),
                                cex.axis = 0.618,
                                las = 1
                        )

                        graphics::image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        batch
                                ) - 1)),
                                labels = rownames(batch),
                                cex.axis = 0.618,
                                las = 1
                        )

                        graphics::image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        error
                                ) - 1)),
                                labels = rownames(error),
                                cex.axis = 0.618,
                                las = 1
                        )

                        graphics::image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (
                                        nrow(datacor) - 1
                                )),
                                labels = rownames(datacor),
                                cex.axis = 0.618,
                                las = 1
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::par(mar = c(3, 0, 2, 3))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                        li <-
                                list(datacor,
                                     data,
                                     pValuesSv < pt &
                                             qValuesSv < qt)
                        names(li) <-
                                c('dataCorrected', 'data', 'pqvalues')
                        return(li)
                }
                else{
                        message('SVs corrected while p-values and q-values have no results')
                        graphics::layout(matrix(rep(
                                c(1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6), 13
                        ), 13, 13, byrow = TRUE))
                        graphics::par(mar = c(3, 6, 2, 3))
                        zlim <-
                                range(c(signal, data, batch, error, datacor))

                        graphics::image(
                                t(data),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        data
                                ) - 1)),
                                labels = colnames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::axis(
                                2,
                                at = seq(0, 1, 1 / (nrow(
                                        data
                                ) - 1)),
                                labels = rownames(data),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(signal),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-signal',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        signal
                                ) - 1)),
                                labels = colnames(signal),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(batch),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-batch',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        batch
                                ) - 1)),
                                labels = colnames(batch),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 3, 2, 1))
                        graphics::image(
                                t(error),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-error',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (ncol(
                                        error
                                ) - 1)),
                                labels = colnames(error),
                                cex.axis = 0.618,
                                las = 2
                        )
                        graphics::par(mar = c(3, 6, 2, 3))
                        graphics::image(
                                t(datacor),
                                col = icolors,
                                xlab = 'samples',
                                main = 'peaks-corrected',
                                xaxt = "n",
                                yaxt = "n",
                                zlim = zlim
                        )
                        graphics::axis(
                                1,
                                at = seq(0, 1, 1 / (
                                        ncol(datacor) - 1
                                )),
                                labels = colnames(datacor),
                                cex.axis = 0.618,
                                las = 2
                        )

                        breaks <-
                                seq(zlim[1], zlim[2], round((zlim[2] - zlim[1]) / 10))
                        poly <-
                                vector(mode = "list", length(icolors))
                        graphics::par(mar = c(3, 0, 2, 3))
                        graphics::plot(
                                1,
                                1,
                                t = "n",
                                xlim = c(0, 1),
                                ylim = zlim,
                                xaxt = 'n',
                                yaxt = 'n',
                                xaxs = "i",
                                yaxs = "i",
                                ylab = '',
                                xlab = 'intensity',
                                frame.plot = F
                        )
                        graphics::axis(
                                4,
                                at = breaks,
                                labels = round(breaks),
                                las = 1,
                                pos = 0.4
                        )
                        bks <-
                                seq(zlim[1], zlim[2], length.out = (length(icolors) + 1))
                        for (i in seq(poly)) {
                                graphics::polygon(
                                        c(0.1, 0.1, 0.3, 0.3),
                                        c(bks[i], bks[i + 1], bks[i + 1], bks[i]),
                                        col = icolors[i],
                                        border = NA
                                )
                        }
                }
        }
}

#' Get the corrected data after SVA for metabolanalyst
#' @param xset xcmsset object
#' @param lv group information
#' @return csv files for both raw and corrected data for metabolanalyst if SVA could be applied
#' @export
svaupload <- function(xset,
                      lv = NULL) {
        raw <- svacor(xset,
                      lv = lv,
                      annotation = F,
                      nSlaves = 12)
        if (is.null(raw$dataCorrected)) {
                data <- raw$data
                data <-
                        rbind(group = as.character(xset@phenoData[, 1]), data)
                utils::write.csv(data, file = 'Peaklist.csv')
                return(data)
        }
        else{
                datacor <- raw$dataCorrected
                data <- raw$data
                data <-
                        rbind(group = as.character(xset@phenoData[, 1]), data)
                utils::write.csv(datacor, file = 'Peaklistcor.csv')
                utils::write.csv(data, file = 'Peaklist.csv')
                return(list(data, datacor))
        }
}
