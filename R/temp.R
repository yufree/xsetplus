# Function for SPME obitrap
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
# Function for SPME qToF
getopqtofdata <- function(path,
                          index = F,
                          xsmethod = "centWave",
                          peakwidth = c(6.75, 46.5),
                          ppm = 29,
                          noise = 0,
                          snthresh = 10,
                          mzdiff = 0.00692,
                          prefilter = c(0, 0),
                          mzCenterFun = "wMean",
                          integrate = 1,
                          fitgauss = FALSE,
                          verbose.columns = FALSE,
                          BPPARAM = BiocParallel::SnowParam(workers = 12),
                          rmethod = "obiwarp",
                          plottype       = "none",
                          distFunc       = "cor_opt",
                          profStep       = 1,
                          center         = 1,
                          response       = 1,
                          gapInit        = 0.48,
                          gapExtend      = 2.4,
                          factorDiag     = 2,
                          factorGap      = 1,
                          localAlignment = 0,
                          gmethod = "density",
                          bw      = 1.74,
                          mzwid   = 0.01905,
                          minfrac = 1,
                          minsamp = 1,
                          gmax     = 20,
                          lockmass = T,
                          ...){
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
                lockMassFreq = lockmass,
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
#' Output pmd results into csv(s) for MS/MS validation
#' @param list list from getstd or globalstd with independant peaks index
#' @param Dppm ppm for targeted analysis
#' @param Drt retention time drift for targeted analysis
#' @param ce collision energy for MS/MS analysis
#' @param name file name
#' @param n max ions numbers within retention time drift windows
#' @return NULL
#' @export
pmdtarget <- function(list,Dppm = 20,Drt = 0.5,ce = NA, name = 'target',n=NULL){
        head <-  c('On', 'Prec. m/z', 'Delta m/z (ppm)','Z', 'Prec. Type', 'Ret. Time (min)', 'Delta Ret. Time (min)', 'Iso. Width', 'Collision Energy')
        mz <- list$mz[list$stdmassindex]
        rt <- round(list$rt[list$stdmassindex]/60,3)
        temp = cbind('TRUE',mz,Dppm,1,'Preferred',rt,Drt,'Narrow (~1.3 m/z)',ce)
        data <- rbind(head,temp)
        colnames(data) <- c('AutoPreferredExcludeMSMSTable',rep('',8))

        if(is.null(n)){
                name2 <- paste0(name,'.csv')
                utils::write.csv(data,file = name2,row.names = F)

        }else{
                idx <- targetsep(list$rt[list$stdmassindex],Drt,n)
                for(i in 1:length(table(idx))){
                        namei <- paste0(name,i,'.csv')
                        idx2 <- idx == i
                        idx3 <- c(T,idx2)
                        datai <- data[idx3,]
                        utils::write.csv(datai,file = namei,row.names = F)
                }
        }

        return(data)
}

#' Output pmd results into csv(s) for MS/MS validation
#' @param rt retention time vector for peaks
#' @param Drt retention time drift for targeted analysis
#' @param n max ions numbers within retention time drift windows
#' @return index for each injection
#' @export
targetsep <- function(rt,Drt,n=6){
        D <- Drt*60
        dis <- stats::dist(rt, method = "manhattan")
        fit <- stats::hclust(dis)
        inji <- rtcluster <- stats::cutree(fit, h = D)
        maxd <- max(table(rtcluster))
        m <- length(unique(rtcluster))
        inj <- ceiling(maxd/n)
        message(paste('You need',inj,'injections!'))
        for(i in c(1:m)) {
                z = 1:inj
                x <- rt[rtcluster==i]
                while(length(x) > inj & length(x)>n){
                        t <- sample(x,n)
                        w <- sample(z,1)
                        inji[rt %in% t] <- w
                        z <- z[!(z%in%w)]
                        x <- x[!(x %in% t)]
                }
                inji[rtcluster==i & rt %in% x] <- sample(z,sum(rtcluster==i & rt %in% x),replace = T)
        }
        return(inji)
}
#' Plot circle connection of metabolites from list
#' @param list list with data as peaks list, mz, rt and group information
#' @param class metabolites or compounds class, default NULL
#' @param cutoff cutoff of the absolute value of correlation coefficient
#' @return NULL
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' list <- getmr(cdfpath, pmethod = ' ')
#' data <- list$data
#' lv <- as.character(list$group$class)
#' plotcon(list, lv)
#' }
#' @export

plotcon <- function(list,class=NULL,cutoff = 0.6){
        cmat <- stats::cor(t(list$data))
        compounds <- rownames(list$data)
        if(is.null(class)){
                edges=data.frame(from="origin", to=compounds)
        }else{
                d1=data.frame(from="origin", to=class)
                d2=data.frame(from=class, to=compounds)
                edges=rbind(d1, d2)
        }
        vertices <- cbind.data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
        vertices$group = edges$from[ match( vertices$name, edges$to ) ]

        df <- data.frame(from = compounds[which(lower.tri(cmat), arr.ind = T)[, 1]],
                         to = compounds[which(lower.tri(cmat),arr.ind = T)[, 2]],cor = c(cmat[lower.tri(cmat)]))

        dft <- df[abs(df$cor)>cutoff,]

        from = match( dft$from, vertices$name)
        to = match( dft$to, vertices$name)

        vertices$id=NA
        myleaves=which(is.na( match(vertices$name, edges$from) ))
        nleaves=length(myleaves)
        vertices$id[ myleaves ] = seq(1:nleaves)
        vertices$angle= 90 - 360 * vertices$id / nleaves

        vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
        vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

        mygraph <- igraph::graph_from_data_frame(edges, vertices=vertices )
        ggraph::ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
                ggraph::geom_conn_bundle(data = ggraph::get_con(from = from, to = to, value = abs(dft$cor)), alpha=0.07, colour="skyblue",ggplot2::aes(width=value))  +
                ggraph::geom_node_point(ggplot2::aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, alpha=0.2)) +
                ggraph::geom_node_text(ggplot2::aes(x = x*1.1, y=y*1.1, filter = leaf, label=name, angle = angle, hjust=hjust), size=1.2, alpha=1)+
                ggplot2::theme_void() +
                ggplot2::theme(
                        legend.position="none",
                        plot.margin=grid::unit(c(0,0,0,0),"cm")
                )+
                ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

}
