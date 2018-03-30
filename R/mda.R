#' Group peaks by the mass defect interval group for different substructures
#' @param list a peaks list with mass to charge
#' @param submass mass vector of sub structure of homologous series
#' @param mdgn mass defect groups numbers for interval, 20 means 0.05 inteval on mass defect scale from -0.5 to 0.5
#' @param lv group info for the data
#' @return list with mass defect analysis dataframe.
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{plotstd}},\code{\link{plotstdmd}},\code{\link{plotstdrt}}
#' @export
getmdg <- function(list, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdgn = 20, lv = NULL){
        if(is.null(list$stdmass)&is.null(list$paired)){
                mz <- list$mz
                rt <- list$rt
                data <- list$data
                colnames(data) <- lv
        }else if(is.null(list$stdmass)){
                mz <- list$mz[list$pairedindex]
                rt <- list$rt[list$pairedindex]
                data <- list$data[list$pairedindex,]
                colnames(data) <- lv
        }else{
                mz <- list$mz[list$stdmassindex]
                rt <- list$rt[list$stdmassindex]
                data <- list$data[list$stdmassindex,]
                colnames(data) <- lv
        }
        # perform mass defect analysis for std mass
        mda <- cbind.data.frame(mz = mz, rt = rt, data)
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(mz*mdst) - mz*mdst

                mdg <- cut(msdefect, seq(from = -.5, to = .5, by = 1/mdgn),include.lowest = T)
                mdg2 <- mdg[!is.na(mdg)]
                index <- mdg %in% Mode(mdg2)

                name <- c(submass[i],paste0(submass[i],'gi'), 'majormdg')
                md <- cbind.data.frame(msdefect,mdg,index)
                colnames(md) <- name
                mda <- cbind.data.frame(mda,md)
        }
        # get the data
        list$mda <- mda
        return(list)
}
#' Paired correlationship among peak list based on cluster analysis
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in cluster
#' @param submass mass vector of sub structure of homologous series
#' @param mdcutoff mass defect cluster cutoff
#' @return list with retention time cluster, std mass defect analysis dataframe based on max average correlation
getcorstd <- function(list, rtcutoff = 9, submass = c(15.9949,14.003074,26.01568,14.01565,43.00581,30.01056,34.96885,78.91834), mdcutoff = 0.02){

        # paired mass diff analysis
        groups <- cbind.data.frame(mz = list$mz, rt = list$rt, list$data)
        resultstd <- resultsolo <- resultiso <- result <- NULL

        dis <- stats::dist(list$rt, method = "manhattan")
        fit <- stats::hclust(dis)
        rtcluster <- stats::cutree(fit, h=rtcutoff)
        n <- length(unique(rtcluster))
        message(paste(n, 'retention time cluster found.'))
        # search:
        for (i in 1:length(unique(rtcluster))) {
                # find the mass within RT
                rtxi <- list$rt[rtcluster == i]
                bin = groups[groups$rt %in% rtxi, ]
                medianrtxi <- stats::median(rtxi)

                if (nrow(bin) > 1) {
                        # get mz diff
                        cor <- stats::cor(t(bin[,-c(1,2)]))
                        cormean <- apply(cor,1,mean)
                        corindex <- which.max(cormean)
                        df <- cbind(bin[corindex,],rtg = i)
                        result <- rbind(result,df)
                }else{
                        solo <- cbind(bin,rtg = i)
                        resultsolo <- rbind(solo,resultsolo)
                }
        }

        resultstd <- rbind(result,resultsolo)
        resultstd <- unique(resultstd)

        # perform mass defect analysis for std mass
        for(i in 1:length(submass)){
                mdst <- round(submass[i])/submass[i]
                msdefect <- round(resultstd$mz*mdst) - resultstd$mz*mdst
                dis <- stats::dist(msdefect, method = "manhattan")
                fit <- stats::hclust(dis)
                mdcluster <- stats::cutree(fit, h=mdcutoff)
                n <- length(unique(mdcluster))
                message(paste(n, 'mass defect clusters found for mass', submass[i], 'substructures' ))
                name <- c(colnames(resultstd),submass[i],paste0(submass[i],'g'))
                resultstd <- cbind.data.frame(resultstd,msdefect,mdcluster)
                colnames(resultstd) <- name
        }

        # filter the list
        list$rtcluster <- rtcluster
        list$stdmassindex <- (round(list$mz,4) %in% round(resultstd$mz,4))
        list$stdmass <- resultstd
        return(list)
}
