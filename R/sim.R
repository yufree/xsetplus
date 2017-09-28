#' Generate simulated count data with batch effects for npeaks
#'
#' @param npeaks Number of genes to simulate
#' @param nbatch Number of batches to simulate
#' @param ncond Number of conditions to simulate
#' @param npercond Number of samples per condition per batch to simulate
#' @param basemean Base mean
#' @param ppstep peak to peak step variation
#' @param bbstep Batch to Batch step variation
#' @param ccstep Condition to Condition step variation
#' @param basedisp Base Dispersion
#' @param bdispstep Batch to Batch Dispersion step variation
#' @param swvar Sample-wise extra variation
#' @param seed Random seed for reproducibility
#' @return rtmz data matrix
#' @export
#' @examples
#' rtmzsim()
rtmzsim <- function(npeaks = 100, nbatch = 3, ncond = 2, npercond = 10,
                       basemean = 10000, ppstep = 50, bbstep = 2000, ccstep = 800,
                       basedisp = 100, bdispstep = 10, swvar = 1000, seed = 42) {
        set.seed(seed)
        mu <- seq(0, length.out = npeaks, by = ppstep)
        bmu <- seq(0, length.out = nbatch, by = bbstep)
        cmu <- seq(0, length.out = ncond, by = ccstep)
        bsize <- seq(basedisp, length.out = nbatch, by = bdispstep)
        ncol <- nbatch * ncond * npercond
        A.matrix <- matrix(0, nrow = npeaks, ncol = ncol)
        samplewisevar <- swvar*rbeta(ncol,2,2)
        for (i in 1:npeaks) {
                peaki <- c()
                for (j in 1:nbatch) {
                        for (k in 1:ncond) {
                                for (l in 1:npercond) {
                                        peaki <- c(peaki, rnbinom(1, size = bsize[j],
                                                                  mu = basemean+mu[i]+bmu[j]+cmu[k]+samplewisevar[j*k*l]))
                                }
                        }
                }
                A.matrix[i, ] <- peaki
        }
        return(A.matrix)
}
