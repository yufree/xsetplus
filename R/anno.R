suppressPackageStartupMessages(library(MAIT))
suppressPackageStartupMessages(library(xMSannotator))
#' Annotation with MAIT package
#' @param path the path to your data
#' @param name the name of your project
#' @param mode defalut is for positive mode, otherwise use 'negAdducts' for negative
#' @return as shown in MAIT package
#' @references Fernández-Albert, F.; Llorach, R.; Andrés-Lacueva, C.; Perera, A. Bioinformatics 2014, 30 (13), 1937–1939.
#' @examples
#' \dontrun{
#' path <- "./data/"
#' name <- "fishriver"
#' anno(path,name)
#' }
#' @export

anno <- function(path, name, mode = NULL) {
        MAIT <- MAIT::sampleProcessing(dataDir = path, project = name)
        MAIT <-
                MAIT::peakAnnotation(
                        MAIT.object = MAIT,
                        corrWithSamp = 0.7,
                        corrBetSamp = 0.75,
                        perfwhm = 0.6,
                        adductTable = mode
                )
        MAIT <-
                MAIT::spectralSigFeatures(
                        MAIT.object = MAIT,
                        pvalue = 0.05,
                        p.adj = "BH",
                        scale = FALSE
                )
        signTable <-
                MAIT::sigPeaksTable(MAIT.object = MAIT, printCSVfile = T)
        MAIT::Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)
        if (mode == 'negAdducts') {
                modename = 'negative'
        } else{
                modename = 'positive'
        }
        MAIT <-
                MAIT::identifyMetabolites(
                        MAIT.object = MAIT,
                        peakTolerance = 0.005,
                        polarity = modename
                )
        metTable <- MAIT::metaboliteTable(MAIT)
        return(list(signTable, metTable))
}

#' Annotation with xMSannotator package
#' @param xset a xcmsset object to be annotated
#' @param outloc the path for your result
#' @param mode defalut is for positive mode, otherwise use 'neg' for negative
#' @param db_name default is 'HMDB', other database options: 'KEGG', 'LipidMaps', 'T3DB'
#' @param num_nodes default 10
#' @param ppm default 5 for mass accuracy
#' @param ... parameters for multilevelannotation function in xMSannotator
#' @return as shown in xMSannotator package
#' @references Uppal, K.; Walker, D. I.; Jones, D. P. Anal. Chem. 2017, 89 (2), 1063–1067.
#' @examples
#' \dontrun{
#' path <- "./data/"
#' xset <- getdata(path)
#' result <- fanno(xset)
#' }
#' @export

fanno <-
        function(xset,
                 outloc = "./result/",
                 mode = 'pos',
                 db_name = 'HMDB', num_nodes = 10,ppm=5,...) {
                data <- xcms::groupval(xset, 'medret', "into")
                adduct_weights = cbind.data.frame(Adduct = c('M+H','M-H'),Weight = c(5,5))
                mz <- xcms::groups(xset)[, 1]
                time <- xcms::groups(xset)[, 4]
                data <- as.data.frame(cbind(mz, time, data))
                data <- unique(data)
                if ( mode == 'neg') {
                        annotres <-
                                xMSannotator::multilevelannotation(
                                        dataA = data,
                                        max.mz.diff = ppm,
                                        mode = mode,
                                        outloc = outloc,
                                        db_name = db_name,
                                        adduct_weights = adduct_weights,
                                        filter.by = c("M-H"),
                                        mass_defect_mode = mode,
                                        num_nodes = num_nodes,
                                        ...
                                )
                }else{
                        annotres <-
                                xMSannotator::multilevelannotation(
                                        dataA = data,
                                        max.mz.diff = ppm,
                                        mode = mode,
                                        outloc = outloc,
                                        db_name = db_name,
                                        adduct_weights = adduct_weights,
                                        filter.by = c("M+H"),
                                        mass_defect_mode = mode,
                                        num_nodes = num_nodes,
                                        ...
                                )
                }
                return(annotres)
        }
#' Annotation with xMSannotator package for list
#' @param list a list to be annotated
#' @param outloc the path for your result
#' @param mode defalut is for positive mode, otherwise use 'neg' for negative
#' @param db_name default is 'HMDB', other database options: 'KEGG', 'LipidMaps', 'T3DB'
#' @param num_nodes default 10
#' @param ppm for mass accuracy
#' @param ... parameters for multilevelannotation function in xMSannotator
#' @return as shown in xMSannotator package
#' @references Uppal, K.; Walker, D. I.; Jones, D. P. Anal. Chem. 2017, 89 (2), 1063–1067.
#' @examples
#' \dontrun{
#' path <- "./data/"
#' xset <- getdata(path)
#' list <- getmzrt(xset)
#' result <- fanno2(list)
#' }
#' @export
fanno2 <-
        function(list,
                 outloc = "./result/",
                 mode = 'pos',
                 db_name = 'HMDB', num_nodes = 10,ppm=5,...) {
                data <- list$data
                mz <- list$mz
                time <- list$rt
                adduct_weights = cbind.data.frame(Adduct = c('M+H','M-H'),Weight = c(5,5))
                data <- as.data.frame(cbind(mz, time, data))
                data <- unique(data)
                if ( mode == 'neg') {
                        annotres <-
                                xMSannotator::multilevelannotation(
                                        dataA = data,
                                        max.mz.diff = ppm,
                                        mode = mode,
                                        outloc = outloc,
                                        db_name = db_name,
                                        adduct_weights = adduct_weights,
                                        filter.by = c("M-H"),
                                        mass_defect_mode = mode,
                                        num_nodes = num_nodes,
                                        ...
                                )
                }else{
                        annotres <-
                                xMSannotator::multilevelannotation(
                                        dataA = data,
                                        mode = mode,
                                        max.mz.diff = ppm,
                                        outloc = outloc,
                                        db_name = db_name,
                                        adduct_weights = adduct_weights,
                                        filter.by = c("M+H"),
                                        mass_defect_mode = mode,
                                        num_nodes = num_nodes,
                                        ...
                                )
                }
                return(annotres)
        }
