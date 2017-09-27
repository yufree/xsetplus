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
                 db_name = 'HMDB', num_nodes = 10,...) {
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

#' Output pathway annotation data for Mummichog algorithm
#' @param xset a xcmsset object to be annotated
#' @param lv the group factor
#' @param name the file name
#' @param method parameter for groupval function
#' @param intensity parameter for groupval function
#' @details use 'python2 main.py -c 0.05 -f test.txt -p 100 -m negative -o myoutput' to get the Mummichog annotation
#' @return text file for Mummichog algorithm
#' @references Li, S.; Park, Y.; Duraisingham, S.; Strobel, F. H.; Khan, N.; Soltow, Q. A.; Jones, D. P.; Pulendran, B. PLOS Computational Biology 2013, 9 (7), e1003123.

#' @export
mumdata <-
        function(xset,
                 lv = NULL,
                 name = 'test',
                 method = "medret",
                 intensity = 'inio') {
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
                fstats <- sva::fstats(data, mod, mod0)
                pvalue <- sva::f.pvalue(data, mod, mod0)
                df <- cbind.data.frame(mz, rt, pvalue, fstats)
                filename <- paste0(name, '.txt')
                utils::write.table(df,
                            file = filename,
                            sep = "\t",
                            row.names = F)
                return(df)
        }

#' Annotation with MAIT package of list from svacor function
#' @param raw the list from svacor function
#' @param lv group information
#' @param polarity defalut is for positive mode, otherwise use 'negative' for negative
#' @param projectname the name of your project
#' @return as shown in MAIT package
#' @references Fernández-Albert, F.; Llorach, R.; Andrés-Lacueva, C.; Perera, A. Bioinformatics 2014, 30 (13), 1937–1939.
#' @export
svaanno <-
        function(raw,
                 lv,
                 polarity = "positive",
                 projectname = "test") {
                if (is.null(raw$dataCorrected)) {
                        table <- MAIT::MAITbuilder(
                                data = raw$data,
                                spectraID = NULL,
                                masses = raw$mz,
                                rt = raw$rt,
                                spectraEstimation = TRUE,
                                significantFeatures = T,
                                classes = lv,
                                rtRange = 0.2,
                                corThresh = 0.7
                        )
                } else{
                        table <- MAIT::MAITbuilder(
                                data = raw$dataCorrected,
                                spectraID = NULL,
                                masses = raw$mz,
                                rt = raw$rt,
                                spectraEstimation = TRUE,
                                significantFeatures = T,
                                classes = lv,
                                rtRange = 0.2,
                                corThresh = 0.7
                        )
                }
                if (polarity == 'positive') {
                        importMAIT <- MAIT::Biotransformations(
                                MAIT.object = table,
                                adductAnnotation = TRUE,
                                peakPrecision = 0.005,
                                adductTable = NULL
                        )
                } else{
                        importMAIT <- MAIT::Biotransformations(
                                MAIT.object = table,
                                adductAnnotation = TRUE,
                                peakPrecision = 0.005,
                                adductTable = "negAdducts"
                        )
                }

                importMAIT <- MAIT::identifyMetabolites(
                        MAIT.object = importMAIT,
                        peakTolerance = 0.005,
                        polarity = polarity,
                        projectname = projectname
                )
                return(importMAIT)
        }
#' Annotation with xMSannotator package
#' @param raw a xcmsset object to be annotated
#' @param outloc the path for your result
#' @param mode defalut is for positive mode, otherwise use 'neg' for negative
#' @param db_name default is 'HMDB', other database options: 'KEGG', 'LipidMaps', 'T3DB'
#' @param ... parameters for multilevelannotation function in xMSannotator
#' @return as shown in xMSannotator package
#' @references Uppal, K.; Walker, D. I.; Jones, D. P. Anal. Chem. 2017, 89 (2), 1063–1067.
#' @export
svafanno <- function(raw,
                     outloc = "./result/",
                     mode = 'pos',
                     db_name = 'HMDB',num_nodes = 10,...) {
        adduct_weights = cbind.data.frame(Adduct = c('M+H','M-H'),Weight = c(5,5))
        if (is.null(raw$dataCorrected)) {
                data <- raw$data
        }
        else{
                data <- raw$dataCorrected
        }
        mz <- raw$mz
        time <- raw$rt
        data <- as.data.frame(cbind(mz, time, data))
        if( mode == 'neg' ){
                annotres <-
                        xMSannotator::multilevelannotation(
                                dataA = data,
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
