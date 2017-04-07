#' @description annotation with MAIT package
#' @param path the path to your data
#' @param name the name of your project
#' @param ion defalut is for positive mode, otherwise use 'negAdducts' for negative
#' @return as shown in MAIT package
#' @references Fernández-Albert, F.; Llorach, R.; Andrés-Lacueva, C.; Perera, A. Bioinformatics 2014, 30 (13), 1937–1939.
#' @examples
#' path <- "./data/"
#' name <- "fishriver"
#' anno(path,name)
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

#' @description annotation with xMSannotator package
#' @param xset a xcmsset object to be annotated
#' @param outloc the path for your result
#' @param mode defalut is for positive mode, otherwise use 'neg' for negative
#' @param list adductlist for annotation. the default adductlist is queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") . You might use other options for negative mode: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
#' @param db_name default is 'HMDB', other database options: 'KEGG', 'LipidMaps', 'T3DB'
#' @return as shown in xMSannotator package
#' @references Uppal, K.; Walker, D. I.; Jones, D. P. Anal. Chem. 2017, 89 (2), 1063–1067.
#' @examples
#' path <- "./data/"
#' xset <- getdata(path)
#' result <- fanno(xset)
#' @export

fanno <-
        function(xset,
                 outloc = "./result/",
                 mode = 'pos',
                 list = c(
                         "M+2H",
                         "M+H+NH4",
                         "M+ACN+2H",
                         "M+2ACN+2H",
                         "M+H",
                         "M+NH4",
                         "M+Na",
                         "M+ACN+H",
                         "M+ACN+Na",
                         "M+2ACN+H",
                         "2M+H",
                         "2M+Na",
                         "2M+ACN+H",
                         "M+2Na-H",
                         "M+H-H2O",
                         "M+H-2H2O"
                 ),
                 db_name = 'HMDB') {
                adduct_weights = cbind.data.frame(Adduct = c('M+H','M-H'),Weight = c(5,5))
                mz <- xcms::groups(xset)[, 1]
                time <- xcms::groups(xset)[, 4]
                data <- as.data.frame(cbind(mz, time, data))
                data <- unique(data)
                annotres <-
                        xMSannotator::multilevelannotation(
                                dataA = data,
                                max.mz.diff = 5,
                                max.rt.diff = 10,
                                cormethod = "pearson",
                                num_nodes = 12,
                                queryadductlist = list,
                                mode = mode,
                                outloc = outloc,
                                db_name = db_name,
                                adduct_weights = adduct_weights,
                                num_sets = 1000,
                                allsteps = TRUE,
                                corthresh = 0.7,
                                NOPS_check = TRUE,
                                customIDs = NA,
                                missing.value = NA,
                                hclustmethod = "complete",
                                deepsplit = 2,
                                networktype = "unsigned",
                                minclustsize = 10,
                                module.merge.dissimilarity = 0.2,
                                filter.by = c("M+H"),
                                biofluid.location = NA,
                                origin = NA,
                                status = "Detected and Quantified",
                                boostIDs = NA,
                                max_isp = 5,
                                HMDBselect = "union",
                                mass_defect_window = 0.01,
                                pathwaycheckmode = "pm",
                                mass_defect_mode = mode
                        )
                return(annotres)
        }

#' @description output pathway annotation data for Mummichog algorithm
#' @param xset a xcmsset object to be annotated
#' @param lv the group factor
#' @param name the file name
#' @param method parameter for groupval function
#' @param instensity parameter for groupval function
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
                mod <- model.matrix(~ lv)
                mod0 <- as.matrix(c(rep(1, ncol(data))))
                fstats <- sva::fstats(data, mod, mod0)
                pvalue <- sva::f.pvalue(data, mod, mod0)
                df <- cbind.data.frame(mz, rt, pvalue, fstats)
                filename <- paste0(name, '.txt')
                write.table(df,
                            file = filename,
                            sep = "\t",
                            row.names = F)
                return(df)
        }

#' @description annotation with MAIT package of list from svacor function
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
#' @description annotation with xMSannotator package
#' @param xset a xcmsset object to be annotated
#' @param outloc the path for your result
#' @param mode defalut is for positive mode, otherwise use 'neg' for negative
#' @param list adductlist for annotation. the default adductlist is queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O") . You might use other options for negative mode: c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H"); c("positive"); c("negative"); c("all");see data(adduct_table) for complete list
#' @param db_name default is 'HMDB', other database options: 'KEGG', 'LipidMaps', 'T3DB'
#' @return as shown in xMSannotator package
#' @references Uppal, K.; Walker, D. I.; Jones, D. P. Anal. Chem. 2017, 89 (2), 1063–1067.
#' @export
svafanno <- function(raw,
                     outloc = "./result/",
                     mode = 'pos',
                     list = c(
                             "M+2H",
                             "M+H+NH4",
                             "M+ACN+2H",
                             "M+2ACN+2H",
                             "M+H",
                             "M+NH4",
                             "M+Na",
                             "M+ACN+H",
                             "M+ACN+Na",
                             "M+2ACN+H",
                             "2M+H",
                             "2M+Na",
                             "2M+ACN+H",
                             "M+2Na-H",
                             "M+H-H2O",
                             "M+H-2H2O"
                     ),
                     db_name = 'HMDB') {
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
        annotres <-
                xMSannotator::multilevelannotation(
                        dataA = data,
                        max.mz.diff = 5,
                        max.rt.diff = 10,
                        cormethod = "pearson",
                        num_nodes = 12,
                        queryadductlist = list,
                        mode = mode,
                        outloc = outloc,
                        db_name = db_name,
                        adduct_weights = adduct_weights,
                        num_sets = 1000,
                        allsteps = TRUE,
                        corthresh = 0.7,
                        NOPS_check = TRUE,
                        customIDs = NA,
                        missing.value = NA,
                        hclustmethod = "complete",
                        deepsplit = 2,
                        networktype = "unsigned",
                        minclustsize = 10,
                        module.merge.dissimilarity = 0.2,
                        filter.by = c("M+H"),
                        biofluid.location = NA,
                        origin = NA,
                        status = "Detected and Quantified",
                        boostIDs = NA,
                        max_isp = 5,
                        HMDBselect = "union",
                        mass_defect_window = 0.01,
                        pathwaycheckmode = "pm",
                        mass_defect_mode = mode
                )
        return(annotres)
}
