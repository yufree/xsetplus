% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch.R
\name{limmapca}
\alias{limmapca}
\title{Principal component analysis(PCA) for limma result with batch}
\usage{
limmapca(list, center = T, scale = T, lv = NULL)
}
\arguments{
\item{list}{results from `limmafit` function or related batch correction method}

\item{center}{parameters for PCA}

\item{scale}{parameters for scale}

\item{lv}{factor vector for the group infomation}
}
\value{
plot
}
\description{
Principal component analysis(PCA) for limma result with batch
}
\examples{
\dontrun{
library(faahKO)
library(enviGCMS)
cdfpath <- system.file('cdf', package = 'faahKO')
list <- getmr(cdfpath, pmethod = ' ')
li <- svacor(list$data,list$group$class)
limmapca(li)
}
}
\seealso{
\code{\link{isvacor}},\code{\link{svacor}}, \code{\link{pcacor}},\code{\link{limmafit}}
}
