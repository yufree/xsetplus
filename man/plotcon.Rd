% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temp.R
\name{plotcon}
\alias{plotcon}
\title{Plot circle connection of metabolites from list}
\usage{
plotcon(list, class = NULL, cutoff = 0.6)
}
\arguments{
\item{list}{list with data as peaks list, mz, rt and group information}

\item{class}{metabolites or compounds class, default NULL}

\item{cutoff}{cutoff of the absolute value of correlation coefficient}
}
\description{
Plot circle connection of metabolites from list
}
\examples{
\dontrun{
library(faahKO)
cdfpath <- system.file('cdf', package = 'faahKO')
list <- getmr(cdfpath, pmethod = ' ')
data <- list$data
lv <- as.character(list$group$class)
plotcon(list, lv)
}
}
