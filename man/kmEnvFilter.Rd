% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmEnvFilter.R
\name{kmEnvFilter}
\alias{kmEnvFilter}
\title{Faster Kmeans-based environmental subsampling of species records}
\usage{
kmEnvFilter(
  coords,
  envSpace,
  k = 30,
  clustFrac = 0.05,
  doPCA = FALSE,
  centerPCA = TRUE,
  scalePCA = TRUE,
  nPC = NULL,
  percVar = 90,
  doPlot = TRUE,
  ...
)
}
\arguments{
\item{coords}{A two-column matrix or dataframe with coordinates for 
the species records (either longitude/laitude or x/y).}

\item{envSpace}{A matrix, dataframe or RasterStack with environmental 
variables used to create the strata used in subsampling. If a matrix or 
dataframe then each row must match the coordinates in \code{coords}. 
If \code{envSpace} is a \code{RasterStack} object then the function 
\code{raster::extract} is used to get environmental data.
The environmental variables, related to species niche/distributions 
should be as meaningful and uncorrelated as possible.}

\item{k}{Number of strata/clusters to obtain using K-means (default: 30).}

\item{clustFrac}{The fraction of observations to randomly subsample from each 
strata/cluster. Values are within ]0,1[ (default: 0.05).}

\item{doPCA}{Do Principal Components Analysis? (default: FALSE). This will 
use \code{\link[stats:prcomp]{stats::prcomp}}.}

\item{centerPCA}{Center data before PCA? (default: TRUE).}

\item{scalePCA}{Scale data before PCA (default: TRUE).}

\item{nPC}{An integer number of Principal Components to extract 
from PCA (default: NULL). If defined this parameter will override \code{percVar}.}

\item{percVar}{Minimum percentage amount of variance to retain 
in Principal Components (default: 90).}

\item{doPlot}{Do plot of selected data points? (default: TRUE).}

\item{...}{Additional parameters passed to \code{link[stats:kmeans]{kmeans}}.}
}
\value{
A set of subsampled records as a dataframe object containing 
the following columns:

\itemize{
   \item ID - the record unique ID
   \item coords - two cloulmns with lon/lat or x/y coordinates
   \item clust - the strata/cluster assigned to the observation
   \item ... - a set of columns with environmental data used in k-means
}
}
\description{
The function uses a sort of stratified random sampling with strata formed 
by k-means clustering and based on environmental variables.
}
