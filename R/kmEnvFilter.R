

#' Faster Kmeans-based environmental subsampling of species records
#' 
#' The function uses a sort of stratified random sampling with strata formed 
#' by k-means clustering and based on environmental variables. 
#' 
#' @param coords A two-column matrix or dataframe with coordinates for 
#' the species records (either longitude/laitude or x/y).
#' 
#' @param envSpace A matrix, dataframe or RasterStack with environmental 
#' variables used to create the strata used in subsampling. If a matrix or 
#' dataframe then each row must match the coordinates in \code{coords}. 
#' If \code{envSpace} is a \code{RasterStack} object then the function 
#' \code{raster::extract} is used to get environmental data.
#' The environmental variables, related to species niche/distributions 
#' should be as meaningful and uncorrelated as possible.
#' 
#' @param k Number of strata/clusters to obtain using K-means (default: 30). 
#' 
#' @param clustFrac The fraction of observations to randomly subsample from each 
#' strata/cluster. Values are within ]0,1[ (default: 0.05).
#' 
#' @param doPCA Do Principal Components Analysis? (default: FALSE). This will 
#' use \code{\link[stats:prcomp]{stats::prcomp}}.
#' 
#' @param centerPCA Center data before PCA? (default: TRUE).
#' 
#' @param scalePCA Scale data before PCA (default: TRUE). 
#' 
#' @param nPC An integer number of Principal Components to extract 
#' from PCA (default: NULL). If defined this parameter will override \code{percVar}.
#' 
#' @param percVar Minimum percentage amount of variance to retain 
#' in Principal Components (default: 90).
#' 
#' @param doPlot Do plot of selected data points? (default: TRUE).
#' 
#' @param ... Additional parameters passed to \code{link[stats:kmeans]{kmeans}}.
#' 
#' @return A set of subsampled records as a dataframe object containing 
#' the following columns:
#' 
#' \itemize{
#'    \item ID - the record unique ID
#'    \item coords - two cloulmns with lon/lat or x/y coordinates
#'    \item clust - the strata/cluster assigned to the observation
#'    \item ... - a set of columns with environmental data used in k-means
#' }
#' 
#' @importFrom raster extract
#' @importFrom stats prcomp
#' @importFrom stats kmeans
#' @importFrom dplyr filter
#' @importFrom dplyr sample_frac
#' @importFrom dplyr bind_rows
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom maps map
#' @importFrom graphics par 
#' @importFrom graphics points
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' 
#' @export
#' 
#' 

kmEnvFilter <- function(coords, envSpace, k = 30, clustFrac = 0.05, 
                        doPCA = FALSE, centerPCA = TRUE, scalePCA = TRUE, 
                        nPC = NULL, percVar = 90,
                        doPlot = TRUE, ...){
  
  if(nrow(coords) != nrow(envSpace)){
    stop("Different number of rows between coords and envSpace!")
  }
  if(clustFrac<=0 | clustFrac>=1){
    stop("clustFrac must be within: ]0, 1[")
  }
  
  if(inherits(envSpace,"RasterStack")){
    rstStack <- envSpace
    envSpace <- raster::extract(rstStack, coords)
  }
  
  if(doPCA){
    pcaObj <- try(stats::prcomp(envSpace, 
                                center = centerPCA, 
                                scale. = scalePCA))
    
    if(inherits(pcaObj,"try-error")){
      stop("An error occurred while trying to perform PCA transformation.")
    }
    
    if(is.null(nPC)){
      expVar <- round((pcaObj$sdev / sum(pcaObj$sdev))*100, 2)
      expVarCumSum <- cumsum(expVar)
      # Number of PCs to extract that have the defined amount of variance
      nPC <- min((1:length(expVarCumSum))[expVarCumSum >= percVar])
    }

    # Retrieve the PCA space for the defined number of PCs or the amount 
    # of variance set
    envSpace_original <- envSpace
    envSpace <- pcaObj$x[,1:nPC]
  }
  
  # Run kmeans
  km <- stats::kmeans(envSpace, centers = k, ...)
  
  # Build dataset with k-means clusters
  selDF <- data.frame(ID = 1:nrow(coords), 
                      coords, 
                      clust = km$cluster, 
                      envSpace)
  
  for(i in 1:k){
    
    tmpSelRows <- selDF %>% 
      dplyr::filter(.data$clust==i) %>% 
      dplyr::sample_frac(size = clustFrac)
    
    if(i==1){
      final_points <- tmpSelRows
    }else{
      final_points <- dplyr::bind_rows(final_points, tmpSelRows)
    }
  }
  
  if(doPlot){
    
    rndCols <- randomcoloR::distinctColorPalette(k)
    
    graphics::par(mfrow=c(1,2), mar=c(4,4,0,0.5))
    
    plot (envSpace[,1], envSpace[,2], pch=19, 
          col=rndCols[final_points$clust], 
          xlab="Dim. 1", ylab="Dim. 2")
    
    graphics::points(final_points[,5], final_points[,6], 
            pch=19, col="black")
    
    plot(coords, pch=19, col="grey50")
    maps::map(add=T)
    graphics::points(final_points[,2:3], pch=19, col="#88000090")
    
  }
  return(final_points)
}

