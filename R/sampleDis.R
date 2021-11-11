#' Calculate sample dissimilarities
#'
#' Function to calculate dissimilarities between samples.
#' If only sample data is supplied Bray-Curtis dissimilarities are calculated.
#' If also a compound dissimilarity matrix is supplied, Generalized UniFrac
#' dissimilarities are calculated.
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param compDisMat Compound dissimilarity matrix. If this is supplied UniFrac
#' dissimilarities are calcualted, otherwise Bray-Curtis distances
#' are calculated
#' @param alpha Alpha value used to calculate UniFracs
#'
#' @return Sample dissimilarity matrix with Bray-Curtis or
#' UniFrac dissimilarities
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' sampleDis(minimalSampData)
#' sampleDis(minimalSampData, minimalCompDis, alpha = 0.5)
sampleDis <- function(sampleData,
                      compDisMat = NULL,
                      alpha = 1) {

  # If no compound dissimilarity matrix is supplied, and data is not
  # proportional a warning is produced. (If compound dissimilarity matrix
  # is supplied, UniFracs are calculated, which are always based on
  # proportions, so no warning then)
  if(is.null(compDisMat) & !sum(rowSums(sampleData)) == nrow(sampleData)) {
    warning("sampleData does not contain proportion data.
            This will affect NMDS on Bray-Curtis dissimilarities")
  }
  if(!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &
        all(colnames(sampleData) == rownames(compDisMat)))){
      stop("The sampleData column names, compDisMat column names and
           compDisMat row names must all be identical.")
    }
  }

  if(!is.null(compDisMat)) {

    print("Calculating UniFrac dissimilarities")

    # Note method is now "average" (UPGMA). This is used in Chemoecology paper
    disClust <- stats::hclust(stats::as.dist(compDisMat), method="average")

    # Does this require ape-package all of a sudden?
    disClustPhylo <- ape::as.phylo(disClust)

    # GUniFrac() returns a list with an array. Seem unnecessarily complicated,
    # so syntax becomes a bit strange (but I think it's correct)
    # Doesn't seem to warn about the NaNs produced for d_VAW anymore.

    uniFracs <- GUniFrac::GUniFrac(sampleData, disClustPhylo, alpha = alpha)

    sampleDisUniFrac <- uniFracs$unifracs[, , paste("d_",alpha,sep="")]

    colnames(sampleDisUniFrac) <- 1:ncol(sampleDisUniFrac)
    rownames(sampleDisUniFrac) <- 1:nrow(sampleDisUniFrac)

    return(sampleDisUniFrac)

  } else {

    print("Calculating Bray-Curtis dissimilarities")

    sampleDisBrayCurtis <- as.matrix(vegan::vegdist(sampleData, method = "bray"))

    return(sampleDisBrayCurtis)

  }
}
