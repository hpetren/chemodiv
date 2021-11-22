#' Calculate sample dissimilarities
#'
#' Function to calculate dissimilarities between samples.
#' If only sample data is supplied Bray-Curtis dissimilarities are calculated.
#' If also a compound dissimilarity matrix is supplied, Generalized UniFrac
#' dissimilarities are calculated.
#'
#' @param sampleData Data frame with samples as rows and compounds as columns.
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. If this is supplied UniFrac dissimilarities are
#' calculated, otherwise Bray-Curtis distances are calculated
#' @param alpha Alpha value used to calculate UniFracs. alpha can be set
#' between 0 and 1. With alpha = 0 equal emphasis is put on every branch,
#' With values closer to 1, more emphasis is put on high abundance branches.
#' alpha = 0.5 strikes a balance. alpha 0.5 or 1 is recommended
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

  if(!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &
        all(colnames(sampleData) == rownames(compDisMat)))){
      stop("The sampleData column names, compDisMat column names and
           compDisMat row names must all be identical.")
    }
  }

  if(!is.null(compDisMat)) {

    message("Calculating UniFrac dissimilarities")

    # Note method is now "average" (UPGMA). This is used in Chemoecology paper
    disClust <- stats::hclust(stats::as.dist(compDisMat), method = "average")
    disClustPhylo <- ape::as.phylo(disClust)

    # GUniFrac() returns a list with an array. Seem unnecessarily complicated,
    # so syntax becomes a bit strange
    # Suppressing the warnings about NaNs produced for d_VAW

    uniFracs <- suppressWarnings(GUniFrac::GUniFrac(otu.tab = sampleData,
                                                    tree = disClustPhylo,
                                                    alpha = alpha))
    sampleDisUniFrac <- uniFracs$unifracs[, , paste0("d_", alpha)]
    colnames(sampleDisUniFrac) <- 1:ncol(sampleDisUniFrac)
    rownames(sampleDisUniFrac) <- 1:nrow(sampleDisUniFrac)
    return(sampleDisUniFrac)

  } else {

    message("Calculating Bray-Curtis dissimilarities")
    # If data is not proportional a message is produced. This only
    # applies if we're working with Bray-Curtis, UniFracs are always on props
    if(!sum(rowSums(sampleData)) == nrow(sampleData)) {
      message("sampleData does not contain proportion data.
              This will affect NMDS on Bray-Curtis dissimilarities")
    }
    sampleDisBrayCurtis <- as.matrix(vegan::vegdist(sampleData,
                                                    method = "bray"))
    return(sampleDisBrayCurtis)
  }
}
