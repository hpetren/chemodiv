#' Calculate sample dissimilarities
#'
#' Function to calculate dissimilarities between samples.
#' This is either Bray-Curtis dissimilarities and/or Generalized UniFrac
#' dissimilarities.
#'
#' \code{sampleDis} calculates a dissimilarity matrix for all the samples
#' in \code{sampleData}. Bray-Curtis dissimilarities are calculated using only
#' the \code{sampleData}. If a compound dissimilarity matrix, \code{compDisMat},
#' is supplied, Generalized UniFrac dissimilarities can also be calculated,
#' which uses the compound dissimilarity matrix (transformed into a dendrogram)
#' in the dissimilarity calculations. Thereby, the biosynthetic/structural
#' properties of the compounds is taken into account for the calculations
#' of sample dissimilarities.
#'
#'
#' @param sampleData Data frame with samples as rows and compounds as columns.
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. If this is supplied, Generalized UniFrac
#' dissimilarities can be calculated.
#' @param type Type of sample dissimilarity to be calculated. This is
#' Bray-Curtis dissimilarity, \code{BrayCurtis} and/or Generalized UniFrac
#' dissimilarity, \code{GenUniFrac}.
#' @param alpha Alpha value used to calculate Generalized UniFracs. alpha can
#' be set between 0 and 1. With alpha = 0 equal emphasis is put on every
#' branch in the dendrogram. With values closer to 1, more emphasis is put
#' on high abundance branches. alpha = 0.5 strikes a balance.
#' alpha 0.5 or 1 is recommended, with alpha = 1 as default.
#'
#' @return List with sample dissimilarity matrices. A list is always
#' outputted, even if only one matrix is calculated.
#'
#' @references Bray-Curtis ref?
#' Chen, J., K. Bittinger, E. S. Charlson, C. Hoffmann, J. Lewis, G. D. Wu,
#' R. G. Collman, F. D. Bushman, and H. Li. 2012.
#' Associating microbiome composition with environmental covariates using
#' generalized UniFrac distances. Bioinformatics 28:2106–2113.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' sampleDis(minimalSampData)
#' sampleDis(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = c("BrayCurtis", "GenUniFrac"), alpha = 0.5)
sampleDis <- function(sampleData,
                      compDisMat = NULL,
                      type = "BrayCurtis",
                      alpha = 1) {

  if (!(any(c("BrayCurtis", "GenUniFrac") %in% type))) {
    stop("Provide at least one type of dissimilarity to calculate:
         BrayCurtis or GenUniFrac")
  }
  if ("GenUniFrac" %in% type && is.null(compDisMat)) {
    stop("A compound dissimilarity matrix must be supplied when calculating
         Generalized UniFrac dissimilarities")
  }
  if (!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &&
        all(colnames(sampleData) == rownames(compDisMat)))) {
      stop("The sampleData column names, compDisMat column names and
           compDisMat row names must all be identical.")
    }
  }
  if ("BrayCurtis" %in% type && length(type) == 1 && !is.null(compDisMat)) {
    message("Note that Bray-Curtis dissimilarity calculations do not use the
            compound dissimilarity matrix.")
  }

  sampleDisMatList <- list()

  if ("BrayCurtis" %in% type) {

    # If data is not proportional a message is produced. This only
    # applies if we're working with Bray-Curtis, UniFracs are always on props
    if(!sum(rowSums(sampleData)) == nrow(sampleData)) {
      message("sampleData does not contain proportion data. Bray-Curtis
              dissimilarities calculated on absolute and proportional data
              are not identical.")
    }
    sampleDisBrayCurtis <- as.matrix(vegan::vegdist(sampleData,
                                                    method = "bray"))
    sampleDisMatList[["BrayCurtis"]] <- sampleDisBrayCurtis
  }

  if ("GenUniFrac" %in% type) {

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
    sampleDisMatList[["GenUniFrac"]] <- sampleDisUniFrac
  }
  return(sampleDisMatList)
}
