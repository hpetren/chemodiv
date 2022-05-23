#' Calculate sample dissimilarities
#'
#' Function to calculate dissimilarities between samples.
#' Either Bray-Curtis dissimilarities and/or Generalized UniFrac
#' dissimilarities are calculated.
#'
#' The function calculates a dissimilarity matrix for all the samples
#' in \code{sampleData}, for the given dissimilarity index/indices.
#' Bray-Curtis dissimilarities are calculated using only
#' the \code{sampleData}. This is the most commonly calculated dissimilarity
#' index used for phytochemical data (other types of dissimilarities are
#' easily calculated using the \code{\link[vegan]{vegdist}} function in
#' the \code{vegan} package).
#'
#' If a compound dissimilarity matrix, \code{compDisMat}, is supplied,
#' Generalized UniFrac dissimilarities can be calculated, which also
#' use the compound dissimilarity matrix for the sample dissimilarity
#' calculations. For the calculation of Generalized UniFrac
#' dissimilarities (Chen et al. 2012), the compound dissimilarity matrix is
#' transformed into a dendrogram using hierarchical clustering (with the
#' UPGMA method). Calculations of UniFrac dissimilarities quantifies the
#' fraction of the total branch length of the dendrogram that leads to
#' compounds present in either sample, but not both. The (weighted) Generalized
#' UniFrac dissimilarities implemented here additionally take compound
#' abundances into account. In this way, both the relative proportions of
#' compounds and the biosynthetic/structural dissimilarities of the compounds
#' are accounted for in the calculations of sample dissimilarities, such that
#' two samples containing more biosynthetically/structurally different
#' compounds have a higher pairwise dissimilarity than two samples
#' containing more biosynthetically/structurally similar compounds.
#' As with Bray-Curtis dissimilarities, Generalized UniFrac dissimilarities
#' range in value from 0 to 1.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. If this is supplied, Generalized UniFrac
#' dissimilarities can be calculated.
#' @param type Type of sample dissimilarities to be calculated. This is
#' Bray-Curtis dissimilarities, \code{type = "BrayCurtis"}, and/or
#' Generalized UniFrac dissimilarities, \code{type = "GenUniFrac"}.
#' @param alpha Parameter used in calculations of Generalized UniFrac
#' dissimilarities. alpha can be set between 0 and 1.
#' With \code{alpha = 0}, equal weight is put on every
#' branch in the dendrogram. With \code{alpha = 1}, branches are
#' weighted by their abundance, and hence more emphasis is put on high
#' abundance branches. \code{alpha = 0.5} strikes a balance between the two.
#' alpha 0.5 or 1 is recommended, with \code{alpha = 1} as default.
#' See Chen et al. 2012 for details.
#'
#' @return List with sample dissimilarity matrices. A list is always
#' outputted, even if only one matrix is calculated.
#'
#' @references
#' Bray JR, Curtis JT. 1957. An Ordination of the Upland Forest Communities
#' of Southern Wisconsin. Ecological Monographs 27: 325-349.
#'
#' Chen J, Bittinger K, Charlson ES, et al. 2012. Associating microbiome
#' composition with environmental covariates using generalized
#' UniFrac distances. Bioinformatics 28: 2106-2113.
#'
#' Lozupone C, Knight R. 2005. UniFrac: a New Phylogenetic Method for Comparing
#' Microbial Communities. Applied and Environmental Microbiology 71: 8228-8235.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' sampDis(minimalSampData)
#' sampDis(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = c("BrayCurtis", "GenUniFrac"), alpha = 0.5)
#'
#' data(alpinaSampData)
#' data(alpinaCompDis)
#' sampDis(sampleData = alpinaSampData, compDisMat = alpinaCompDis,
#' type = "GenUniFrac")
sampDis <- function(sampleData,
                      compDisMat = NULL,
                      type = "BrayCurtis",
                      alpha = 1) {

  if (!(any(c("BrayCurtis", "GenUniFrac") %in% type))) {
    stop("Provide at least one type of dissimilarity to calculate:
         BrayCurtis or GenUniFrac.")
  }
  if (("GenUniFrac" %in% type) && is.null(compDisMat)) {
    stop("A compound dissimilarity matrix must be supplied when calculating
         Generalized UniFrac dissimilarities.")
  }
  if (!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &&
        all(colnames(sampleData) == rownames(compDisMat)))) {
      stop("The name and order of the columns in sampleData should be identical
         to the name and order of the columns/rows in compDisMat.")
    }
  }
  if (("BrayCurtis" %in% type) && (length(type) == 1) && !is.null(compDisMat)) {
    message("Note that Bray-Curtis dissimilarity calculations do not use the
            compound dissimilarity matrix.")
  }

  sampDisMatList <- list()

  if ("BrayCurtis" %in% type) { # Bray-Curtis

    # Turn non-proportion data into proportion-data
    if(!sum(rowSums(sampleData)) == nrow(sampleData)) {
      sampleData <- sampleData / rowSums(sampleData)

      message("sampleData appears to not contain proportion data. Data is made
              into proportions before dissimiliarty calcualtions, as Bray-Curtis
              dissimilarities calculated on absolute and proportional data
              are not identical.")
    }
    sampDisBrayCurtis <- as.matrix(vegan::vegdist(sampleData,
                                                    method = "bray"))
    sampDisMatList[["BrayCurtis"]] <- sampDisBrayCurtis
  }

  if ("GenUniFrac" %in% type) { # Generalized UniFracs

    # Note method "average" (UPGMA)
    disClust <- stats::hclust(stats::as.dist(compDisMat), method = "average")
    disClustPhylo <- ape::as.phylo(disClust)

    # Suppressing warnings about NaNs produced for d_VAW (not relevant)
    uniFracs <- suppressWarnings(GUniFrac::GUniFrac(otu.tab = sampleData,
                                                    tree = disClustPhylo,
                                                    alpha = alpha))
    sampDisUniFrac <- uniFracs$unifracs[, , paste0("d_", alpha)]
    colnames(sampDisUniFrac) <- 1:ncol(sampDisUniFrac)
    rownames(sampDisUniFrac) <- 1:nrow(sampDisUniFrac)
    sampDisMatList[["GenUniFrac"]] <- sampDisUniFrac
  }
  return(sampDisMatList)
}
