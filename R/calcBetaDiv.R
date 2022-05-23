#' Calculate beta diversity
#'
#' Function to calculate beta diversity in the Hill diversity framework.
#' This can be calculated as Hill beta diversity or
#' Functional Hill beta diversity.
#'
#' The function calculates a single beta diversity value for the supplied
#' \code{sampleData}. This is calculated as *beta = gamma / alpha*. Gamma
#' diversity represents the diversity of the pooled data set, alpha diversity
#' represents the mean diversity across individual samples, and
#' beta diversity represents turnover or variability among samples.
#' With \code{type = "HillDiv"} and \code{q = 0} the calculated beta diversity
#' is equal to the well-known and most simple measure of beta diversity
#' introduced by Whittaker 1960, where *beta = gamma / alpha*, based only
#' on the number of species (here compounds).
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. Has to be supplied for calculations of
#' Functional Hill beta diversity.
#' @param type Type(s) of Hill beta diversity to calculate. \code{"HillDiv"}
#' and/or \code{"FuncHillDiv"}.
#' @param q Diversity order to use for the calculation of beta diversity.
#' See \code{\link{calcDiv}} for further details on *q*.
#'
#' @return Data frame with type of Hill beta diversity calculated, *q*, and
#' values for gamma diversity, mean alpha diversity and beta diversity.
#'
#' @export
#'
#' @references
#' Chao A, Chiu C-H, Jost L. 2014. Unifying Species Diversity,
#' Phylogenetic Diversity, Functional Diversity, and Related Similarity and
#' Differentiation Measures Through Hill Numbers.
#' Annual Review of Ecology, Evolution, and Systematics 45: 297-324.
#'
#' Jost L. 2007. Partitioning diversity into independent alpha and
#' beta components. Ecology 88: 2427-2439.
#'
#' Whittaker RH. 1960. Vegetation of the Siskiyou Mountains, Oregon
#' and California. Ecological Monographs 30: 279-338.
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcBetaDiv(sampleData = minimalSampData)
#' calcBetaDiv(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = c("HillDiv", "FuncHillDiv"), q = 2)
#'
#' data(alpinaSampData)
#' data(alpinaCompDis)
#' calcBetaDiv(sampleData = alpinaSampData, compDisMat = alpinaCompDis,
#' type = "FuncHillDiv")
calcBetaDiv <- function(sampleData,
                        compDisMat = NULL,
                        type = "HillDiv",
                        q = 1) {

  if (!(any(c("HillDiv", "FuncHillDiv") %in% type))) {
    stop("Provide at least one type of beta-diversity to calculate:
         HillDiv or FuncHillDiv.")
  }
  if(is.null(compDisMat) && ("FuncHillDiv" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill beta-diversity.")
  }
  if(q < 0) {
    stop("q must be >= 0")
  }
  if (!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &&
         all(colnames(sampleData) == rownames(compDisMat)))) {
      stop("The name and order of the columns in sampleData should be identical
         to the name and order of the columns/rows in compDisMat.")
    }
  }
  if (!is.null(compDisMat) && ("HillDiv" %in% type) && length(type) == 1) {
    message("Note that the calculated beta-diveristy does not use the
            compound dissimilarity matrix.")
  }

  betaDiv <- as.data.frame(matrix(data = NA,
                                  nrow = length(type),
                                  ncol = 5))
  colnames(betaDiv) <- c("Type", "q", "gamma", "alpha", "beta")
  betaDiv$Type <- type

  if ("HillDiv" %in% type) {
    betaDiv[betaDiv$Type == "HillDiv", 2:5] <-
      hillR::hill_taxa_parti(comm = sampleData, q = q)[, 1:4]
  }
  if ("FuncHillDiv" %in% type) {
    betaDiv[betaDiv$Type == "FuncHillDiv", 2:5] <-
      hillR::hill_func_parti(comm = sampleData,
                             traits = compDisMat,
                             traits_as_is = TRUE,
                             q = q)[, c(1, 3:5)]
  }
  return(betaDiv)
}
