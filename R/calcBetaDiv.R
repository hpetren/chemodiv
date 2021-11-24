#' Calculate beta diversity
#'
#' Function to calculate Hill beta-diversity or Functional Hill beta-diversity.
#' Done with hillR:hill_taxa_parti, see that for details. In essence,
#' beta = gamma/alpha, and one value is generated for the whole data set.
#'
#' Hill beta-diversity is calculated if compDisMat is not supplied,
#' Functional Hill beta-diversity is calculated if compDisMat is supplied.
#' The gamma- and alpha-diversity values used to calculate beta are also output.
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param compDisMat Compound distance matrix, as calculated by
#' \code{\link{NPCTable}}. Has to be supplied for
#' calculations of Functional Hill beta-diversity.
#' @param q Diversity order to use for (Functional) Hill diversity.
#'
#' @return Data frame with alpha-, beta- and gamma-diversity.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcBetaDiv(sampleData = minimalSampData)
#' calcDiv(sampleData = minimalSampData, compDisMat = minimalCompDis, q = 0)
#'
calcBetaDiv <- function(sampleData,
                        compDisMat = NULL,
                        q = 1) {
  if(q < 0) {
    stop("q must be >= 0")
  }

  if (is.null(compDisMat)) { # If no compound matrix, normal Hill beta
    betaDiv <- hillR::hill_taxa_parti(comm = sampleData,
                                      q = q)
    betaDivOut <- betaDiv[, 2:4]
    colnames(betaDivOut) <- c("gamma", "alpha", "beta")
    return(betaDivOut)
  } else { # If compound matrix, functional Hill beta
    betaDiv <- hillR::hill_func_parti(comm = sampleData,
                                      traits = compDisMat,
                                      traits_as_is = TRUE,
                                      q = q)
    betaDivOut <- betaDiv[, 3:5]
    colnames(betaDivOut) <- c("gamma", "alpha", "beta")
    return(betaDivOut)
  }
}
