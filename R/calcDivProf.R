#' Calculate a diversity profile
#'
#' Function to calculate a diversity profile, i.e. calculate Hill diversity
#' or Functional Hill Diversity for a range of *q* values.
#'
#' The function calculates a diversity profile for each sample
#' in \code{sampleData}. A diversity profile is a calculation of
#' Hill Diversity or Functional Hill Diversity for a range of
#' different values of *q*. This function performs the calculations,
#' while \code{\link{chemoDivPlot}} can be used to conveniently
#' create the diversity profile plot, where Hill Diversity is
#' plotted as a function of *q* within the chosen range.
#' The shape of the diversity profile curve reflects the evenness
#' of compound proportions in the sample. For a perfectly even sample
#' the curve is flat. The more uneven the compound proportions are,
#' the more steep is the decline of the curve. A common range,
#' used as default, of *q* values is between \code{qMin = 0} and
#' \code{qMax = 3}, as diversity should change little beyond \code{qMax = 3}.
#' See \code{\link{calcDiv}} for further details on *q*.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compDisMat Compound distance matrix, as calculated by
#' \code{\link{compDis}}. Has to be supplied for calculations of
#' Functional Hill diversity.
#' @param type Type of Hill Diversity to calculate for the diversity profile.
#' \code{"HillDiv"} or \code{"FuncHillDiv"}.
#' @param qMin Minimum value of *q*.
#' @param qMax Maximum value of *q*.
#' @param step Increment by which *q* will be calculated between \code{qMin}
#' and \code{qMax}.
#'
#' @return List with a diversity profile data frame with samples as rows
#' and the Hill diversity or Functional Hill diversity for different *q*
#' values as columns; and values for \code{type}, \code{qMin},
#' \code{qMax} and \code{step}.
#'
#' @export
#'
#' @references
#' Chao A, Chiu C-H, Jost L. 2014. Unifying Species Diversity,
#' Phylogenetic Diversity, Functional Diversity, and Related Similarity and
#' Differentiation Measures Through Hill Numbers.
#' Annual Review of Ecology, Evolution, and Systematics 45: 297-324.
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcDivProf(sampleData = minimalSampData)
#' calcDivProf(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = "FuncHillDiv")
#'
#' data(alpinaCompData)
#' data(alpinaCompDis)
#' calcDivProf(sampleData = alpinaSampData, compDisMat = alpinaCompDis,
#' type = "FuncHillDiv", qMin = 1, qMax = 2, step = 0.2)
calcDivProf <- function(sampleData,
                        compDisMat = NULL,
                        type = "HillDiv",
                        qMin = 0,
                        qMax = 3,
                        step = 0.1) {

  if (length(type) > 1) {
    stop("Provide only one type of diversity profile to calculate.")
  }
  if (!(type == "HillDiv" || type == "FuncHillDiv")) {
    stop("type should be HillDiv or FuncHillDiv.")
  }
  if(is.null(compDisMat) && ("FuncHillDiv" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity.")
  }
  if (!is.null(compDisMat) && type == "HillDiv") {
    message("Note that the calculated diveristy profile does not use the
            compound dissimilarity matrix.")
  }
  if (!is.null(compDisMat)) {
    if(!(all(colnames(sampleData) == colnames(compDisMat)) &&
         all(colnames(sampleData) == rownames(compDisMat)))) {
      stop("The name and order of the columns in sampleData should be identical
         to the name and order of the columns/rows in compDisMat.")
    }
  }
  if(qMin < 0) {
    stop("qMin should be >= 0.")
  }
  if(qMin > qMax) {
    stop("qMin should be less than qMax.")
  }
  if((step <= 0) || (step > (qMax - qMin))) {
    stop("step must be > 0 and less than qMax - qMin.")
  }

  qAll <- seq(from = qMin, to = qMax, by = step)
  divProf <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(sampleData),
                                  ncol = length(qAll)))
  colnames(divProf) <- paste0("q", qAll)

  for (c in 1:length(qAll)) { # For each column (q)
    for (r in 1:nrow(divProf)) { # For each row (sample)
      if (type == "FuncHillDiv") {
        divProf[r, c] <- funcHillDiv(data = sampleData[r,],
                                     Dij = compDisMat,
                                     q = qAll[c])
      } else {
        divProf[r, c] <- hillR::hill_taxa(sampleData[r,],
                                          q = qAll[c])
      }
    }
  }

  diversityProfile <- list("divProf" = divProf,
                           "type" = type,
                           "qMin" = qMin,
                           "qMax" = qMax,
                           "step" = step)
  return(diversityProfile)
}

