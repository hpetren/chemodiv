#' Calculate a diversity profile
#'
#' Function to calculate a diversity profile, i.e. calculate Hill diversity
#' or Functional Hill Diversity for a range of q values.
#'
#' A diversity profile is a calculation of Hill Diversity or Functional Hill
#' Diversity at different values of q. This function performs the calculations,
#' while \code{\link{chemDivPlot}} can be used to conveniently create the
#' diversity profile plot, where Hill Diversity is plotted as a function of q
#' within the chosen range. The shape of the diversity profile curve
#' reflects  the unevenness of compound proportions in the sample. For a
#' perfectly even sample the curve is flat. The more uneven the compound
#' proportions are, the more steep is the decline of the curve.
#' A common range, used as default, of q values is between \code{qMin = 0} and
#' \code{qMax = 3} as diversity should change little beyond \code{qMax = 3}.
#' See \code{\link{calcDiv}} for further details on q.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compDisMat Compound distance matrix, as calculated by
#' \code{\link{compDis}}. Has to be supplied for calculations of
#' Functional Hill diversity.
#' @param type Type of Hill Diversity to calculate for the diversity profile.
#' \code{"HillDiv"} or \code{"FuncHillDiv"}.
#' @param qMin Minimum value of q.
#' @param qMax Maximum value of q.
#' @param step Increment by which q will be calculated between qMin and qMax.
#'
#' @return List with a diversity profile data frame with samples as rows
#' and the Hill diversity or Functional Hill diversity for different q values
#' as columns; and values for type, qMin, qMax and step.
#'
#' @export
#'
#' @references Who came up with divprofs? But e.g. Jost 2010,
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcDivProf(sampleData = minimalSampData)
#' calcDivProf(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = "FuncHillDiv", qMin = 1, qMax = 2, step = 0.2)
calcDivProf <- function(sampleData,
                        compDisMat = NULL,
                        type = "HillDiv",
                        qMin = 0,
                        qMax = 3,
                        step = 0.1) {

  if (!(type == "HillDiv" || type == "FuncHillDiv") || length(type) > 1) {
    stop("Provide a single type of diversity profile to calculate:
         HillDiv or FuncHillDiv.")
  }
  if(is.null(compDisMat) && ("FuncHillDiv" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity.")
  }
  if (!is.null(compDisMat) && type == "HillDiv") {
    message("Note that the calculated diveristy profile does not use the
            compound dissimilarity matrix.")
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

  # Creating datasets to store diversity values
  qAll <- seq(from = qMin, to = qMax, by = step)
  divProf <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(sampleData),
                                  ncol = length(qAll)))
  colnames(divProf) <- paste0("q", qAll)

  for (c in 1:length(qAll)) { # For each column
    for (r in 1:nrow(divProf)) { # For each row
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

