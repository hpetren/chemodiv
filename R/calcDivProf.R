#' Calculate diversity profile
#'
#' Function to calculate a "diversity profile", i.e. calculate (Functional)
#' Hill diveristy for multiple values of q.
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param compDisMat Compound distance matrix. Has to be supplied for
#' calculations of Functional Hill diversity or Rao's Q.
#' @param qMin Minimume value of q
#' @param qMax Maximum value of q
#' @param step Step by which q will be calculated between qMin and qMax
#'
#' @return List with a diversity profile data frame with samples as rows
#' and the (Functional) Hill diversity at different q values as columns;
#' and values for qMin, qMax and step.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcDivProf(sampleData = minimalSampData)
#' calcDivProf(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' qMin = 1, qMax = 2, step = 0.2)
calcDivProf <- function(sampleData,
                        compDisMat = NULL,
                        qMin = 0,
                        qMax = 3,
                        step = 0.1) {

  if(qMin < 0) stop("qMin should be >= 0")

  if(qMin > qMax) stop("qMin should be smaller than qMax")

  if((step <= 0) | (step > (qMax - qMin)))
    stop("step must be > 0 and less than qMax - qMin")



  # Creating datasets to store diversity values
  qAll <- seq(from = qMin, to = qMax, by = step)

  divHillProf <- as.data.frame(matrix(data = NA,
                                      nrow = nrow(sampleData),
                                      ncol = length(qAll)))

  colnames(divHillProf) <- paste0("q", qAll)


  for (c in 1:length(qAll)) {

    # Calculate diversity for each row (each sample) in that column
    for (r in 1:nrow(divHillProf)) {

      if (!is.null(compDisMat)) {
        divHillProf[r,c] <- funcHillDiv(data = sampleData[r,], Dij = compDisMat, q = qAll[c])
      } else {
        divHillProf[r,c] <- hillR::hill_taxa(sampleData[r,], q = qAll[c])
      }

    }
  }

  diversityProfile <- list("divHillProf" = divHillProf, "qMin" = qMin, "qMax" = qMax, "step" = step)

  return(diversityProfile)

}

