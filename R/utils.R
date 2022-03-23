#' Functional Hill diversity
#'
#' Function to calculate Functional Hill diversity, as in equations
#' 4b and 6b in Chiu & Chao 2014 Plos One.
#' Note that this is done for one sample (row).
#'
#' @param data Dataset. Can be proportional or absolute data.
#' @param Dij Dissimilarity matrix.
#' @param q Diversity order.
#'
#' @return Functional Hill diversity (FD).
funcHillDiv <- function(data, Dij, q) {

  # Rao's Q
  pAbs <- data[data != 0]
  p <- pAbs / sum(pAbs)
  dij <- Dij[data != 0, data != 0]
  Q <- sum(dij * (p %*% t(p)))

  temp <- 0

  if(length(p) > 1) { # If more than one compound

    # FD for q != 1
    if (q != 1) {
      for (i in 1:length(p)) {
        for (j in 1:length(p)) {
          temp <- temp + dij[i, j] * (p[i] * p[j] / Q)^q
        }
      }
      FD <- temp^(1 / (1 - q))

    } else { # FD for q == 1
      for (i in 1:length(p)) {
        for (j in 1:length(p)) {
          temp <- temp + dij[i, j] * p[i] * p[j] / Q * log(p[i] * p[j] / Q)
        }
      }
      FD <- exp(-temp)
    }

  } else {
    # FD is mathematically not defined with one compound, so setting to NA
    FD <- NA
  }
  return(FD)
}

#' Rao's Q
#'
#' Function to calculate Rao's Q.
#' Note that this is done for one sample (row).
#'
#' @param data Dataset. Can be proportional or absolute data.
#' @param Dij Dissimilarity matrix.
#'
#' @return Value of Rao's Q.
calculateQ = function(data, Dij) {

  pAbs <- data[data != 0]
  p <- pAbs / sum(pAbs)
  dij <- Dij[data != 0, data != 0]
  Q <- sum(dij * (p %*% t(p)))

  return(Q)
}
