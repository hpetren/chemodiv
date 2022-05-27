#' Functional Hill diversity
#'
#' Function to calculate Functional Hill diversity, as in equations
#' 4b and 6b in Chiu & Chao 2014 PLoS One.
#' Note that this is done for one sample (row).
#'
#' @param data Dataset. Single row in dataframe. Proportional or absolute data.
#' @param Dij Dissimilarity matrix.
#' @param q Diversity order.
#'
#' @noRd
#'
#' @return Functional Hill diversity.
funcHillDiv <- function(data, Dij, q) {

  # Rao's Q
  pAbs <- data[data != 0]
  p <- pAbs / sum(pAbs)
  dij <- Dij[data != 0, data != 0]
  Q <- sum(dij * (p %*% t(p)))

  temp <- 0

  if(length(p) > 1) { # If more than one compound

    # FHD for q != 1
    if (q != 1) {
      for (i in 1:length(p)) {
        for (j in 1:length(p)) {
          temp <- temp + dij[i, j] * (p[i] * p[j] / Q)^q
        }
      }
      FD <- temp^(1 / (1 - q))

    } else { # FHD for q == 1
      for (i in 1:length(p)) {
        for (j in 1:length(p)) {
          temp <- temp + dij[i, j] * p[i] * p[j] / Q * log(p[i] * p[j] / Q)
        }
      }
      FD <- exp(-temp)
    }

  } else {
    # FHD not mathematically defined with one compound
    FD <- NA
  }
  return(FD)
}

#' Rao's Q
#'
#' Function to calculate Rao's Q.
#' Note that this is done for one sample (row).
#'
#' @param data Dataset. Single row in dataframe. Proportional or absolute data.
#' @param Dij Dissimilarity matrix.
#'
#' @noRd
#'
#' @return Value of Rao's Q.
calculateQ <- function(data, Dij) {

  pAbs <- data[data != 0]
  p <- pAbs / sum(pAbs)
  dij <- Dij[data != 0, data != 0]
  Q <- sum(dij * (p %*% t(p)))

  return(Q)
}
