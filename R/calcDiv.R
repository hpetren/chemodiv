#' Calculate various diversity/evenness measures
#'
#' Function to calculate Hill diveristy, Functional Hill diversity,
#' Shannon's diveristy, Simpson diveristy, Pielou's evenness, Hill evenness
#' and Rao's Q.
#'
#' Add details here
#' Hill diversity:
#' Functional Hill diversity:
#' Shannon's diversity:
#' Simpson diversity:
#' Pielou's evenness:
#' Hill evenness:
#' Rao's Q:
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param compDisMat Compound distance matrix. Has to be supplied for
#' calculations of Functional Hill diversity or Rao's Q.
#' @param type Type of diversity or evenness to calculate. Any of
#' \code{"HillDiv", "FuncHillDiv", "Shannon", "Simpson", "PielouEven",
#' "HillEven", "RaoQ"}.
#' @param q Diversity order to use for (Functional) Hill diversity
#'
#' @return Data frame with calculated diversity/evenness for each sample
#'
#' @export
#'
#' @examples
#' \dontrun{Add minimal executable example when PubChem PUG works}
calcDiv <- function(sampleData,
                    compDisMat = NULL,
                    type = "HillDiv",
                    q = 1) {

  if(is.null(compDisMat) & ("FuncHillDiv" %in% type | "RaoQ" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity or Rao's Q")
  }

  divData <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(sampleData),
                                  ncol = length(type)))
  colnames(divData) <- type


  if ("HillDiv" %in% type) {

    # This function works row-wise without loop
    divData$HillDiv <- hillR::hill_taxa(comm = sampleData, q = q)

  }

  if ("FuncHillDiv" %in% type) {
    # Probably move this function into a separate (hidden) one when
    # making this into R-package. See section 6.5 in R-package book.
    # Or put it separately with in the same script as the main
    # functions that use this. Section 7.1. Or in utils.R. Same section.

    # My new function to calculate functional Hill diversity.
    # See "FUNCTION DEVELOPMENT" in v.4 for comments on code.
    # (I don't have to make this into a separate function, but I wrote
    # it that way, so keeping it like this now)
    funcHillDiv <- function(data, Dij, q) {

      # First we calculate Rao's Q
      pAbs <- data[data!=0]
      p <- pAbs/sum(pAbs)
      dij <- Dij[data!=0, data!=0]
      Q = sum(dij*(p %*% t(p)))

      temp <- 0


      if(length(p) > 1) { # If we have more than one compound

        # FD for q!=1
        if (q != 1) {
          for (i in 1:length(p)) {
            for (j in 1:length(p)) {
              temp <- temp + dij[i,j]*(p[i]*p[j]/Q)^q
            }
          }
          FD <- temp^(1/(1-q))

        } else { # FD for q==1
          for (i in 1:length(p)) {
            for (j in 1:length(p)) {
              temp <- temp + dij[i,j] * p[i]*p[j]/Q * log(p[i]*p[j]/Q)
            }
          }
          FD <- exp(-temp)
        }

      } else {FD <- NA} # FD is mathematically not defined for one compound
      # so setting NA as that for now, but maybe one could define it to
      # be = 1, which would match normal Hill diversity

      return(FD)
    }

    # My functions needs a loop to work on dataframe
    for (i in 1:nrow(sampleData)) {

      divData$FuncHillDiv[i] <- funcHillDiv(data = sampleData[i,],
                                            Dij = compDisMat,
                                            q = q)
    }

  }

  if ("Shannon" %in% type) {

    divData$Shannon <- vegan::diversity(sampleData, index = "shannon", base = exp(1))

  }

  if ("Simpson" %in% type) { # Note that vegan uses invsimpson for this

    divData$Simpson <- vegan::diversity(sampleData, index = "invsimpson", base = exp(1))

  }

  if ("PielouEven" %in% type) {

    divData$PielouEven <- vegan::diversity(sampleData, index = "shannon", base = exp(1)) /
      log(vegan::specnumber(sampleData))

  }

  if ("HillEven" %in% type) {

    tempHillDiv <- hillR::hill_taxa(comm = sampleData, q = q)
    divData$HillEven <- tempHillDiv / vegan::specnumber(sampleData)

  }

  if ("RaoQ" %in% type) {

    # Function for calculating Raos Q is from DivAnalysesSimulated.R.
    # Not 100% sure how it needs to be, so might have to edit this
    calculateQ = function(data, Dij) {

      Xi <- data[data!=0]
      distance <- Dij[data!=0, data!=0]
      a <- Xi/sum(Xi)
      Q = sum(distance*(a %*% t(a)))

      return(Q)
    }

    for (i in 1:nrow(sampleData)) {

      divData$RaoQ[i] <- calculateQ(data = sampleData[i,],
                                    Dij = compDisMat)

    }
  }

  return(divData)

}
