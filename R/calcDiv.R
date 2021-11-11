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
#' @param compDisMat Compound dissimilarity matrix. Has to be supplied for
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
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcDiv(sampleData = minimalSampData)
#' calcDiv(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = "FuncHillDiv", q = 2)
#'
calcDiv <- function(sampleData,
                    compDisMat = NULL,
                    type = "HillDiv",
                    q = 1) {

  if(length(type) > 1) {
    stop("Provide only one type of diveristy/evenness to calculate")
  }
  if (!(any(c("HillDiv", "FuncHillDiv", "Shannon", "Simpson",
              "PielouEven", "HillEven", "RaoQ") %in% type))) {
    stop("Provide one type of diversity/evenness to calculate:
         HillDiv, FuncHillDiv, Shannon, Simpson, PielouEven, HillEven or RaoQ")
  }
  if(is.null(compDisMat) & ("FuncHillDiv" %in% type | "RaoQ" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity or Rao's Q")
  }



  divData <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(sampleData),
                                  ncol = length(type)))
  colnames(divData) <- type


  if ("HillDiv" %in% type) {
    if(q < 0) stop("q must be > 0")

    # This function works row-wise without loop
    divData$HillDiv <- hillR::hill_taxa(comm = sampleData, q = q)

  }

  if ("FuncHillDiv" %in% type) {

    if(q < 0) stop("q must be >= 0")

    # My functions needs a loop to work on dataframe.
    for (i in 1:nrow(sampleData)) {
      # Function in utils.R. The chemdiv:: is not needed, in fact this
      # causes all kinds of problems with check(). No idea how R knows
      # how to get funcHillDiv though (but it works)...
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


    for (i in 1:nrow(sampleData)) {
      # function in utils.R
      divData$RaoQ[i] <- calculateQ(data = sampleData[i,],
                                    Dij = compDisMat)

    }
  }

  return(divData)

}

