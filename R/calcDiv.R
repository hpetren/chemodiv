#' Calculate various diversity/evenness measures
#'
#' Function to calculate Hill diveristy, Functional Hill diversity,
#' Shannon's diveristy, Simpson diveristy, Pielou's evenness, Hill evenness,
#' Rao's Q and Functional Attribute Diversity (FAD).
#'
#' Add details here
#' Hill diversity:
#' Functional Hill diversity: With compound matrix.
#' Shannon's diversity:
#' Simpson diversity:
#' Pielou's evenness:
#' Hill evenness:
#' Rao's Q: With compound matrix.
#' FAD: Functional Attribute Diversity.
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. Has to be supplied for
#' calculations of Functional Hill diversity, Rao's Q and FAD.
#' @param type Type(s) of diversity or evenness to calculate. Any of
#' \code{"HillDiv", "FuncHillDiv", "Shannon", "Simpson", "PielouEven",
#' "HillEven", "RaoQ", "FAD"}.
#' @param q Diversity order to use for (Functional) Hill diversity. q should
#' be equal to or larger than zero. This parameter determines the sensitivity
#' of the (Functional) Hill diversity measure to the relative frequencies
#' of compounds. For q = 0 compound proportions are not taken into account.
#' For q = 1 compounds are weighed according to their proportion in the sample.
#' For q = 2, more weight is put on compounds with high proportions.
#'
#' @return Data frame with calculated diversity/evenness for each sample.
#'
#' @export
#'
#' @references Chao's and other papers
#' Chao, A., C.-H. Chiu, and L. Jost. 2014. Unifying Species Diversity,
#' Phylogenetic Diversity, Functional Diversity, and Related Similarity
#' and Differentiation Measures Through Hill Numbers. Annual Review of
#' Ecology, Evolution, and Systematics 45:297–324.
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' calcDiv(sampleData = minimalSampData)
#' calcDiv(sampleData = minimalSampData, type = c("HillDiv", "HillEven"))
#' calcDiv(sampleData = minimalSampData, compDisMat = minimalCompDis,
#' type = "FuncHillDiv", q = 2)
calcDiv <- function(sampleData,
                    compDisMat = NULL,
                    type = "HillDiv",
                    q = 1) {

  if (!(any(c("HillDiv", "FuncHillDiv", "Shannon", "Simpson",
              "PielouEven", "HillEven", "RaoQ", "FAD") %in% type))) {
    stop("Provide at least one type of diversity/evenness to calculate:
         HillDiv, FuncHillDiv, Shannon, Simpson, PielouEven, HillEven or RaoQ")
  }
  if(is.null(compDisMat) && ("FuncHillDiv" %in% type ||
                             "RaoQ" %in% type ||
                             "FAD" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity, Rao's Q, or FAD.")
  }
  if(q < 0) {
    stop("q must be >= 0")
  }

  divData <- as.data.frame(matrix(data = NA,
                                  nrow = nrow(sampleData),
                                  ncol = length(type)))
  colnames(divData) <- type

  if ("HillDiv" %in% type) {
    # This function works row-wise without loop
    divData$HillDiv <- hillR::hill_taxa(comm = sampleData,
                                        q = q)
  }
  if ("FuncHillDiv" %in% type) {
    # My utils.R function needs a loop to work on dataframe
    for (i in 1:nrow(sampleData)) {
      divData$FuncHillDiv[i] <- funcHillDiv(data = sampleData[i,],
                                            Dij = compDisMat,
                                            q = q)
    }
  }
  if ("Shannon" %in% type) {
    divData$Shannon <- vegan::diversity(sampleData,
                                        index = "shannon",
                                        base = exp(1))
  }
  if ("Simpson" %in% type) { # Note that vegan uses invsimpson for this
    divData$Simpson <- vegan::diversity(sampleData,
                                        index = "invsimpson",
                                        base = exp(1))
  }
  if ("PielouEven" %in% type) {
    divData$PielouEven <- vegan::diversity(sampleData,
                                           index = "shannon",
                                           base = exp(1)) /
      log(vegan::specnumber(sampleData))
  }
  if ("HillEven" %in% type) {
    tempHillDiv <- hillR::hill_taxa(comm = sampleData,
                                    q = q)
    divData$HillEven <- tempHillDiv / vegan::specnumber(sampleData)
  }
  if ("RaoQ" %in% type) {
    for (i in 1:nrow(sampleData)) {
      # function in utils.R
      divData$RaoQ[i] <- calculateQ(data = sampleData[i,],
                                    Dij = compDisMat)
    }
  }
  if ("FAD" %in% type) {
    divData$FAD <- sum(compDisMat)
  }
  return(divData)
}

