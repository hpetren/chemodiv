#' Calculate various diversity/evenness measures
#'
#' Function to calculate different common measures of diversity and evenness.
#' This includes Hill diversity, Functional Hill diversity,
#' Shannon's diversity, Simpson diversity, Rao's Q, Pielou's evenness and
#' Hill evenness.
#'
#' The function can calculate the following values:
#' \itemize{
#' \item \code{Shannon}. Shannon's Diversity.
#' \item \code{Simpson}. Simpson Diversity, often referred to as
#' the Inverse Simpson Index.
#' \item \code{HillDiv}. Hill Diversity. Also referred to as the Hill number
#' or the effective number of species (here compounds). The parameter q
#' determines the sensitivity of the measure to the relative frequencies of
#' compounds (see above for details). For q = 0, this equals the number of
#' compounds in a sample. For q = 1, this equals the exponential of
#' Shannon's Diversity. For q = 2, this equals the (inverse) Simpson Diversity.
#' \item \code{FuncHillDiv}. Functional Hill Diversity. Requires a compound
#' dissimilarity matrix. Functional Hill Diversity quantifies the effective
#' total dissimilarity between compounds in the sample. The parameter q
#' determines the  sensitivity of the measure to the relative frequencies
#' of compounds (see above for details). For q = 0, this is
#' equal to Functional Attribute Diversity (FAD) which is the sum of all
#' dissimilarities in the dissimilarity matrix.
#' \item \code{RaoQ}. Rao's Q. The perhaps most common measure of functional
#' diversity. Requires a compound dissimilarity matrix. This represents the
#' average dissimilarity of two randomly selected (weighed by their proportions)
#' compounds in the sample.
#' \item \code{PielouEven}. Pielou's Evenness, also referred to as
#' Shannon's equitability. This is perhaps the most common evenness
#' measure. Equal to the Shannon's Diversity divided by the natural
#' logarithm of the number of compounds. In other words, it expresses evenness
#' with the observed Shannon's diversity as a proportion of the maximum
#' Shannon's diversity where all compounds are equally abundant. Therefore,
#' this is a relative measure with a minimum value of 0 and a maximum value
#' of 1. This measure of evenness is not replication invariant.
#' \item \code{HillEven}. Hill Evenness, as defined by equation 8 in
#' Tuomisto 2012 . This is equal to the Hill Diversity
#' (at the given value of q) divided by the number of compounds, and therefore
#' has a minimum value of 1 / number of compounds and maximum value of 1.
#' This measure of evenness is replication invariant. This measure can be
#' normalized to range from 0 to 1 (eq. 13 in Tuomisto 2012).
#' }
#'
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}. Has to be supplied for calculations of
#' Functional Hill Diversity and Rao's Q.
#' @param type Type(s) of diversity or evenness to calculate. Any of
#' \code{"HillDiv", "FuncHillDiv", "Shannon", "Simpson", "PielouEven",
#' "HillEven", "RaoQ"}.
#' @param q Diversity order to use for Hill diversity, Functional Hill
#' Diversity and Hill Evenness. q should be equal to or larger than zero.
#' This parameter determines the sensitivity of the (Functional) Hill Diversity
#' measure to the relative frequencies of compounds. Commonly set to 0, 1 or 2,
#' although any value > 0 may be used. For q = 0 compound proportions are not
#' taken into account. For q = 1 (default) compounds are weighed according to
#' their proportion in the sample. For q = 2, more weight is put on compounds
#' with high proportions.
#'
#' @return Data frame with calculated diversity/evenness values for each sample.
#'
#' @export
#'
#' @references Chao's and other papers, e.g. Botta-Dukat 2005, Hill 1973,
#' Pielous paper, Tuomisto oikos, Chiu & Chao, Shannon/Simpson papers?
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
              "PielouEven", "HillEven", "RaoQ") %in% type))) {
    stop("Provide at least one type of diversity/evenness to calculate:
         HillDiv, FuncHillDiv, Shannon, Simpson, PielouEven, HillEven or RaoQ.")
  }
  if(is.null(compDisMat) && ("FuncHillDiv" %in% type || "RaoQ" %in% type)) {
    stop("A compound dissimilarity matrix must be supplied
         when calculating Functional Hill diversity or Rao's Q.")
  }
  if(q < 0) {
    stop("q must be >= 0")
  }
  if (any(c("HillDiv", "Shannon", "Simpson",
            "PielouEven", "HillEven") %in% type) &&
      !any(c("FuncHillDiv", "RaoQ") %in% type) &&
      !is.null(compDisMat)) {
    message("Note that the calculated diveristy measures do not use the
            compound dissimilarity matrix.")
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
  return(divData)
}

