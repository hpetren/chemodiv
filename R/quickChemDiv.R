#' Quickly calculate/plot chemodiversity
#'
#' This function is a shortcut that makes use of many of the other functions in
#' the package. In one simple step chemodiversity is calculated, and if
#' requested also visualized, for users wanting to quickly explore their
#' data using standard parameters.
#'
#' The function takes compound data and sample data as input. See
#' \code{\link{chemdiv}} for details on data formatting. This function then
#' uses the following other functions in the package:
#' \itemize{
#' \item \code{\link{compDis}} is used to calculate compound dissimilarities
#' using PubChem Fingerprints.
#' \item \code{\link{calcDiv}} is used to calculate Functional Hill Diversity
#' for q = 1.
#' \item \code{\link{calcDivProf}} is used to calculate a diversity profile
#' with Functional Hill Diversity for q = 0-3.
#' \item \code{\link{sampleDis}} is used to calculate Generalized UniFrac
#' dissimilarities between samples.
#' \item \code{\link{chemDivPlot}} is used to create four different
#' chemodiversity plots if requested.
#' }
#'
#' \code{quickChemDiv} is designed to provide an easy way to visualize
#' the most central measures of phytochemical diversity. It uses
#' default parameters to do so, which should be reasonable in most cases.
#' However, for detailed analyses we recommend using the separate
#' functions to allow for full control of function parameters and output.
#'
#' @param compoundData Data frame with the chemical compounds of interest.
#' Should have a column named "compound" with common names, a column named
#' "smiles" with SMILES IDs for the compounds, and a column named "inchikey"
#' with the InChIKey IDs for the compounds. See \code{\link{chemdiv}} for
#' details on obtaining SMILES and InChIKey IDs.
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param groupData Grouping data.
#' @param outputType Type of output that should be returned: either
#' \code{"data"} to output a list with different types of chemodiversity data,
#' or \code{"plots"} to instead produce standard plots of this data.
#'
#' @return Four types of chemodiversity measures, as either elements in a list
#' or separate plots. If \code{outputType = "data"}, function returns
#' a compound dissimilarity matrix, a data frame with Functional Hill
#' Diversity at q = 1, a data frame with a Functional Hill Diversity profile
#' for q = 0-3, and a sample dissimilarity matrix. If
#' \code{outputType = "plots"}, these data sets are plotted as a dendrogram,
#' a boxplot, a diversity profile plot and an NMDS plot, respectively.
#'
#' @export
#'
#' @examples
#' data(minimalCompData)
#' data(minimalSampData)
#' groups <- c("A", "A", "B", "B")
#' quickChemDiv(compoundData = minimalCompData, sampleData = minimalSampData,
#' groupData = groups, outputType = "plots")
quickChemDiv <- function(compoundData,
                         sampleData,
                         groupData = NULL,
                         outputType = "plots") {

  if (length(outputType) > 1 ||
      !(outputType == "plots" || outputType == "data")) {
    stop('outputType must be either "data" or "plots".')
  }

  quickCompDis <- compDis(compoundData = compoundData,
                          type = "PubChemFingerprint",
                          npcTable = NULL,
                          unknownCompoundsMean = FALSE)

  quickDiv <- calcDiv(sampleData = sampleData,
                      compDisMat = quickCompDis$fingerDisMat,
                      type = "FuncHillDiv",
                      q = 1)

  quickDivProf <- calcDivProf(sampleData = sampleData,
                              compDisMat = quickCompDis$fingerDisMat,
                              type = "FuncHillDiv",
                              qMin = 0,
                              qMax = 3,
                              step = 0.1)

  quickSampleDis <- sampleDis(sampleData = sampleData,
                              compDisMat = quickCompDis$fingerDisMat,
                              type = "GenUniFrac",
                              alpha = 1)

  if (outputType == "data") {
    chemDivData <- list(compoundDissimilarity = quickCompDis$fingerDisMat,
                        funcHillDiv = quickDiv,
                        funcHillDivProf = quickDivProf$divProf,
                        sampleDissimilarity = quickSampleDis$GenUniFrac)
  } else {
    chemDivPlot(compDisMat = quickCompDis$fingerDisMat,
                divData = quickDiv,
                divProfData = quickDivProf,
                sampleDisMat = quickSampleDis$GenUniFrac,
                groupData = groupData)
  }
}
