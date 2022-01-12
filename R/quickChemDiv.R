#' Quickly calculate and plot chemodiversity
#'
#' This function is a shortcut that makes use of many of the other functions in
#' the package. In one simple step chemodiversity is visualized for
#' users wanting to quickly explore their data using standard parameters.
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
#' chemodiversity plots, including (1) a plot for compound dissimilarities,
#' (2) a plot of Functional Hill Diversity at q = 1, (3) a Functional Hill
#' Diversity profile for q = 0-3, and (4) a NMDS plot of the Generalized
#' UniFrac dissimilarities.
#' }
#'
#' @param compoundData Data frame with the chemical compounds of interest.
#' Should have a column named "compound" with common names, a column named
#' "smiles" with SMILES IDs for the compounds, and a column named "inchikey"
#' with the InChIKey IDs for the compounds. See \code{\link{chemdiv}} for
#' details on obtaining SMILES and InChIKey IDs.
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param groupData Grouping data.
#'
#' @return Four plots visualizing chemodiversity.
#'
#' @export
#'
#' @examples
#' data(minimalCompData)
#' data(minimalSampData)
#' groups <- c("A", "A", "B", "B")
#' quickChemDiv(compoundData = minimalCompData, sampleData = minimalSampData,
#' groupData = groups)
quickChemDiv <- function(compoundData,
                         sampleData,
                         groupData = NULL) {

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

  chemDivPlot(compDisMat = quickCompDis$fingerDisMat,
              divData = quickDiv,
              divProfData = quickDivProf,
              sampleDisMat = quickSampleDis$GenUniFrac,
              groupData = groupData)
}
