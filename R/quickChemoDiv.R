#' Quickly calculate or plot chemodiversity
#'
#' This function is a shortcut that makes use of many of the other functions
#' in the package. In one simple step chemodiversity is calculated, and if
#' requested also plotted, for users wanting to quickly explore their
#' data using standard parameters.
#'
#' The function requires sample data as input, and can also include
#' compound data. \code{\link{chemoDivCheck}} can be used to ensure these
#' datasets are correctly formatted, see \code{\link{chemodiv}} for further
#' details on data formatting. If only sample data is supplied,
#' phytochemical diversity and dissimilarity will be calculated
#' as Hill diversity and Bray-Curtis dissimilarity, respectively.
#' If sample data and compound data is supplied, phytochemical diversity
#' and dissimilarity will be calculated as Functional Hill diversity
#' and Generalized UniFrac dissimilarity, respectively.
#' This function then uses the following other functions in the package:
#' \itemize{
#' \item \code{\link{compDis}} is used to calculate compound dissimilarities
#' using PubChem Fingerprints, if compound data is supplied.
#' \item \code{\link{calcDiv}} is used to calculate (Functional) Hill
#' Diversity for q = 1.
#' \item \code{\link{calcDivProf}} is used to calculate a diversity profile
#' with (Functional) Hill Diversity for *q* = 0-3.
#' \item \code{\link{sampDis}} is used to calculate Bray-Curtis or
#' Generalized UniFrac dissimilarities between samples.
#' \item \code{\link{chemoDivPlot}} is used to create different
#' chemodiversity plots if requested.
#' }
#'
#' \code{quickChemoDiv} is designed to provide an easy way to calculate and
#' visualize the most central measures of phytochemical diversity. It uses
#' default parameters to do so, which should be reasonable in most cases.
#' However, for detailed analyses it is recommended to use the separate
#' functions to allow for full control of function input, arguments and output.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compoundData Data frame with the compounds in \code{sampleData}
#' as rows. Should have a column named "compound" with common names of
#' the compounds, a column named "smiles" with SMILES IDs of the compounds,
#' and a column named "inchikey" with the InChIKey IDs for the compounds.
#' See \code{\link{chemodiv}} for details on obtaining SMILES and InChIKey IDs.
#' @param groupData Grouping data (e.g. population, species etc.).
#' Should be either a vector, or a data frame with a single column.
#' @param outputType Type of output that should be returned: either
#' \code{data} to output a list with different types of chemodiversity data,
#' or \code{plots} to instead produce standard plots of this data.
#'
#' @return Different types of chemodiversity measures, either as elements in
#' a list or as separate plots. If \code{outputType = "data"}, function returns
#' a compound dissimilarity matrix (if compound data was supplied),
#' a data frame with (Functional) Hill Diversity at *q* = 1,
#' a data frame with a (Functional) Hill Diversity profile for *q* = 0-3,
#' and a sample dissimilarity matrix. If \code{outputType = "plots"},
#' these data sets are plotted as a dendrogram (if compound data was supplied),
#' a boxplot, a diversity profile plot and an NMDS plot, respectively.
#'
#' @export
#'
#' @examples
#' data(minimalCompData)
#' data(minimalSampData)
#' groups <- c("A", "A", "B", "B")
#' quickChemoDiv(sampleData = minimalSampData, groupData = groups,
#' outputType = "data") # Without compound data
#'
#' data(alpinaSampData)
#' data(alpinaPopData)
#' quickChemoDiv(sampleData = alpinaSampData, outputType = "plots",
#' groupData = alpinaPopData) # Without compound data
quickChemoDiv <- function(sampleData,
                          compoundData = NULL,
                          groupData = NULL,
                          outputType = "plots") {

  if (length(outputType) > 1) {
    stop('outputType must be either data or plots.')
  }
  if (!(outputType == "plots" || outputType == "data")) {
    stop('outputType must be either data or plots.')
  }

  if (is.null(compoundData)) {

    quickDiv <- calcDiv(sampleData = sampleData,
                        type = "HillDiv",
                        q = 1)

    quickDivProf <- calcDivProf(sampleData = sampleData,
                                type = "HillDiv",
                                qMin = 0,
                                qMax = 3,
                                step = 0.1)

    quickSampDis <- sampDis(sampleData = sampleData,
                            type = "BrayCurtis")

    if (outputType == "data") {
      chemoDivData <- list(HillDiv = quickDiv,
                           HillDivProf = quickDivProf$divProf,
                           SampleDissimilarity = quickSampDis$BrayCurtis)
      return(chemoDivData)
    } else {
      chemoDivPlot(divData = quickDiv,
                   divProfData = quickDivProf,
                   sampDisMat = quickSampDis$BrayCurtis,
                   groupData = groupData)
    }

  } else {
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

    quickSampDis <- sampDis(sampleData = sampleData,
                            compDisMat = quickCompDis$fingerDisMat,
                            type = "GenUniFrac",
                            alpha = 1)

    if (outputType == "data") {
      chemoDivData <- list(CompoundDissimilarity = quickCompDis$fingerDisMat,
                           FuncHillDiv = quickDiv,
                           FuncHillDivProf = quickDivProf$divProf,
                           SampleDissimilarity = quickSampDis$GenUniFrac)
      return(chemoDivData)
    } else {
      chemoDivPlot(compDisMat = quickCompDis$fingerDisMat,
                   divData = quickDiv,
                   divProfData = quickDivProf,
                   sampDisMat = quickSampDis$GenUniFrac,
                   groupData = groupData)
    }
  }
}
