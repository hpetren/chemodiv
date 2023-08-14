#' Check data formatting
#'
#' Function to check that the datasets used by other functions in
#' the *chemodiv* package are correctly formatted.
#'
#' The function performs a number of checks on the two main datasets used
#' as input data, to make sure datasets are formatted in a way suitable
#' for the other functions in the package. This should make it easier for
#' users to correctly construct datasets before starting with analyses.
#'
#' Two datasets are needed to use the full set of analyses included in
#' the package, and these can be checked for formatting issues.
#' The first dataset should contain data on the proportions
#' of different compounds (columns) in different samples (rows).
#' Note that all calculations of diversity, and most calculations of
#' dissimilarity, are only performed on relative, rather than absolute,
#' values. The second dataset should contain, in each of three
#' columns in a data frame, the compound name, SMILES and InChIKey IDs of
#' all the compounds present in the first dataset. See
#' \code{\link{chemodiv}} for details on obtaining SMILES and InChIKey IDs.
#' Avoid including Greek letters or other special characters in the
#' compound names.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param compoundData Data frame with the compounds in \code{sampleData}
#' as rows. Should have a column named "compound" with common names of
#' the compounds, a column named "smiles" with SMILES IDs of the compounds,
#' and a column named "inchikey" with the InChIKey IDs for the compounds.
#'
#' @return One or several messages pointing out problems with data formatting,
#' or a message informing that the datasets appear to be correctly formatted.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompData)
#' chemoDivCheck(minimalSampData, minimalCompData) # Correct format
#' chemoDivCheck(minimalSampData, minimalCompData[c(2,3,1),]) # Incorrect format
#'
#' data(alpinaSampData)
#' data(alpinaCompData)
#' chemoDivCheck(sampleData = alpinaSampData, compoundData = alpinaCompData)
chemoDivCheck <- function(sampleData,
                          compoundData) {

  formatProblem <- FALSE

  colnames(compoundData) <- tolower(colnames(compoundData))

  if(!is.data.frame(compoundData)) {
    stop("compoundData should be data frame.")
  }
  if(!is.data.frame(sampleData)) {
    stop("sampleData should be data frame.")
  }
  if (!all(c("compound", "smiles", "inchikey") %in% colnames(compoundData))) {
    message("compoundData should include columns with the folowing column names:
            compound, smiles, inchikey.")
    formatProblem <- TRUE
  }
  if (round(sum(rowSums(sampleData))) != nrow(sampleData)) {
    message("Not all row sums in sampleData are equal to 1. Dataset may
            consist of absolute rather than relative values.")
    formatProblem <- TRUE
  }
  if (!all(colnames(sampleData) == compoundData$compound)) {
    message("The name and order of the compounds in the compound column in
            compoundData should be identical to the name and order of the
            columns in sampleData.")
    formatProblem <- TRUE
  }
  if (formatProblem == FALSE) {
    message("The two datasets appear to be correctly formatted.")
  }
}




