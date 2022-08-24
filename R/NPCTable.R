#' Generate NPClassifier classification
#'
#' Function to classify compounds with *NPClassifier*, and put the results in
#' a data frame containing the pathway, superclass and class for each compound.
#'
#' *NPClassifier* (Kim et al. 2021) is a deep-learning tool that automatically
#' classifies natural products (i.e. phytochemical compounds) into a
#' hierarchical classification of three levels: pathway, superclass and class.
#' This classification largely corresponds to the biosynthetic groups/pathways
#' the compounds are produced in. The \code{NPCTable} function
#' conveniently performs this classification directly in R on the
#' compounds in \code{compoundData}, by accessing the tool
#' at \url{https://npclassifier.ucsd.edu/} and downloading the classifications.
#'
#' @param compoundData Data frame with the chemical compounds of interest,
#' usually the compounds found in the sample dataset.
#' Should include a column named "compound" with common names of
#' the compounds and a column named "smiles" with SMILES IDs of the compounds.
#'
#' @return Data frame with the NPClassifier classification for each compound
#' as pathway, superclass and class. Note that compounds may be classified
#' in more than one group, or no group, at each level of classification.
#'
#' @export
#'
#' @references
#' Kim HW, Wang M, Leber CA, Nothias L-F, Reher R, Kang KB,
#' van der Hooft JJJ, Dorrestein PC, Gerwick WH, Cottrell GW. 2021.
#' NPClassifier: A Deep Neural Network-Based Structural Classification
#' Tool for Natural Products. Journal of Natural Products 84: 2795-2807.
#'
#' @examples
#' data(minimalCompData)
#' NPCTable(minimalCompData)
#'
#' data(alpinaCompData)
#' NPCTable(compoundData = alpinaCompData[1:3,]) # First three compounds only
NPCTable <- function(compoundData) {

  httr::set_config(httr::config(http_version = 0))
  if(!curl::has_internet()) {
    message("No internet connection available to download compound data.")
    return(invisible(NULL))
  } else if (httr::GET("https://npclassifier.ucsd.edu/")$status_code != 200) {
    message("NPClassifier, https://npclassifier.ucsd.edu/, is unavailable.")
    return(invisible(NULL))
  }

  colnames(compoundData) <- tolower(colnames(compoundData))
  npcTab <- compoundData

  # Assuming no compound has >3 belongings pathways/superclasses/classes
  npcTab[c("pathway", "pathway2", "pathway3",
           "superclass", "superclass2", "superclass3",
           "class", "class2", "class3")] <- NA

  # Get NPC for each SMILE
  for (i in 1:nrow(npcTab)) {

    if (!is.na(npcTab$smiles[i])) {

      # (May produce strange error if run separately)
      npcClass1 <- httr::GET("https://npclassifier.ucsd.edu/classify",
                             query = list(smiles = npcTab$smiles[i]))

      npcClass2 <- httr::content(npcClass1, as = "text")

      # If request fails, npcClass2 does not begin with "{", check for that
      if (substring(npcClass2, 1, 1) == "{") {

        npcClass3 <- jsonlite::fromJSON(npcClass2)

        # Check if there are any classifications, then put data in data frame
        if (is.character(npcClass3$pathway_results)) {

          npcTab$pathway[i] <- npcClass3$pathway_results[1]
          npcTab$pathway2[i] <- npcClass3$pathway_results[2]
          npcTab$pathway3[i] <- npcClass3$pathway_results[3]

        } else {
          message(paste("NPClassifier has no classification for compound", i))
        }

        if (is.character(npcClass3$superclass_results)) {
          npcTab$superclass[i] <- npcClass3$superclass_results[1]
          npcTab$superclass2[i] <- npcClass3$superclass_results[2]
          npcTab$superclass3[i] <- npcClass3$superclass_results[3]
        }

        if (is.character(npcClass3$class_results)) {
          npcTab$class[i] <- npcClass3$class_results[1]
          npcTab$class2[i] <- npcClass3$class_results[2]
          npcTab$class3[i] <- npcClass3$class_results[3]
        }
        # If the output from NPClassifier API is not as expected
      } else {
        message(paste0("NPClassifier gave error output for compound ", i,". ",
                      "Is the SMILES correct?"))
      }
    }
  }
  if(all(is.na(npcTab$pathway))) {
    stop("No compounds could be classified.")
  }
  # Removes columns with only NA
  npcTab <- npcTab[, colSums(is.na(npcTab)) < nrow(npcTab)]
  return(npcTab)
}
