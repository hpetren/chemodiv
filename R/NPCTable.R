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
#' compounds in \code{compoundData}, by accessing the
#' tool at \url{https://npclassifier.ucsd.edu/}.
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
#' \dontrun{
#' data(alpinaCompData)
#' NPCTable(compoundData = alpinaCompData)
#' }
NPCTable <- function(compoundData) {

  colnames(compoundData) <- tolower(colnames(compoundData))
  npcTab <- compoundData

  # Assuming for now that a compound never has more than 3 belongings
  # in a pathway/superclass/class
  npcTab[c("pathway", "pathway2", "pathway3",
           "superclass", "superclass2", "superclass3",
           "class", "class2", "class3")] <- NA

  for (i in 1:nrow(npcTab)) {

    # No NA
    if (!is.na(npcTab$smiles[i])) {

      # Get the NPC for the SMILES
      # NOTE! This gives an error, but still works...
      npcClass1 <- httr::GET("https://npclassifier.ucsd.edu/classify",
                             query = list(smiles = npcTab$smiles[i]))

      # Retrieving contents
      npcClass2 <- httr::content(npcClass1, as = "text")

      # If the request for some reason fails then npcClass2
      # is not in proper format, and does not begin with "{",
      # so we check for that (we have to do that, because otherwise
      # fromJSON() throws an error and loop stops)
      if (substring(npcClass2, 1, 1) == "{") {

        # Parsing the long json character string, into a list
        npcClass3 <- jsonlite::fromJSON(npcClass2)

        # Check if there is any classifications (otherwise empty list), then
        # put in correct place in df
        if (is.character(npcClass3$pathway_results)) {

          # Adding pathway. In cases there is only 1 pathway, trying to
          # take an index that is higher then what is in the vector
          # just returns NA, which is appropriate
          npcTab$pathway[i] <- npcClass3$pathway_results[1]
          npcTab$pathway2[i] <- npcClass3$pathway_results[2]
          npcTab$pathway3[i] <- npcClass3$pathway_results[3]

          # If the compound is not classified, but there is nothing wrong
          # (i.e. the pathway is NA)
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
  # If no compounds were classified
  if(all(is.na(npcTab$pathway))) {
    stop("No compounds could be classified.")
  }
  # Removes columns that are only NA
  npcTab <- npcTab[, colSums(is.na(npcTab)) < nrow(npcTab)]
  return(npcTab)
}
