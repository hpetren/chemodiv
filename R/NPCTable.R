#' NPClassifier Data Table
#'
#' This function classifies compounds with NPClassifier,
#' and puts the results in a dataframe containing the pathway,
#' superclass and class for each compound.
#'
#' Connects to \url{https://npclassifier.ucsd.edu/}
#'
#' @usage NPCTable(compoundData)
#'
#' @param compoundData Data frame with the chemical compounds of interest.
#' Should have a column named "compound" with common names, and a column named
#' "smiles" with SMILES IDs for the compounds.
#'
#' @return Data frame with the NPClassifier classification for each compound
#' as pathways, superclass and class. Compounds may be classified in more than
#' one group, or no group, at each level of classification.
#'
#' @export
#'
#' @references
#' Kim, H. W., M. Wang, C. A. Leber, L.-F. Nothias, R. Reher, K. B. Kang,
#' J. J. J. van der Hooft, P. C. Dorrestein, W. H. Gerwick,
#' and G. W. Cottrell. 2021. NPClassifier: A Deep Neural Network-Based
#' Structural Classification Tool for Natural Products.
#' Journal of Natural Products:acs.jnatprod.1c00399.
#' UPDATE WITH FULL REF WHEN AVAILABLE
#'
#' @examples
#' data(minimalCompData)
#' NPCTable(minimalCompData)
NPCTable <- function(compoundData) {


  npcTable <- compoundData

  # Assuming for now that a compound never has more than 3 belongings
  # in a certain pathway/superclass/class (automating this further
  # down seems tricky, so sticking with this for now)
  npcTable$pathway <- NA
  npcTable$pathway2 <- NA
  npcTable$pathway3 <- NA
  npcTable$superclass <- NA
  npcTable$superclass2 <- NA
  npcTable$superclass3 <- NA
  npcTable$class <- NA
  npcTable$class2 <- NA
  npcTable$class3 <- NA


  for (i in 1:nrow(npcTable)) {

    # Avoiding NA
    if(!is.na(npcTable$smiles[i])) {

      # Get the NPC for the SMILES and get into correct format
      npcclass_prel <- httr::GET("https://npclassifier.ucsd.edu/classify",
                                query = list(smiles = npcTable$smiles[i]))


      # This is just one long character string
      npcclass_real <- httr::content(npcclass_prel, as = "text")



      # If the request for some reason fails (which it seems to do for
      # some compounds, although I don't know why since the SMILES work
      # when doing it manually on the webpage) then npcclass_real
      # is not in proper json format, and does not begin with "{",
      # so we check for that (we have to do that, because otherwise
      # fromJSON() throws an error and loop stops)
      if(substring(npcclass_real, 1, 1) == "{") {

        #print(paste0("Classifying compound: ", i))

        # Parsing the long character string, which is formatted as json,
        # into a list
        npcclass_real_correct <- jsonlite::fromJSON(npcclass_real)


        # Putting in correct place in df with putting in alt if there
        # are 2 records for a level (and an extra check which I think
        # not is unnecessary)
        # Note order
        if (is.character(npcclass_real_correct$pathway_results)) {

          # Adding pathway. In cases there is only 1 pathway, trying to
          # take in index that is higher then what is in the vector
          # just returns NA, which is appropriate (maybe not the fastest
          # solution, but it works and I don't have to construct more loops)
          npcTable$pathway[i] <- npcclass_real_correct$pathway_results[1]
          npcTable$pathway2[i] <- npcclass_real_correct$pathway_results[2]
          npcTable$pathway3[i] <- npcclass_real_correct$pathway_results[3]

          # If the compound is not classified, but there is nothing wrong
          # (i.e. the pathway is NA)
        } else { print(paste0("NPClassifier has no classification for compound ", i)) }


        if (is.character(npcclass_real_correct$superclass_results)) {

          npcTable$superclass[i] <- npcclass_real_correct$superclass_results[1]
          npcTable$superclass2[i] <- npcclass_real_correct$superclass_results[2]
          npcTable$superclass3[i] <- npcclass_real_correct$superclass_results[3]

        }


        if (is.character(npcclass_real_correct$class_results)) {

          npcTable$class[i] <- npcclass_real_correct$class_results[1]
          npcTable$class2[i] <- npcclass_real_correct$class_results[2]
          npcTable$class3[i] <- npcclass_real_correct$class_results[3]

        }

        # If something went wrong
      } else { warning(paste0("Compound ", i, " could not be classified at all. Is the SMILES correct?")) }
    }

  }

  # Removes columns that are only NA
  # (This is purely aesthetic, and does not affect other functions anyways)
  npcTable <- npcTable[,colSums(is.na(npcTable)) < nrow(npcTable)]

  return(npcTable)

}
