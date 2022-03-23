#' Generate molecular network with properties
#'
#' Function which generates a molecular network object, and some
#' basic properties of the network.
#'
#' Note that the supplied dissimilarity matrix is translated to a
#' similarity matrix, and this is what cutOff values are set for.
#'
#' @param compDisMat Compound dissimilarity matrix, as calculated by
#' \code{\link{compDis}}.
#' @param npcTable An already generated \code{\link{NPCTable}} can be supplied
#' for calculations of the number of NPC pathways and network modularity.
#' @param cutOff Cut-off value for compound similarities. Any similarity
#' lower than this value will be set to zero when the network is generated,
#' which strongly affects the look of the network. The value can be set
#' manually 0 <= cutOff <= 1 (maybe mention that others use 0.6 for another
#' kind of molecular similarity measure), to the \code{median} dissimilarity
#' value in the compDisMat, or to \code{minPathway}, the lowest within-pathway
#' similarity (which allows all within-NPC-pathway to be kept) if an NPCTable
#' is supplied.
#'
#' @return List with a graph (tidygrph/igraph?) object, the number of compounds,
#' number of NPC pathways and modularity.
#'
#' @export
#'
#' @examples
#' data(minimalCompDis)
#' data(minimalNPCTable)
#' molNet(minimalCompDis)
#' molNet(minimalCompDis, minimalNPCTable, cutOff = 0)
#'
#' \dontrun{
#' data(alpinaCompData)
#' alpinaNPCTable <- NPCTable(compoundData = alpinaCompData)
#' alpinaCompDis <- compDis(compoundData = alpinaCompData)
#' molNet(compDisMat = alpinaCompDis$fingerDisMat,
#' npcTable = alpinaNPCTable, cutOff = 0.75)
#' }
molNet <- function(compDisMat,
                   npcTable = NULL,
                   cutOff = "median") {

  if (is.numeric(cutOff) && ((cutOff < 0) || (cutOff > 1))) {
    stop("Numeric values for cutOff must be between 0 and 1.")
  }
  if (is.character(cutOff) && !(cutOff == "median" || cutOff == "minPathway")) {
    stop("cutOff must be a value between 0 and 1, median or minPathway.")
  }
  if (!all(colnames(compDisMat) == rownames(compDisMat))) {
    stop("compDisMat must be a square matrix with identical row and column names.")
  }

  # From dissimilarity to similarity matrix
  compSimMat <- 1 - compDisMat

  # Keep diagnoal being 0
  diag(compSimMat) <- 0

  # Keep the full similarity matrix
  compSimMatFull <- compSimMat

  # Set cut-off point. Median, manual or lowest within-pathway similarity
  if (cutOff == "median") {

    message("Using median similarity as cut-off.")
    medianSim <- stats::median(compSimMat[lower.tri(compSimMat)])
    compSimMat[compSimMat < medianSim] <- 0

  } else if (is.numeric(cutOff)) {

    message(paste("Using cut-off value =", cutOff))
    compSimMat[compSimMat < cutOff] <- 0

  } else if (cutOff == "minPathway" && !is.null(npcTable)) {

    message("Using lowest within-pathway similarity as cut-off.")
    minIntraPath <- max(compSimMat)
    simValue <- NA

    for (col in 1:(ncol(compSimMat)-1)) {
      for (row in (col+1):nrow(compSimMat)) {
        # Lower triangle
        if (!is.na(npcTable$pathway[col]) && !is.na(npcTable$pathway[row]) &&
            npcTable$pathway[col] == npcTable$pathway[row]) {
          # If both are not NA and compounds are from the same pathway
          simValue <- compSimMat[row,col]
          if (simValue < minIntraPath) {
            # If this value is smaller than previous one, replace
            minIntraPath <- simValue
          }
        }
      }
    }
    compSimMat[compSimMat < minIntraPath] <- 0

  } else if (cutOff == "minPathway" && is.null(npcTable)) {
    stop("Using minPathway as cut-off requires npcTable.")
  }

  # Creating network
  networkObject <- tidygraph::as_tbl_graph(compSimMat)

  # Number of compounds
  nCompounds <- nrow(compSimMatFull)

  # Maybe remove this since we added FAD.
  # Mean similarity. This is done on FULL matrix (not cut-off) as
  # that makes more sense (you want the mean similarity of the full bouquet,
  # not an arbitrarily cut-off one). Kept as similarity instead of
  # dissimilarity as that makes more sense for network plot
  meanSimilarity <- mean(compSimMatFull[lower.tri(compSimMatFull)])

  # Modularity. Done on the FULL matrix, as we don't want it to depend
  # on arbitrarily chosen cut-off. Done with pathways as groups,
  # as algorithmically finding group structure seems like overkill
  # (and I don't understand it). This then measures how separate/modular
  # molecules in different pathways are, which is a bit strange measure,
  # but I can't figure out how modularity would be more useful in another way.
  # We also calculate the number of pathways. Both these are only done if
  # npcTable is supplied.
  nNpcPathways <- NA
  modularity <- NA
  if(!is.null(npcTable)) {
    nNpcPathways <- length(unique(npcTable$pathway[!is.na(npcTable$pathway)]))

    modNet <- igraph::graph.adjacency(compSimMatFull,
                                      mode = "undirected",
                                      weighted = TRUE)
    membership <- as.numeric(as.factor(npcTable$pathway))
    # Give NA their own category, otherwise modularity() crashes the session
    membership[is.na(membership)] <- max(membership, na.rm = TRUE) + 1
    modularity <- igraph::modularity(modNet, membership)
  }

  networkOutput <- list("networkObject" = networkObject,
                        "nCompounds" = nCompounds,
                        "meanSimilarity" = meanSimilarity,
                        "nNpcPathways" = nNpcPathways,
                        "modularity" = modularity)
  return(networkOutput)
}
