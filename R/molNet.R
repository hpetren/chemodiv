#' Generate molecular network with properties
#'
#' Function which generates a molecular network object, and some
#' basic properties of the network/dissimilarity matrix
#'
#' Note that the supplied dissimilarity matrix is translated to a
#' similarity matrix, and this is what cutOff values are set for.
#'
#' @param compDisMat Compound dissimilarity matrix.
#' @param npcTable An NPCTable can optionally be supplied generate the
#' number of NPC pathways and calculate network modularity
#' @param cutOff Cut-off value for compound similarities. Any similarity
#' lower than this value will be set to zero when the network is generated,
#' which strongly affects the look of the network. The value can be set
#' manually 0 <= cutOff <= 1 (maybe mention that others use 0.6 for another
#' kind of molecular similarity measure), to the \code{median} dissimilarity
#' value in the compDisMat, or to \code{minPathway} the lowest within-pathway
#' similarity (which allows all within-NPC-pathway to be kept) if an NPCTable
#' is supplied.
#'
#' @return List with a graph (tidygrph/igraph?) object, the number of compounds,
#' the mean similarity of the compounds, number of NPC pathways and modularity.
#'
#' @export
#'
#' @examples
#' data(minimalCompDis)
#' data(minimalNPCTable)
#' molNet(minimalCompDis)
#' molNet(minimalCompDis, minimalNPCTable, cutOff = 0)
molNet <- function(compDisMat,
                   npcTable = NULL,
                   cutOff = "median") {

  if (is.numeric(cutOff) & ((cutOff < 0) | (cutOff > 1))){
    stop("Numeric values for cutOff must be between 0 and 1")
  }
  if (is.character(cutOff) & !(cutOff == "median" | cutOff == "minPathway")){
    stop("cutOff must be a value between 0 and 1, median or minPathway")
  }

  # From dissimilarity to similarity matrix
  compSimMat <- 1 - compDisMat

  # Set diagonal to 0
  diag(compSimMat) <- 0

  # Keep the full similarity matrix
  compSimMatFull <- compSimMat


  # Now we set cut-off point. This is median, manual,
  # or the lowest within-pathway similarity

  if (cutOff == "median") {

    medianSim <- stats::median(compSimMat[lower.tri(compSimMat)])
    compSimMat[compSimMat < medianSim] <- 0

    message("Using median similarity as cut-off")

  } else if (is.numeric(cutOff)) { # Should maybe test if 0 <= cutOff <= 1

    compSimMat[compSimMat < cutOff] <- 0

    message(paste0("Using cut-off value = ", cutOff))

  } else if (cutOff == "minPathway" & !is.null(npcTable)) {

    minIntraPath <- max(compSimMat)
    simValue <- NA

    for (col in 1:(ncol(compSimMat)-1)) {
      for (row in (col+1):nrow(compSimMat)) {
        # Going through lower triangle

        if (!is.na(npcTable$pathway[col]) & !is.na(npcTable$pathway[row])) {
          # Only if both are not NA

          if (npcTable$pathway[col] == npcTable$pathway[row]) {
            # If these compound are from the same pathway

            simValue <- compSimMat[row,col]

            if (simValue < minIntraPath) {
              # If this value is smaller than previous one, replace

              minIntraPath <- simValue

            }
          }
        }
      }
    }

    compSimMat[compSimMat < minIntraPath] <- 0

    message("Using lowest within-pathway similarity as cut-off")
  } else if (cutOff == "minPathway" & is.null(npcTable)) {
    # If trying to use minPathway but not having npcTable
    stop("Using minPathway as cut-off requires npcTable")
  }


  # Creating network
  networkObject <- tidygraph::as_tbl_graph(compSimMat)


  # Calculating network indices and other stuff ###

  # Number of compounds
  nCompounds <- nrow(compSimMatFull)

  # Mean similarity. This is done on FULL matrix (not cut-off) as
  # that makes more sense (you want the mean similarity of the full bouquet,
  # not an arbitarily cut-off one).
  # Should we chagne this to mean DISsimilarity?
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

    modNet <- igraph::graph.adjacency(compSimMatFull, mode = "undirected", weighted = TRUE)
    membership <- as.numeric(as.factor(npcTable$pathway))
    # Give NA their own category, otherwise modularity() crashes the session
    membership[is.na(membership)] <- max(membership, na.rm = TRUE) + 1
    modularity <- igraph::modularity(modNet, membership)

  }

  # MODULARITY IS REEEEEEEALLY TRICKY
  # I didn't understand it correctly at first. Newman 2006 PNAS eq. 11
  # actually requires that one defines if nodes are in the same
  # or different groups. So it seems modularity requires defintion
  # of groups.
  # The bipartite function computeModules is probably not useful here,
  # as it is for a bipartite network, and I have no clue of what it
  # actually does.
  # multilevel.community() can apply groups by algorithm and then
  # do modularity, but I don't understand it and it gives two values sometimes
  # which is strange.
  # So doing modularity() with pathway groups, as that is atleast understandable
  # Modularity links:
  # https://bookdown.org/markhoff/social_network_analysis/finding-groups-in-networks.html
  # https://igraph.org/r/pdf/1.2.7/igraph.pdf
  # https://stackoverflow.com/questions/5519349/network-modularity-calculations-in-r
  # https://dshizuka.github.io/networkanalysis/05_community.html
  # http://networksciencebook.com/chapter/9#modularity


  # Make and return list
  networkOutput <- list("networkObject" = networkObject,
                        "nCompounds" = nCompounds,
                        "meanSimilarity" = meanSimilarity,
                        "nNpcPathways" = nNpcPathways,
                        "modularity" = modularity)

  return(networkOutput)

}
