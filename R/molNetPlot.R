#' Plot molecular network
#'
#' Function to conveniently create a basic plot of the molecular network
#' created by the \code{\link{molNet}} function
#' (based on the compound dissimilarity matrix). Nodes represent compounds,
#' with size proportional to proportions, and edge widths represent
#' compound similarity.
#'
#' The network object from \code{\link{molNet}}  and sample dataset have to
#' be supplied. In addition, groupData and/or \code{\link{NPCTable}} can
#' be supplied. If groupData is supplied,
#' one network will be created for each group. If and NPCTable is supplied,
#' node colours will represent NPC pathways, and node size the proportional
#' concentration of the compounds. Edge widths represent compound similarity,
#' and only edges with similarity values above the \code{\link{molNet}}
#' function's cut-off will be plotted.
#'
#' @param sampleData Data frame with samples as rows and compounds as columns.
#' @param groupData Grouping data. If supplied, a separate network will be
#' created for each group
#' @param networkObject tidygraph network object, as created by the
#' \code{\link{molNet}} function. Note that this is only the network object,
#' which has to be taken from the list.
#' @param npcTable An \code{\link{NPCTable}} can optionally be supplied.
#' This will result in network nodes being coloured by their
#' NPC pathway classification
#' @param plotNames Indicates if compounds names should be included
#' in the molecular network plot.
#'
#' @return Molecular network(s) created with ggraph
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalMolNet)
#' data(minimalNPCTable)
#' groups <- c("A", "A", "B", "B")
#' molNetPlot(minimalSampData, minimalMolNet)
#' molNetPlot(minimalSampData, minimalMolNet, groups)
#' molNetPlot(minimalSampData, minimalMolNet, minimalNPCTable)
#' molNetPlot(minimalSampData, minimalMolNet, plotNames = TRUE)
molNetPlot <- function(sampleData,
                       networkObject,
                       groupData = NULL,
                       npcTable = NULL,
                       plotNames = FALSE) {

  # If one does not use groupData = , or npcTable = , in function input
  # then these can become mixed as you only have to supply one.
  # I.e. third input could be either the grouping data or the NPC,
  # and if not specific, this will be assigned to groupData
  # This checks what kind was inputted, and makes it correct. Although
  # there should be a better solution for this.
  if((is.null(groupData) | is.null(npcTable)) &
     "pathway" %in% colnames(groupData)) {
    npcTable <- groupData
    groupData <- NULL
  }

  if (!is.null(groupData) & plotNames) {
    stop("Names can only be plotted without grouping data")
  }

  if (is.null(groupData) & is.null(npcTable) & !plotNames) {
    # One network

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
      ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = compoundMean), size = 8) +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(color = "Proportion", width = "Molecular similarity")

    networkList <- list(p1)

  } else if (is.null(groupData) & is.null(npcTable) & plotNames) {
    # One network with names

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
      ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = compoundMean), size = 8) +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
      ggraph::geom_node_label(aes_(label = ~name), nudge_x = 0, nudge_y = 0.1)

    networkList <- list(p1)

  } else if (is.null(groupData) & !is.null(npcTable) & !plotNames) {
    # One network with NPC

    # One of these is probably the best one
    # ggraph(graph = networkObject, layout = "stress")
    # ggraph(graph = networkObject, layout = "igraph", algorithm = "kk")
    # ggraph(graph = networkObject, layout = "igraph", algorithm = "nicely")

    p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
      ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway), size = 8) +
      ggplot2::labs(color = "Pathway", width = "Molecular similarity")

    networkList <- list(p1)

  } else if (is.null(groupData) & !is.null(npcTable) & plotNames) {
    # One network with names and NPC

    p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
      ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway), size = 8) +
      ggplot2::labs(color = "Pathway", width = "Molecular similarity") +
      ggraph::geom_node_label(aes_(label = ~name), nudge_x = 0, nudge_y = 0.1)

    networkList <- list(p1)

  } else if (!is.null(groupData) & !is.null(npcTable)) {
    # Both group and NPC

    # Calculating average per compound per group. Note that this is
    # the average of proportions (as that is indata)
    compoundMean <- stats::aggregate(sampleData,
                                     by = list(Group = groupData),
                                     mean)
    compoundMeanTrans <- t(compoundMean[,2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()

    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
          ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                              size = compoundMeanTrans[,j])) +
          ggplot2::scale_size(range = c(2, 10)) +
          ggplot2::labs(color = "Pathway", width = "Molecular similarity", size = "Proportion") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j])

        print(p1)
      })
    }

  } else {
    # Only group data

    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), mean)
    compoundMeanTrans <- t(compoundMean[,2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()

    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
          ggraph::geom_edge_link(ggplot2::aes_(width = ~weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = compoundMeanTrans[,j]), size = 8) +
          ggplot2::scale_colour_viridis_c() +
          ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j])

        print(p1)
      })
    }
  }

  # Can't figure out how to get this to behave as ggplot plot, i.e.
  # you can save output into variable and when that variable is run the
  # plot is shown. arrangeGrob makes it so that nothing is plotted when
  # function is run, but then have to save function output into variable
  # and run plot(variable) to actually see the plot. Not doing that.
  #a <- gridExtra::arrangeGrob(grobs = networkList, nrow = ceiling(sqrt(length(networkList))))
  #return(a)
  # Sticking with the normal grid.arrange, which always make so a plot
  # is made, even if function output is saved as variable (so function
  # behaves like plot() does).
  return(gridExtra::grid.arrange(grobs = networkList, nrow = ceiling(sqrt(length(networkList)))))
}
