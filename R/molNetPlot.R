#' Plot molecular network
#'
#' Function to conveniently create a basic plot of the molecular network
#' created by the \code{\link{molNet}} function. Nodes are compounds,
#' with node size or node colour representing proportional concentration
#' of the compounds. Edge widths represent compound similarity.
#'
#' The network object from \code{\link{molNet}} and sample data frame have to
#' be supplied. In addition, groupData and/or \code{\link{NPCTable}} can
#' be supplied. If groupData is supplied, one network will be created
#' for each group. If an NPCTable is supplied (which is recommended),
#' node colours will represent NPC pathways, and node size the proportional
#' concentration of the compounds. Edge widths represent compound similarity,
#' and only edges with similarity values above the \code{\link{molNet}}
#' function's cut-off will be plotted.
#'
#' @param sampleData Data frame with samples as rows and compounds as columns.
#' @param groupData Grouping data. If supplied, a separate network will be
#' created for each group.
#' @param networkObject tidygraph network object, as created by the
#' \code{\link{molNet}} function. Note that this is only the network object,
#' which is one of the elements in the list outputted by \code{\link{molNet}}.
#' The network is extracted as molNetOutput$networkObject.
#' @param npcTable It is recommended but optional to supply an
#' \code{\link{NPCTable}} This will result in network nodes being coloured
#' by their NPC pathway classification.
#' @param plotNames Indicates if compounds names should be included
#' in the molecular network plot.
#' @param layout Layout used by \code{\link[ggraph]{ggraph}} when creating
#' the network. The default chosen here, \code{"kk"}, is the the Kamada-Kawai
#' layout algorithm which in most cases should produce a visually
#' pleasing network. Another useful option is \code{"circle"}, which puts all
#' nodes in a circle, for easier comparisons between different networks.
#'
#' @return Molecular network(s) created with ggraph.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' data(minimalNPCTable)
#' groups <- c("A", "A", "B", "B")
#' minimalMolNet <- molNet(minimalCompDis)
#' molNetPlot(minimalSampData, minimalMolNet$networkObject)
#' molNetPlot(minimalSampData, minimalMolNet$networkObject, groups)
#' molNetPlot(minimalSampData, minimalMolNet$networkObject, npcTable = minimalNPCTable)
#' molNetPlot(minimalSampData, minimalMolNet$networkObject, plotNames = TRUE)
molNetPlot <- function(sampleData,
                       networkObject,
                       groupData = NULL,
                       npcTable = NULL,
                       plotNames = FALSE,
                       layout = "kk") {

  if (!is.null(groupData) && plotNames) {
    stop("Names can only be plotted without grouping data.")
  }
  if (is.null(npcTable)) {
    message("It is recommended to include an npcTable for an improved
            network visualization.")
  }

  if (is.null(groupData) && is.null(npcTable) && !plotNames) {
    # One network

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = compoundMean), size = 16) +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14))

    networkList <- list(p1)

  } else if (is.null(groupData) && is.null(npcTable) && plotNames) {
    # One network with names

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = compoundMean), size = 16) +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
      ggraph::geom_node_label(ggplot2::aes(label = .data$name), nudge_x = 0, nudge_y = 0.2) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14))

    networkList <- list(p1)

  } else if (is.null(groupData) && !is.null(npcTable) && !plotNames) {
    # One network with NPC

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                                           size = compoundMean)) +
      ggplot2::scale_size(range = c(8, 16)) +
      ggplot2::labs(color = "Pathway", width = "Molecular similarity", size = "Proportion") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14))

    networkList <- list(p1)

  } else if (is.null(groupData) && !is.null(npcTable) && plotNames) {
    # One network with names and NPC

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                                           size = compoundMean)) +
      ggplot2::scale_size(range = c(8, 16)) +
      ggplot2::labs(color = "Pathway", width = "Molecular similarity", size = "Proportion") +
      ggraph::geom_node_label(ggplot2::aes(label = .data$name), nudge_x = 0, nudge_y = 0.2) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14))

    networkList <- list(p1)

  } else if (!is.null(groupData) && !is.null(npcTable)) {
    # Both group and NPC

    # Calculating average per compound per group. Note that this is
    # the average of proportions (as that is indata)
    compoundMean <- stats::aggregate(sampleData,
                                     by = list(Group = groupData),
                                     mean)
    compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()

    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
          ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                              size = compoundMeanTrans[,j])) +
          ggplot2::scale_size(range = c(4, 10)) +
          ggplot2::labs(color = "Pathway", width = "Molecular similarity", size = "Proportion") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) +
          ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 8))


        print(p1)
      })
    }

  } else {
    # Only group data

    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), mean)
    compoundMeanTrans <- t(compoundMean[, 2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()

    for (j in 1:ncol(compoundMeanTrans)) {
      networkList[[colnames(compoundMeanTrans)[j]]] <- local({
        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
          ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = compoundMeanTrans[,j]), size = 10) +
          ggplot2::scale_colour_viridis_c() +
          ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) +
          ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 8))

        print(p1)
      })
    }
  }
  return(gridExtra::grid.arrange(grobs = networkList,
                                 nrow = ceiling(sqrt(length(networkList)))))
}
