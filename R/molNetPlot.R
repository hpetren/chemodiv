#' Plot molecular network(s)
#'
#' Function to conveniently create a basic plot of the molecular network
#' created by the \code{\link{molNet}} function. Molecular networks can be
#' used to illustrate the biosynthetic/structural similarity of
#' phytochemical compounds in a sample, while simultaneously visualizing
#' their relative concentrations. In the network, nodes are compounds,
#' with node sizes or node colours representing the relative concentrations
#' of compounds. Edges connects nodes, with edge widths representing
#' compound similarity.
#'
#' The network object from \code{\link{molNet}} and \code{sampleData} have to
#' be supplied. In addition, \code{groupData} and/or an \code{\link{NPCTable}}
#' can be supplied. If an NPCTable is supplied, which is recommended,
#' node colours will represent NPC pathways, and node sizes the relative
#' concentration of the compounds. Edge widths represent compound similarity,
#' and only edges with similarity values above the \code{cutOff} value
#' in the \code{\link{molNet}} function will be plotted.
#' If \code{groupData} is supplied, one network will be created for each group.
#' When \code{groupData} but not an \code{\link{NPCTable}} is supplied,
#' compounds missing (i.e. having a mean of 0) from
#' specific groups are plotted as triangles. When \code{groupData}
#' and an \code{\link{NPCTable}} is supplied, compounds missing from
#' specific groups have a white fill. This is done so that networks are more
#' easy to compare across groups.
#'
#' @param sampleData Data frame with the relative concentration of each
#' compound (column) in every sample (row).
#' @param groupData Grouping data (e.g. population, species etc.).
#' If supplied, a separate network will be created for each group.
#' Should be either a vector, or a data frame with a single column.
#' @param networkObject A network object, as created by the
#' \code{\link{molNet}} function. Note that this is only the network object,
#' which is one of the elements in the list outputted by \code{\link{molNet}}.
#' The network is extracted as \code{molNetOutput$networkObject}.
#' @param npcTable It is optional but recommended to supply an
#' \code{\link{NPCTable}}. This will result in network nodes being coloured
#' by their NPC pathway classification.
#' @param plotNames Indicates if compounds names should be included
#' in the molecular network plot.
#' @param layout Layout used by \code{\link[ggraph]{ggraph}} when creating
#' the network. The default chosen here, \code{"kk"}, is the the Kamada-Kawai
#' layout algorithm which in most cases should produce a visually
#' pleasing network. Another useful option is \code{"circle"}, which puts all
#' nodes in a circle, for easier comparisons between different networks.
#'
#' @return A plot with one or more molecular networks.
#'
#' @export
#'
#' @examples
#' data(minimalSampData)
#' data(minimalCompDis)
#' data(minimalNPCTable)
#' data(minimalMolNet)
#' groups <- c("A", "A", "B", "B")
#' molNetPlot(minimalSampData, minimalMolNet)
#' molNetPlot(minimalSampData, minimalMolNet, groups)
#' molNetPlot(minimalSampData, minimalMolNet, npcTable = minimalNPCTable)
#' molNetPlot(minimalSampData, minimalMolNet, plotNames = TRUE)
#'
#' \dontrun{
#' data(alpinaCompData)
#' data(alpinaSampData)
#' data(alpinaPopData)
#' alpinaNPCTable <- NPCTable(compoundData = alpinaCompData)
#' alpinaCompDis <- compDis(compoundData = alpinaCompData)
#' alpinaMolNet <- molNet(compDisMat = alpinaCompDis$fingerDisMat,
#' npcTable = alpinaNPCTable, cutOff = 0.75)
#' molNetPlot(sampleData = alpinaSampData,
#' networkObject = alpinaMolNet$networkObject,
#' groupData = alpinaPopData, npcTable = alpinaNPCTable)
#' }
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

  if (is.data.frame(groupData)) {
    groupData <- as.vector(groupData[,1])
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
                     legend.text = ggplot2::element_text(size = 14),
                     panel.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())

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
                     legend.text = ggplot2::element_text(size = 14),
                     panel.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())

    networkList <- list(p1)

  } else if (is.null(groupData) && !is.null(npcTable) && !plotNames) {
    # One network with NPC

    # Colours from https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#E6AB02"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                                           fill = npcTable$pathway,
                                           size = compoundMean),
                              shape = 21) +
      ggraph::geom_node_point(ggplot2::aes(size = compoundMean),
                              shape = 21,
                              color = npcTable$col,
                              fill = npcTable$col) +
      ggplot2::scale_size(range = c(8, 16)) +
      ggplot2::labs(fill = "Pathway", width = "Molecular similarity", size = "Proportion") +
      ggplot2::guides(color = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14),
                     panel.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())

    networkList <- list(p1)

  } else if (is.null(groupData) && !is.null(npcTable) && plotNames) {
    # One network with names and NPC

    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#E6AB02"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"

    compoundMean <- colMeans(sampleData)

    p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 3), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                                           fill = npcTable$pathway,
                                           size = compoundMean),
                              shape = 21) +
      ggraph::geom_node_point(ggplot2::aes(size = compoundMean),
                              shape = 21,
                              color = npcTable$col,
                              fill = npcTable$col) +
      ggplot2::scale_size(range = c(8, 16)) +
      ggplot2::guides(color = "none") +
      ggplot2::labs(fill = "Pathway", width = "Molecular similarity", size = "Proportion") +
      ggraph::geom_node_label(ggplot2::aes(label = .data$name), nudge_x = 0, nudge_y = 0.2) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 16),
                     legend.text = ggplot2::element_text(size = 14),
                     panel.background = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())

    networkList <- list(p1)

  } else if (!is.null(groupData) && !is.null(npcTable)) {
    # Both group and NPC

    npcTable$col <- NA
    npcTable$col[is.na(npcTable$pathway)] <- "#666666"
    npcTable$col[npcTable$pathway == "Alkaloids"] <- "#E7298A"
    npcTable$col[npcTable$pathway == "Amino acids and Peptides"] <- "#66A61E"
    npcTable$col[npcTable$pathway == "Carbohydrates"] <- "#E6AB02"
    npcTable$col[npcTable$pathway == "Fatty acids"] <- "#7570B3"
    npcTable$col[npcTable$pathway == "Polyketides"] <- "#A6761D"
    npcTable$col[npcTable$pathway == "Shikimates and Phenylpropanoids"] <- "#1B9E77"
    npcTable$col[npcTable$pathway == "Terpenoids"] <- "#D95F02"

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

        # Fix white colours for 0 values in each group
        groupCol <- as.data.frame(npcTable$col)
        colnames(groupCol) <- "col"
        for(row in 1:nrow(compoundMeanTrans)) {
          if (compoundMeanTrans[row,j] == 0) {
            groupCol$col[row] <- "#FFFFFF"
          }
        }

        p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
          ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                                               fill = npcTable$pathway,
                                               size = compoundMeanTrans[,j]),
                                  shape = 21) +
          ggraph::geom_node_point(ggplot2::aes(size = compoundMeanTrans[,j]),
                                  shape = 21,
                                  color = npcTable$col,
                                  fill = groupCol$col) +
          ggplot2::scale_size(range = c(4, 10)) +
          ggplot2::guides(color = "none") +
          ggplot2::labs(fill = "Pathway", width = "Molecular similarity", size = "Proportion") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) +
          ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 8),
                         panel.background = ggplot2::element_blank(),
                         legend.key = ggplot2::element_blank())

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

        # Fix different shape for 0 values in each group
        groupShape <- data.frame(shape = rep("Pos", nrow(compoundMeanTrans)))
        for(row in 1:nrow(compoundMeanTrans)) {
          if (compoundMeanTrans[row,j] == 0) {
            groupShape$shape[row] <- "Zero"
          }
        }

        p1 <- ggraph::ggraph(graph = networkObject, layout = layout) +
          ggraph::geom_edge_link(ggplot2::aes(width = .data$weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = compoundMeanTrans[,j],
                                               shape = groupShape$shape), size = 10) +
          ggplot2::scale_colour_viridis_c() +
          ggplot2::labs(color = "Proportion", width = "Molecular similarity") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j]) +
          ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
                         legend.text = ggplot2::element_text(size = 8),
                         panel.background = ggplot2::element_blank(),
                         legend.key = ggplot2::element_blank())

        print(p1)
      })
    }
  }
  return(gridExtra::grid.arrange(grobs = networkList,
                                 ncol = ceiling(sqrt(length(networkList)))))
}
