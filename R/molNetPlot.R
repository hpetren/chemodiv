#' Plot molecular network
#'
#' Function to conveniently create a basic plot of the molecular network created
#' by the molNet function (based on the compound dissimilarity matrix).
#'
#' The network object and sample dataset have to be supplied. In addition,
#' groupData and/or NPCTable must be supplied. If groupData is supplied,
#' one network will be created for each group, with node colours represtnging
#' the proportional concentration of the compounds. If and NPCTable is supplied,
#' node colours will represent NPC pathways, and node size the proportional
#' concentration of the compounds. Edge widths represent compound similarity,
#' and only edges with similarity values above the molNet function's cut-off
#' will be plotted
#'
#' @param sampleData Dataframe with samples as rows and compounds as columns.
#' @param groupData Grouping data. If supplied, a separate network will be
#' created for each group
#' @param networkObject tidygraph network object, as created by molNet function
#' @param npcTable An NPCTable can optionally be supplied. This will result
#' in network nodes being coloured by their NPC pathway classification
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
#' molNetPlot(minimalSampData, minimalMolNet, groups)
#' molNetPlot(minimalSampData, minimalMolNet, minimalNPCTable)
#' molNetPlot(minimalSampData, minimalMolNet, groups, minimalNPCTable)
molNetPlot <- function(sampleData,
                       networkObject,
                       groupData = NULL,
                       npcTable = NULL) {

  # If one does not use groupData = , or npcTable = , in function input
  # then these can become mixed as you onle have to supply one.
  # I.e. third input could be either the grouping data or the NPC,
  # and if not specific, this will be assigned to groupData
  # This checks what kind was inputted, and makes it correct. Although
  # there should be a better solution for this.
  if(is.null(groupData) | is.null(npcTable)){
    if("pathway" %in% colnames(groupData)){
      npcTable <- groupData
      groupData <- NULL
    }
  }

  if (is.null(groupData) & is.null(npcTable)) { # Not making network here
    stop("Provide groupData and/or npcTable")

  } else if (is.null(groupData) & !is.null(npcTable)) { # One network

    # One of these is probably the best one
    # ggraph(graph = networkObject, layout = "stress")
    # ggraph(graph = networkObject, layout = "igraph", algorithm = "kk")
    # ggraph(graph = networkObject, layout = "igraph", algorithm = "nicely")

    p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
      ggraph::geom_edge_link(ggplot2::aes(width = networkObject$weight), edge_color = "grey40") +
      ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
      ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway), size = 8) +
      ggplot2::labs(color = "Pathway", width = "Molecular similarity")


    networkList <- list(p1)



  } else if (!is.null(groupData) & !is.null(npcTable)) { # Both group and NPC

    # Calcularing average per compound per group. Note that this is currently
    # the average of proportions (as that is indata)
    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), mean)
    compoundMeanTrans <- t(compoundMean[,2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()


    # For the for-loops below, unless we do some very strange things that I don't
    # understand the last run in the loop overwrites the other ones.
    # https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
    for (j in 1:ncol(compoundMeanTrans)) {

      networkList[[colnames(compoundMeanTrans)[j]]] <- local({

        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
          ggraph::geom_edge_link(ggplot2::aes(width = networkObject$weight), edge_color = "grey40") +
          ggraph::scale_edge_width(range = c(0.3, 2), name = "Molecular similarity") +
          ggraph::geom_node_point(ggplot2::aes(color = npcTable$pathway,
                              size = compoundMeanTrans[,j])) +
          ggplot2::scale_size(range = c(2, 10)) +
          ggplot2::labs(color = "Pathway", width = "Molecular similarity", size = "Proportion") +
          ggplot2::ggtitle(colnames(compoundMeanTrans)[j])

        print(p1)
      })
    }

  } else { # Only group data

    # Calcularing average per compound per group. Note that this is currently
    # the average of proportions (as that is indata)
    compoundMean <- stats::aggregate(sampleData, by = list(Group = groupData), mean)
    compoundMeanTrans <- t(compoundMean[,2:ncol(compoundMean)])
    colnames(compoundMeanTrans) <- compoundMean$Group
    compoundMeanTrans <- as.data.frame(compoundMeanTrans)

    networkList <- list()

    for (j in 1:ncol(compoundMeanTrans)) {

      networkList[[colnames(compoundMeanTrans)[j]]] <- local({

        j <- j

        p1 <- ggraph::ggraph(graph = networkObject, layout = "igraph", algorithm = "kk") +
          ggraph::geom_edge_link(ggplot2::aes(width = networkObject$weight), edge_color = "grey40") +
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
  # function is run, but you can see the plot then running a.
  #a <- gridExtra::arrangeGrob(grobs = networkList, nrow = ceiling(sqrt(length(networkList))))
  #return(a)
  # So sticking with the normal grid.arrange, which always make so a plot
  # is made, even if function output is saved as variable (so function
  # behaves more like plot() does).
  return(gridExtra::grid.arrange(grobs = networkList, nrow = ceiling(sqrt(length(networkList)))))
}
