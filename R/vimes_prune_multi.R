#' Prune the graphs using group-specific values
#'
#' @param x A vimes data object
#' @param cutoff The cut-off values to be used for pruning
#' @param group_vect Vector of groups with same length as observations
#' @param graph_opt The settings for the graph. Imported from vimes.
#' @param ... Any additional parameters relating to graph options
#'
#' @return list containing the pruned graph, cluster definition, and cutoff values used
#' @export
#'
#' @import vimes igraph fields
#' @examples
#' g1_obs <- 300 # no of cases observed for species 1
#' g2_obs <- 200 # no of cases observed for species 2
#' obs <- g1_obs + g2_obs
#' group_vect <- as.factor(sample(c(rep("g1", g1_obs), rep("g2", g2_obs)), replace = FALSE))
#'
#' dat_time <- round(runif(obs, min = 0, max = 400)) #days from start of study
#' time_dist <- dist(dat_time)
#'
#' dat_easting <- runif(obs, min = 372900, max = 656900) # this represents the distance data.
#' dat_northing <- runif(obs, min = 8744900, max = 9075900) # this represents the distance data.
#'
#' dat_dist <- as.matrix(c(dat_easting, dat_northing))
#' dat_dist <- fields::rdist(dat_dist)/1000 # convert to km
#'
#' x <- list(time_dist, dat_dist)
#' x <- vimes::vimes_data(x)
#' x <- x[[1]]
#' cutoff <- as.list(c(142, 143, 143))
#'
#' vimes_prune_multi(
#' x=x,
#' cutoff = cutoff,
#' group_vect = group_vect,
#' graph_opt = vimes::vimes_graph_opt()
#' )
#'


vimes_prune_multi <-   function(x, cutoff = NULL, group_vect,
                                graph_opt = vimes::vimes_graph_opt(),...){

  ## CHECKS ##
  if (is.null(x)) {
    stop("input data is NULL")
  }

  ## BUILD GRAPH ##

  ## In the following we create a pruned graph, derive corresponding
  ## clusters, create new graphical attributes for the graph (mostly
  ## coloring clusters).

  # create the species matrix from the species vector

  row_mat <- matrix(group_vect, nrow = length(group_vect), ncol = length(group_vect),
                    byrow = T)


  col_mat <- matrix(group_vect, nrow = length(group_vect), ncol = length(group_vect),
                    byrow = F)


  grp_mat <- matrix(paste0(row_mat, col_mat), nrow = length(group_vect),
                    ncol = length(group_vect),  byrow = F)

  # now we need to change the matrix values to be 1-3

  grp_mat_numbers <- matrix(2,ncol = length(group_vect), nrow = length(group_vect))

  grp_mat_numbers[which(grp_mat == "g1g1")] <- 1
  grp_mat_numbers[which(grp_mat == "g2g2")] <- 3

  cuts_mat <- matrix(cutoff[grp_mat_numbers], ncol = length(group_vect), nrow = length(group_vect), byrow = F)

  new_x <- 1 - (as.matrix(x) > cuts_mat)

  g <- igraph::graph.adjacency(new_x, mode = "undirected",
                               weighted = TRUE, diag = FALSE)

  ## find clusters ##
  groups <- igraph::clusters(g)
  names(groups) <- c("membership", "size", "K")

  ## add cluster colors
  groups$color <- graph_opt$col_pal(groups$K)
  names(groups$color) <- 1:groups$K

  ## Here we add new graphical properties to the graph that will
  ## ultimately be returned.

  g <- vimes:::set_igraph_opt(g, graph_opt)

  ## The returned output should be a self-sufficient list containing
  ## the pruned graph, cluster definition, and cutoff values used.

  out <- list(graph = g, clusters = groups, cutoff = cutoff)
}
