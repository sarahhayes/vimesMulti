vimes_prune_multi <-   function(x, cutoff = NULL,
                                graph_opt = vimes::vimes_graph_opt(), ...){

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

  cuts_mat <- matrix(cutoff[grp_mat_numbers], ncol = length(group_vect),
                     nrow = length(group_vect), byrow = F)

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
