
vimes_trans_multi <- function(x){
  graph <- x[["graph"]]
  graph_df <- igraph::as_data_frame(graph)
  graph_df <- dplyr::select(graph_df, c(from, to))

  graph_df[,"first_grp"] <-  group_vect[as.numeric(graph_df[,"from"])]
  graph_df[,"second_grp"] <-  group_vect[as.numeric(graph_df[,"to"])]

  graph_df[,"time_first_grp"] <-  dat_time[as.numeric(graph_df[,"from"])]
  graph_df[,"time_second_grp"] <-  dat_time[as.numeric(graph_df[,"to"])]


  graph_df[,"trans"] <- paste(graph_df[,"first_grp"], graph_df[,"second_grp"], sep = "")
  graph_df[which(graph_df$trans %in% c("g1g2", "g2g1")), "trans"] <- "mixed"

  return(graph_df)
}

