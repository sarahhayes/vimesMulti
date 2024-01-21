
#' Generate dataframe for use in vimes_multi when user specifying thresholds
#'
#' @param d_type one of 'temporal' or 'spatial'
#' @param g1g1_threshold cutoff distance for group 1 to group 1 transmission
#' @param mixed_threshold cutoff distance for mixed transmission
#' @param g2g2_threshold cutoff distance for group 2 to group 2 transmission
#'
#' @return a dataframe of 3 columns and 3 rows suitable for use in vimes_multi
#' @export
#'
#' @examples
#' set_cutoff("spatial", 4,4,4)



set_cutoff <- function(d_type, g1g1_threshold, mixed_threshold, g2g2_threshold){
  df <- as.data.frame(matrix(ncol = 3, nrow = 3))
  if(d_type != "temporal" & d_type != "spatial"){
    msg_error <- "unsupported value for d_type"
    stop(msg_error)
  }

  if(d_type == "temporal"){
    colnames(df) <- c("trans_type", "threshold_temporal", "proportion_sim_temporal")
    df[,"trans_type"] <- c("g1g1", "mixed", "g2g2")
    df[,"threshold_temporal"] <- c(g1g1_threshold, mixed_threshold, g2g2_threshold)
    df[,"proportion_sim_temporal"] <- NA
  }
  if(d_type == "spatial"){
    colnames(df) <- c("trans_type", "threshold_spatial", "proportion_sim_spatial")
    df[,"trans_type"] <- c("g1g1", "mixed", "g2g2")
    df[,"threshold_spatial"] <- c(g1g1_threshold, mixed_threshold, g2g2_threshold)
    df[,"proportion_sim_spatial"] <- NA
  }

  return(df)
}
