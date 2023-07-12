

get_quantiles_multi <- function(d_type, distrib, obs, rr, n,
                                params, q, assort_mix) {

 ## Calculate the actual number of cases for each species
  g1_n <- obs[1]/rr[1]
  g2_n <- obs[2]/rr[2]

  # use this to calculate the proportions of each species
  g1_prop <- g1_n/(g1_n+g2_n)
  g2_prop <- 1 - g1_prop

  if(g1_prop + g2_prop != 1){
    msg <- "Proportions do not sum to 1"
    stop(msg) # ensure sums to 1.
  }

  #
  g1_sim_n <- round(n*g1_prop)
  g2_sim_n <- n - g1_sim_n  # avoid fractions of individuals and it not adding to n

  if(g1_sim_n + g2_sim_n != n) {
    msg2 <- "Cases do not sum to specified value of n"
    stop(msg2)
  }


  ### Assign the order that the groups will occur in

  x <- data.frame(sp = 'g1', loc = rbeta(n = g1_sim_n, shape1 = 1, shape2 = assort_mix))
  y <- data.frame(sp = 'g2', loc = rbeta(n = g2_sim_n, shape1 = assort_mix, shape2 = 1))

  z <- rbind(x,y)
  z <- z[order(z[,2]),]

  all_ind <- data.frame(grp = z$sp)

  all_ind$detected <- runif(n= n, 0, 1) < ((all_ind$grp == 'g1') * rr[1] +
                                             (all_ind$grp == 'g2') * rr[2])
  all_ind$trans <- c(paste(all_ind$grp[1:n-1], all_ind$grp[2:n], sep = ""),NA)


  if(d_type == "temporal"){

    if(distrib == "gamma"){
      gam_parms_g1g1 <- epitrix::gamma_mucv2shapescale(mu = params[1],
                                                       cv = params[2]/params[1])
      gam_parms_g1g2 <- epitrix::gamma_mucv2shapescale(mu = params[3],
                                                       cv = params[4]/params[3])
      gam_parms_g2g2 <- epitrix::gamma_mucv2shapescale(mu = params[5],
                                                       cv = params[6]/params[5])


      shape_vect <- c(gam_parms_g1g1[1], gam_parms_g1g2[1], gam_parms_g2g2[1]) # s1s1,s1s2,s2s1,s2s2
      scale_vect <- c(gam_parms_g1g1[2], gam_parms_g1g2[2], gam_parms_g2g2[2])

      all_ind$shape <- (all_ind$trans == 'g1g1') * assort_mix[[1]] +
        (all_ind$trans == 'g1g2') * shape_vect[[2]] +
        (all_ind$trans == 'g2g1') * shape_vect[[2]] +
        (all_ind$trans == 'g2g2') * shape_vect[[3]]

      all_ind$scale <- (all_ind$trans == 'g1g1') * scale_vect[[1]] +
        (all_ind$trans == 'g1g2') * scale_vect[[2]] +
        (all_ind$trans == 'g2g1') * scale_vect[[2]] +
        (all_ind$trans == 'g2g2') * scale_vect[[3]]

      all_ind <- all_ind[1:nrow(all_ind)-1,]

      all_ind$distance <- c(rgamma(n-1, shape = all_ind$shape, scale = all_ind$scale))

    }

    if(distrib == "lognormal"){

      meanlog_g1g1 <- log(params[1]/(sqrt(1 + params[2]^2/params[1]^2)))
      sdlog_g1g1 <- sqrt(log(1 + params[2]^2/params[1]^2))

      meanlog_g1g2 <- log(params[3]/(sqrt(1 + params[4]^2/params[3]^2)))
      sdlog_g1g2 <- sqrt(log(1 + params[4]^2/params[3]^2))

      meanlog_g2g2 <- log(params[5]/(sqrt(1 + params[6]^2/params[5]^2)))
      sdlog_g2g2 <- sqrt(log(1 + params[6]^2/params[5]^2))

      meanlog_vect <- c(meanlog_g1g1, meanlog_g1g2,
                        meanlog_g2g2)
      sdlog_vect <- c(sdlog_g1g1, sdlog_g1g2,
                      sdlog_g2g2)

      all_ind$meanlog <- (all_ind$trans == 'g1g1') * meanlog_vect[[1]] +
        (all_ind$trans == 'g1g2') * meanlog_vect[[2]] +
        (all_ind$trans == 'g2g1') * meanlog_vect[[2]] +
        (all_ind$trans == 'g2g2') * meanlog_vect[[3]]

      all_ind$sdlog <- (all_ind$trans == 'g1g1') * sdlog_vect[[1]] +
        (all_ind$trans == 'g1g2') * sdlog_vect[[2]] +
        (all_ind$trans == 'g2g1') * sdlog_vect[[2]] +
        (all_ind$trans == 'g2g2') * sdlog_vect[[3]]

      all_ind <- all_ind[1:nrow(all_ind)-1,] # remove the last row with NA

      all_ind$distance <- c(rlnorm(n-1, meanlog = all_ind$meanlog, sdlog = all_ind$sdlog))

    }

    # work out the cumulative SI
    all_ind$cum = cumsum(all_ind$distance)


    # To match up the types of transmission in later steps we need to shift the measurements down
    # so we are adding up the correct types of transmissions.
    all_ind$cum[2:nrow(all_ind)] <- all_ind$cum[1:nrow(all_ind)-1]
    all_ind[1,"cum"] <- 0


    ## Separate out the different types of transmission
    obs <- all_ind[which(all_ind$detected == TRUE),] # select just the detected cases
    obs$diff <- c(diff(obs$cum), NA) # we don't have the group for the final transmission
    # as don't know what this would be so use NA - assume individual 1 - individual 2 and the SI for 2 relates to this transmission

    obs_n <- nrow(obs)
    obs$obs_trans <- NA
    obs$obs_trans <- c(paste(obs$grp[1:obs_n-1], obs$grp[2:obs_n], sep = ""),NA)

  }

  if(d_type == "spatial") {

    ray_sig_g1g1 <- params[1] / sqrt(acos(-1)/2)
    ray_sig_g1g2 <- params[3] / sqrt(acos(-1)/2)
    ray_sig_g2g2 <- params[5] / sqrt(acos(-1)/2)

    ray_sig_vect <- c(ray_sig_g1g1, ray_sig_g1g2,  ray_sig_g2g2)

    all_ind$ray_sig <- (all_ind$trans == 'g1g1') * ray_sig_vect[[1]] +
      (all_ind$trans == 'g1g2') * ray_sig_vect[[2]] +
      (all_ind$trans == 'g2g1') * ray_sig_vect[[2]] +
      (all_ind$trans == 'g2g2') * ray_sig_vect[[3]]

    all_ind <- all_ind[1:nrow(all_ind)-1,] # remove the last line as NA for transmission

    all_ind$x_rel <- c(rnorm(n-1, mean = 0, sd = all_ind$ray_sig))
    all_ind$y_rel <- c(rnorm(n-1, mean = 0, sd = all_ind$ray_sig))

    # all distances
    all_ind$True_dist <- sqrt(all_ind$x_rel^2+all_ind$y_rel^2)  # using trigonometry of a^2 + b^2 = c^2

    # absolute positions - work out the actual position of each animal rather than position relative to others which is what we have until now
    all_ind$x_abs = cumsum(all_ind$x_rel)
    all_ind$y_abs = cumsum(all_ind$y_rel)

    all_ind$x_abs[2:nrow(all_ind)] <- all_ind$x_abs[1:nrow(all_ind)-1]
    all_ind[1,"x_abs"] <- 0

    all_ind$y_abs[2:nrow(all_ind)] <- all_ind$y_abs[1:nrow(all_ind)-1]
    all_ind[1,"y_abs"] <- 0

    #Extract the case numbers of those observed
    obs <- all_ind[which(all_ind$detected == TRUE),] # select just the detected cases
    obs$new_relative_position_x <- c(diff(obs$x_abs),NA)
    obs$new_relative_position_y <- c(diff(obs$y_abs), NA)

    obs$diff <- c(sqrt(obs$new_relative_position_x^2+obs$new_relative_position_y^2))
    # as don't know what this would be so use NA - assume animal 1 - animal 2 and the dist for 2 relates to this transmission

    obs_n <- nrow(obs)
    obs$obs_trans <- NA
    obs$obs_trans <- c(paste(obs$grp[1:obs_n-1], obs$grp[2:obs_n], sep = ""),NA)

  }

  f_g1g1<- which(obs$obs_trans == "g1g1")
  f_g2g2<- which(obs$obs_trans == "g2g2")
  f_g1g2<- which(obs$obs_trans == "g1g2"| obs$obs_trans == "g2g1" )

  #work out the proportion of each type of transmission

  prop_g1g1 <- length(f_g1g1)/(obs_n-1)
  prop_g2g2 <- length(f_g2g2)/(obs_n-1)
  prop_g1g2 <- length(f_g1g2)/(obs_n-1)


  threshold_sim <- quantile(obs$diff, q, na.rm = T)
  threshold_g1g1 <- quantile(obs[f_g1g1, "diff"], q, na.rm = T)
  threshold_g2g2 <- quantile(obs[f_g2g2, "diff"], q, na.rm = T)
  threshold_g1g2 <- quantile(obs[f_g1g2, "diff"], q, na.rm = T)

  f_g1g1_below_quant <- which(obs$obs_trans == "g1g1" & obs$diff <= threshold_sim)
  f_g2g2_below_quant <- which(obs$obs_trans == "g2g2" & obs$diff <= threshold_sim)
  f_g1g2_below_quant <- which(obs$obs_trans == "g1g2" & obs$diff <= threshold_sim |
                                obs$obs_trans == "g2g1" & obs$diff <= threshold_sim)

  #work out the proportion of each type of transmission

  n_below_quant <- length(which(obs$diff <= threshold_sim))

  prop_g1g1_below_quant <- length(f_g1g1_below_quant)/(n_below_quant)
  prop_g1g2_below_quant <- length(f_g1g2_below_quant)/(n_below_quant)
  prop_g2g2_below_quant <- length(f_g2g2_below_quant)/(n_below_quant)

  res <- as.data.frame(matrix(ncol = 3, nrow = 3))

  if(d_type == "temporal"){
    colnames(res) <- c("trans_type", "threshold_temporal", "proportion_sim_temporal")
    res[,"trans_type"] <- c("g1g1","mixed","g2g2")
    res[,"threshold_temporal"] <- c(threshold_g1g1, threshold_g1g2, threshold_g2g2)
    res[,"proportion_sim_temporal"] <- c(prop_g1g1_below_quant, prop_g1g2_below_quant, prop_g2g2_below_quant)
  }

  if(d_type == "spatial"){
    colnames(res) <- c("trans_type", "threshold_spatial", "proportion_sim_spatial")
    res[,"trans_type"] <- c("g1g1","mixed","g2g2")
    res[,"threshold_spatial"] <- c(threshold_g1g1, threshold_g1g2, threshold_g2g2)
    res[,"proportion_sim_spatial"] <- c(prop_g1g1_below_quant, prop_g1g2_below_quant, prop_g2g2_below_quant)
  }

  return(res)

}
