# Want to look at how the uncertainty in the dates and locations could affect our results.
# Will start by looking at scenario 7 with the Rayleigh distribution as this is our most supported scenario.

rm(list = ls())

library(devtools)

load_all()

g1_obs <- 313
g2_obs <- 236
g1_rr <- 0.50
g2_rr <- 0.25
n = 10000000
q <- 0.95
si_mean <- 27.8175438596491
si_sd <- 36.8565433014125
params_temporal <- c(si_mean, si_sd, si_mean, si_sd, si_mean, si_sd)

# We want to run a number of different values for the shape_2 parameter.
# number of values that we want to check
# n_grid <- 91
assort_vect <- seq(1,5,0.1)
#assort_vect <- c(1, 2.4)

set.seed(1234)

cuts_list_temporal <- list()

for (i in 1:length(assort_vect)) {
  cuts_list_temporal[[i]] <- get_quantiles_multi(d_type = "temporal", distrib = "lognormal",
                                                 obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                                 params = params_temporal,
                                                 n = n, q = q, assort_mix = assort_vect[i])
}

# cuts_list_temporal

# Repeat for distance
dist_mean <- 0.87
params_spatial <- c(dist_mean, NA, dist_mean*2.5, NA,  dist_mean*5, NA)

set.seed(1234)

cuts_list_spatial <- list()

for (i in 1:length(assort_vect)) {
  cuts_list_spatial[[i]] <- get_quantiles_multi(d_type = "spatial", distrib = "rayleigh",
                                                obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                                params = params_spatial,
                                                n = n, q = q, assort_mix = assort_vect[i])
}

# cuts_list_spatial

## For this we need to generate multiple data sets using the uncertainty associated with times and locations

## Need the case data
case_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_time_diffs.csv")
case_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_dists.csv")
uncert_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/dates_uncertainty.csv")
uncert_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/distance_uncertainty.csv")

spp_vect <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_sp_vect.csv")
spp_vect <- as.factor(spp_vect$x)


## we need to introduce the variation into the dates and locations.
## produce the number of sets required for the bootstrapping

# set number for bootstrap
b = 100

## dates
set.seed(54)
new_dates_list <- list()

for (i in 1:b) {
  new_dates <- case_times$x + round(runif(n = nrow(uncert_times), min = -1*(uncert_times$x), max = uncert_times$x))
  new_dates[which(new_dates<0)] <- 0 # to make sure none are before the start of the study.
  new_dates_list[[i]] <- new_dates
}


# locations
set.seed(43)

new_locations_list <- list()

for (i in 1:b) {
  degrees <- runif(n = nrow(case_dists), min=0, max = 360)
  dist_diff <- round(runif(n = nrow(case_dists), min = uncert_dists$lower_uncert, max = uncert_dists$upper_uncert))
  new_df <- case_dists
  new_df[,"UTM.Easting"] <- new_df[,"UTM.Easting"] + dist_diff[] * cos(degrees[])
  new_df[,"UTM.Northing"] <- new_df[,"UTM.Northing"] + dist_diff[] * sin(degrees[])
  new_locations_list[[i]] <- new_df
}


#############################

# now run the vimesMulti function for each of the datasets using the specified cut-offs generated
#  and for each of the assortativity parameters within that dataset

chisq_vals <- c()
prop_succesful <- c()
min_chisq <- c()

tictoc::tic()
for (k in 1:length(assort_vect)) {
  #  bootstrap_results_list <- list()
  for (i in 1:b) {
    bootstrap_results <- vimes_multi(
      dat_time = new_dates_list[[i]],
      dat_geo1 = new_locations_list[[i]]$UTM.Easting,
      dat_geo2 = new_locations_list[[i]]$UTM.Northing,
      temporal_cut_offs = cuts_list_temporal[[k]],
      distance_cut_offs = cuts_list_spatial[[k]],
      group_vect = spp_vect,
      graph_opt = vimes::vimes_graph_opt())
    combined_results <- bootstrap_results$combined_results
    chisq_vals[i] <- sum((combined_results$data_count - combined_results$sim_count)^2/combined_results$sim_count)
  }
  prop_succesful[k] <- length(which(chisq_vals < 3.841))/b
  min_chisq[k] <- min(chisq_vals)
}
tictoc::toc()

prop_results <- cbind(assort_vect, prop_succesful, min_chisq)
prop_results <- as.data.frame(prop_results)

#write.csv(prop_results, "for_paper/bootstrap_100_scenario_5_rayleigh.csv", row.names = F)

######################################################
#### Now repeat for the gamma distribution

# can use the temporal cuts as above

# Repeat for distance
dist_mean <- 0.87
dist_sd <- 1.50
params_spatial <- c(dist_mean, dist_sd, dist_mean*2.5, dist_sd*2.5,  dist_mean*5, dist_sd*5)

set.seed(1234)

cuts_list_spatial_gamma <- list()

for (i in 1:length(assort_vect)) {
  cuts_list_spatial_gamma[[i]] <- get_quantiles_multi(d_type = "spatial", distrib = "gamma",
                                                      obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                                      params = params_spatial,
                                                      n = n, q = q, assort_mix = assort_vect[i])
}

## we can use the new datasets generated above


chisq_vals_gamma <- c()
prop_succesful_gamma <- c()
min_chisq_gamma <- c()

tictoc::tic()
for (k in 1:length(assort_vect)) {
  #  bootstrap_results_list <- list()
  for (i in 1:b) {
    bootstrap_results <- vimes_multi(
      dat_time = new_dates_list[[i]],
      dat_geo1 = new_locations_list[[i]]$UTM.Easting,
      dat_geo2 = new_locations_list[[i]]$UTM.Northing,
      temporal_cut_offs = cuts_list_temporal[[k]],
      distance_cut_offs = cuts_list_spatial_gamma[[k]],
      group_vect = spp_vect,
      graph_opt = vimes::vimes_graph_opt())
    combined_results <- bootstrap_results$combined_results
    chisq_vals_gamma[i] <- sum((combined_results$data_count - combined_results$sim_count)^2/combined_results$sim_count)
  }
  prop_succesful_gamma[k] <- length(which(chisq_vals_gamma < 3.841))/b
  min_chisq_gamma[k] <- min(chisq_vals_gamma)
}
tictoc::toc()

prop_results_gamma <- cbind(assort_vect, prop_succesful_gamma, min_chisq_gamma)
prop_results_gamma <- as.data.frame(prop_results_gamma)


#write.csv(prop_results_gamma, "for_paper/bootstrap_100_scenario_5_gamma.csv", row.names = F)


