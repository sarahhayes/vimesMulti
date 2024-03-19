rm(list = ls())

library(devtools)

load_all()

g1_obs <- 313
g2_obs <- 236
g1_rr <- 0.25
g2_rr <- 0.10
n = 10000000
q <- 0.95
si_mean <- 27.8175438596491
si_sd <- 36.8565433014125
params_temporal <- c(si_mean, si_sd, si_mean, si_sd, si_mean, si_sd)
assort_mix <- 2.4

set.seed(1234)
# first using the props from all the transmissions
cuts_temporal <-  get_quantiles_multi(d_type = "temporal", distrib = "lognormal",
                                      obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                      params = params_temporal,
                                      n = n, q = q, assort_mix = assort_mix)

dist_mean <- 0.87
params_spatial <- c(dist_mean, NA, dist_mean, NA,  dist_mean, NA)

set.seed(1234)
cuts_spatial <- get_quantiles_multi(d_type = "spatial",distrib = "rayleigh",
                                    obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                    params = params_spatial,
                                    n = n, q = q, assort_mix = assort_mix)

cuts_temporal
cuts_spatial


## Need the case data
case_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_time_diffs.csv")
case_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_dists.csv")
spp_vect <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_sp_vect.csv")
spp_vect <- as.factor(spp_vect$x)

dat_time <- case_times$x
dat_geo1 <- case_dists$UTM.Easting
dat_geo2 <- case_dists$UTM.Northing


vimes_multi_results_list <- vimes_multi(
  dat_time = dat_time,
  dat_geo1 = dat_geo1,
  dat_geo2 = dat_geo2,
  temporal_cut_offs = cuts_temporal,
  distance_cut_offs = cuts_spatial,
  group_vect = spp_vect,
  graph_opt = vimes::vimes_graph_opt())

vimes_multi_results_list

trans_table_res <- vimes_multi_results_list$trans_tab_results

# Create a data set from the inputs and add the cluster membership

combined_results <- vimes_multi_results_list$combined_results
dataset <- vimes_multi_results_list$dataset
cluster_size_df <- vimes_multi_results_list$cluster_size_df

chisq_val <- sum((combined_results$data_count - combined_results$sim_count)^2/combined_results$sim_count)
# df in chi sq is 1 because we know how many of each species so only need 1 parameter to work out the rest
chisq_val
1-pchisq(chisq_val, df = 1)

# Review the outputs (Can use this to compare to original script)

vimes_multi_results_list$vimes_results_list$clusters$K
mean(vimes_multi_results_list$vimes_results_list$clusters$size)
hist(vimes_multi_results_list$vimes_results_list$clusters$size, col = "pink", xlab = "Size of cluster", breaks = seq(-1,15,1),
     main = "Histogram of cluster sizes")
table(vimes_multi_results_list$vimes_results_list$clusters$size)

# plot by cluster composition
csl_1 <- cluster_size_df %>%
  dplyr::group_by(total, trans_type)%>%
  dplyr::count()

csl_1_trios_up <- csl_1[which(csl_1$total >=3),]

sum(csl_1_trios_up$n) # so 24 clusters of >= 3 in total.
sum(csl_1_trios_up[which(csl_1_trios_up$trans_type %in% c("g1g1", "g2g2")), "n"])


# plot of the clusters
own_cols <- c("red", "purple", "blue")
own_labels <- c("Domestic only", "Mixed", "Wildlife only")

library(ggplot2)
gg_res_1 <- ggplot(csl_1, aes(x = total, y = n, fill = trans_type))  +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = own_cols, name = "Transmission type",
                    #guide = guide_legend(reverse=TRUE),
                    labels = own_labels) +
  theme_classic() +
  theme(axis.title.x = element_text(size=12, face="plain"),
        axis.title.y = element_text(size=12, face="plain"),
        axis.text.x  = element_text(size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"),
        axis.ticks.length = unit(.2, "cm"),
        legend.text = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 12, face = "plain"),
        legend.position = c(0.8,0.7)) +
  ylim(0,60) +
  #  xlim(1,14) +
  labs(title = "A", y = "Number of clusters", x = "Size of cluster")
gg_res_1
