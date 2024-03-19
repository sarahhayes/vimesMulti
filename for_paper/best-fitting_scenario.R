### Extracting the cluster results from the best-fitting scenario
### This is scenario 2, gamma distribution and assortativity 2.2

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

assort_vect <- 2.2

set.seed(1234)

cuts_list_temporal <- list()

for (i in 1:length(assort_vect)) {
  cuts_list_temporal[[i]] <- get_quantiles_multi(d_type = "temporal", distrib = "lognormal",
                                                 obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                                 params = params_temporal,
                                                 n = n, q = q, assort_mix = assort_vect[i])
}

# Repeat for distance
dist_mean <- 0.87
dist_sd <- 1.50
params_spatial <- c(dist_mean, dist_sd, dist_mean, dist_sd,  dist_mean, dist_sd)

set.seed(1234)

cuts_list_spatial <- list()

for (i in 1:length(assort_vect)) {
  cuts_list_spatial[[i]] <- get_quantiles_multi(d_type = "spatial", distrib = "gamma",
                                                obs = c(g1_obs, g2_obs), rr = c(g1_rr, g2_rr),
                                                params = params_spatial,
                                                n = n, q = q, assort_mix = assort_vect[i])
}

cuts_list_spatial[[1]]


## Need the case data
case_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_time_diffs.csv")
case_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_dists.csv")
spp_vect <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_sp_vect.csv")
spp_vect <- as.factor(spp_vect$x)

dat_time <- case_times$x
dat_geo1 <- case_dists$UTM.Easting
dat_geo2 <- case_dists$UTM.Northing


# run the vimesMulti function for each of the sets of cut-offs generated
vimes_multi_results_list <- list()

for (i in 1:length(assort_vect)) {
  vimes_multi_results_list[[i]] <- vimes_multi(
    dat_time = dat_time,
    dat_geo1 = dat_geo1,
    dat_geo2 = dat_geo2,
    temporal_cut_offs = cuts_list_temporal[[i]],
    distance_cut_offs = cuts_list_spatial[[i]],
    group_vect = spp_vect,
    graph_opt = vimes::vimes_graph_opt())
}

vimes_multi_results_list[[1]]$combined_results

# We want to calculate the chi-squared value for the different scenarios.

combined_results <- vimes_multi_results_list[[1]]$combined_results
chisq_vals <- sum((combined_results$data_count - combined_results$sim_count)^2/combined_results$sim_count)

# Create a data set from the inputs and add the cluster membership

dataset <- vimes_multi_results_list[[1]]$dataset
cluster_size_df <- vimes_multi_results_list[[1]]$cluster_size_df

# Review the outputs (Can use this to compare to original script)

vimes_multi_results_list[[1]]$vimes_results_list$clusters$K
mean(vimes_multi_results_list[[1]]$vimes_results_list$clusters$size)
hist(vimes_multi_results_list[[1]]$vimes_results_list$clusters$size, col = "pink", xlab = "Size of cluster",
     breaks = seq(-1,20,1),
     main = "Histogram of cluster sizes")
# barplot might be better as discrete values
# however, do barplot later coloured by transmission type which is preferable anyway.


table(vimes_multi_results_list[[1]]$vimes_results_list$clusters$size)

# plot by cluster composition
csl_1 <- cluster_size_df %>%
  dplyr::group_by(total, trans_type)%>%
  dplyr::count()

csl_1_trios_up <- csl_1[which(csl_1$total >=3),]

sum(csl_1_trios_up$n) # so 24 clusters of >= 3 in total.
sum(csl_1_trios_up[which(csl_1_trios_up$trans_type %in% c("g1g1", "g2g2")), "n"])


# plot of the clusters
own_cols <- c("red",  "blue", "purple")
own_labels <- c("Domestic only",  "Wildlife only", "Mixed")

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


## Look at the singletons

singles <- as.data.frame(dataset %>% group_by(cluster_no) %>% dplyr::filter(n() == 1))

singles_counts <- hist(singles$dat_time, breaks =  seq(0, 3450, by = 30))$count
plot(counts)

thirty_day_breaks <- seq(0, 3450, by = 30)
year_breaks <- c(0, 365, 731, 1096, 1461, 1826, 2192, 2557, 2922, 3287)
year_labs <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")

png("for_paper/histogram showing singletons.png", height = 400, width = 540)
hist(dataset$dat_time, breaks =  thirty_day_breaks, ylab = "Number of cases", main = "", xlab = "Year",
     col = "#333333", xaxt = "n", ylim = c(0,25))
axis(1, at = year_breaks, labels = year_labs)
hist(singles$dat_time, breaks =  thirty_day_breaks, col = "orange", add = T, labs = "")
dev.off()


##############
## code to plot the clusters on a map

library(sf)
library(tmap)
library(spData)

study_regs <- c("Mtwara", "Lindi")

district_shp <- st_read("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/gis/TZ_District_2012_pop.shp") # Shape files for the districts of Tanzania
plot(district_shp)

study_dis <- subset(district_shp, district_shp$Region_Nam %in% study_regs)
plot(study_dis)
plot(study_dis$geometry)

study_dis_geom <- study_dis$geometry
class(study_dis_geom)
st_crs(study_dis_geom)

study_map <- st_transform(study_dis_geom, 32737)
st_crs(study_map)
# map is in WGS 84 (lat/long)
# coordinates are UTM easting/northing

case_points <- st_as_sf(dataset, coords = c("dat_geo1", "dat_geo2"), crs = 32737)
plot(case_points$geometry)
st_crs(case_points$geometry)

plot(study_map)
plot(case_points, add = T)

xy <- c(dat_geo1, dat_geo2)
(sites <- st_as_sf(dataset, coords = c("dat_geo1", "dat_geo2"),
                   crs = 32737, agr = "constant"))

ggplot(data = study_map) +
  geom_sf() +
  geom_sf(data = sites, size = 1, shape = 23, fill = "darkred")

######



spdf <- SpatialPointsDataFrame(coords = xy, data = dataset,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
levels(spdf$Species)

spdf$cluster_memb <- droplevels(spdf$cluster_memb)

class(case_points)

clusters_vect <-cs_1$cluster_memb
class(clusters_vect)

cluster_df <- cbind(SE_Tanz, res_1_df)
colnames(cluster_df)
cluster_df$cluster_memb <- as.factor(cluster_df$cluster_memb)

cluster_df <- cluster_df[which(cluster_df$cluster_memb %in% clusters_vect),]

xy <- cluster_df[,c("Longitude", "Latitude")]
spdf <- SpatialPointsDataFrame(coords = xy, data = cluster_df,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
levels(spdf$Species)

spdf$cluster_memb <- droplevels(spdf$cluster_memb)

n_cols <-length(unique(spdf$cluster_memb)) # set how many colours we need if colouring by cluster

# cols = c("red", "blue", "yellow") # use this if colouring by species
# use below if colouring by cluster
cols = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_cols)

#
pts = c(17,19, 15) # setting these to use for species

plot(study_dis)
plot(spdf, add = T, pch = pts[spdf$Species], jitter = T, col = cols[spdf$cluster_memb]) #work out how to change point shape and colour by cluster.

##################################################
### This map shows all the clusters of 2 or more.

### would be good to also plot locations of singletons


singles_xy <- vd_singles[,c("Longitude", "Latitude")]
singles_spdf <- SpatialPointsDataFrame(coords = singles_xy, data = vd_singles,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#
pts = c(17,19, 15) # setting these to use for species

plot(study_dis)
plot(singles_spdf, add = T, pch = pts[spdf$Species], jitter = T)#,
#col = cols[spdf$cluster_memb]) #work out how to change point shape and colour by cluster.

# And maybe a plot of the clusters of >=3

clusters_vect_trios <- unlist(cs_1[which(cs_1$total >2),"cluster_memb"])
class(clusters_vect_trios)

cluster_df_trios <- cluster_df[which(cluster_df$cluster_memb %in% clusters_vect_trios),]

xy_trios <- cluster_df_trios[,c("Longitude", "Latitude")]
spdf_trios <- SpatialPointsDataFrame(coords = xy_trios, data = cluster_df_trios,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
levels(spdf_trios$Species)

spdf_trios$cluster_memb <- droplevels(spdf_trios$cluster_memb)

n_cols_trios <-length(unique(spdf_trios$cluster_memb)) # set how many colours we need if colouring by cluster

# cols = c("red", "blue", "yellow") # use this if colouring by species
# use below if colouring by cluster
cols_trios = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_cols_trios)
cols_trios = colorRampPalette(RColorBrewer::brewer.pal(11, "Paired"))(n_cols_trios)


#
pts = c(17,19, 15) # setting these to use for species

plot(study_dis)
plot(spdf_trios, add = T, pch = pts[spdf_trios$Species], jitter = T,
     col = cols_trios[spdf_trios$cluster_memb]) #work out how to change point shape and colour by cluster.









