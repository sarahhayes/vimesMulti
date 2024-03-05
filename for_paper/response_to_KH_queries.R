

## In the feedback to the initial paper, KH suggested that she'd like to see:
##  a) the frequency distributions of minimal distances between cases, versus the pdf of the Raleigh distribution,
##     with the cutoffs indicated, to see the extent to which the algorithm does not link cases -
##     this would make the method (and the applicability to the data) more transparent (to me).
##  b) the distribution of distances to cases where the exact GPS is known vs where there is 5-10 km of uncertainty to see how this is affecting inference.


rm(list = ls())
library(tidyverse)
library(extraDistr)

## need to establish the minimal distances between cases
case_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_time_diffs.csv")
case_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/case_dists.csv")
uncert_times <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/dates_uncertainty.csv")
uncert_dists <- read.csv("C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/Vimes/vimes_multi_sim/tests_sh/temp_trash/distance_uncertainty.csv")

dat_dist <- as.matrix(cbind(case_dists$UTM.Easting, case_dists$UTM.Northing))
dat_dist <- fields::rdist(dat_dist)/1000 # convert to km

# want to extract the minimum distances between all the observed cases.
# first make a smaller version of the distance matrix so can see if working

mini_mat <- dat_dist[1:10, 1:10]
mini_mat
apply(mini_mat, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row
apply(mini_mat, 2, FUN=function(x) {min(x[x > 0])}) # this does it y comumn - as matrix is symmetrical they should be the same


min_case_dists <- apply(dat_dist, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row
range(min_case_dists)
# going to round to nearest m
min_case_dists <- round(min_case_dists, 3)
range(min_case_dists)

dist_breaks <- seq(0,36,0.2)
hist(min_case_dists, breaks = dist_breaks, main = "histogram of minimum distances betwen cases", col = "light blue")

## produce pdf of the rayleigh dist

x<-seq(from=0,to=36,length.out=100)

rayleigh_pdf <- drayleigh(x, sigma = 0.87)
curve(drayleigh(x, sigma = 0.87), from = 0, to = 36)

qrayleigh(0.95, sigma = 0.87)
qrayleigh(0.99, sigma = 0.87)
qrayleigh(1, sigma = 0.87)


dev.off()
#pdf("for_paper/histogram_min_dist_plus_rayleigh.pdf", height = 8, width = 8)
png("for_paper/histogram_min_dist_plus_rayleigh.png")
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
hist(min_case_dists, breaks = dist_breaks, main = "", col = "light blue", xlab= "Minimum distance between cases (km)")
par(new = TRUE)
plot(dist_breaks, drayleigh(dist_breaks, sigma = 0.87), type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "", lwd= 1.3)
axis(side=4, at = pretty(range(rayleigh_pdf)))
mtext("Density", side=4, line=3)
abline(v = qrayleigh(0.95, sigma = 0.87), col = "red")
dev.off()

## rpt zoomed in on the 10km area to make easier to see

z <-seq(from=0,to=10,length.out=100)
zoomed_rayleigh_pdf <- drayleigh(z, sigma = 0.87)

min_case_dists_zoomed <- min_case_dists[min_case_dists<10]
zoom_breaks <- seq(0,10, 0.1)

dev.off()
png("for_paper/histogram_min_dist_plus_rayleigh_zoomed.png")
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
hist(min_case_dists_zoomed, breaks = zoom_breaks, main = "", col = "light pink", xlab= "Minimum distance between cases (km)")
par(new = TRUE)
plot(zoom_breaks, drayleigh(zoom_breaks, sigma = 0.87), type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "", lwd= 1.3)
axis(side=4, at = pretty(range(rayleigh_pdf)))
mtext("Density", side=4, line=3)
abline(v = qrayleigh(0.95, sigma = 0.87), col = "blue")
dev.off()

## also do a zoomed in plot with the cutoff values from the best-fitting scenarios.

dev.off()
png("for_paper/histogram_min_dist_plus_rayleigh_scenario_cuts.png")
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
hist(min_case_dists, breaks = dist_breaks, main = "", col = "grey", xlab= "Minimum distance between cases (km)")
par(new = TRUE)
plot(dist_breaks, drayleigh(dist_breaks, sigma = 0.87), type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "", lwd= 1.3)
axis(side=4, at = pretty(range(rayleigh_pdf)))
mtext("Density", side=4, line=3)
mtext("Density", side=4, line=3)
abline(v = c(2.7, 3.2, 3.7), col = "blue")
abline(v = c(4.3, 5.0, 5.8), col = "red")
dev.off()


## The next part is to look at "the distribution of distances to cases where the exact GPS is known vs
## where there is 5-10 km of uncertainty to see how this is affecting inference.

locs <- cbind(case_dists, uncert_dists)
exact_locs <- locs[which(locs$lower_uncert == 0 & locs$upper_uncert == 0),]
uncert_locs <- locs[which(locs$upper_uncert != 0),]
uncert_2k <- locs[which(locs$upper_uncert == 2000),]
uncert_5k <- locs[which(locs$upper_uncert == 5000),]
uncert_10k <- locs[which(locs$upper_uncert == 10000),]


cases_with_no_uncert <- min_case_dists[which(uncert_dists$upper_uncert == 0)]
cases_with_2k_uncert <- min_case_dists[which(uncert_dists$upper_uncert == 2000)]
cases_with_5k_uncert <- min_case_dists[which(uncert_dists$upper_uncert == 5000)]
cases_with_10k_uncert <- min_case_dists[which(uncert_dists$upper_uncert == 10000)]

png("for_paper/histograms of minimum distances based on uncertainty.png")
par(mfrow = c(2,2))
hist(cases_with_no_uncert, breaks = dist_breaks, col = "orange", main = "Exact", xlab = "Minimum distance between cases")
hist(cases_with_2k_uncert, breaks = dist_breaks, col = "yellow", main = "2Km", xlab = "Minimum distance between cases")
hist(cases_with_5k_uncert, breaks = dist_breaks, col = "green", main = "5Km", xlab = "Minimum distance between cases")
hist(cases_with_10k_uncert, breaks = dist_breaks, col = "purple", main = "10Km", xlab = "Minimum distance between cases")
dev.off()

# We could do distance matrices for all of these, but not sure how much this will bring?

dist_exact <- as.matrix(cbind(exact_locs$UTM.Easting, exact_locs$UTM.Northing))
dist_exact <- fields::rdist(dist_exact)/1000 # convert to km
min_dist_exact <- apply(dist_exact, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row

dist_2K <- as.matrix(cbind(uncert_2k$UTM.Easting, uncert_2k$UTM.Northing))
dist_2K <- fields::rdist(dist_2K)/1000 # convert to km
min_dist_2K <- apply(dist_2K, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row

dist_5K <- as.matrix(cbind(uncert_5k$UTM.Easting, uncert_5k$UTM.Northing))
dist_5K <- fields::rdist(dist_5K)/1000 # convert to km
min_dist_5K <- apply(dist_5K, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row

dist_10K <- as.matrix(cbind(uncert_10k$UTM.Easting, uncert_10k$UTM.Northing))
dist_10K <- fields::rdist(dist_10K)/1000 # convert to km
min_dist_10K <- apply(dist_10K, 1, FUN=function(x) {min(x[x > 0])}) # this extracts the minimum value from each row

range(dist_exact)

par(mfrow = c(2,2))
hist(dist_exact, breaks = dist_breaks, col = "orange", main = "Exact")
