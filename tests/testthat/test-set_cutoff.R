# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
#
#
# test_that("inputs are correct format", {
#
# })


test_that("dataframe is produced", {
  df <- as.data.frame(matrix(ncol = 3, nrow = 3))
  colnames(df) <- c("trans_type", "threshold_temporal", "proportion_sim_temporal")
  df[,"trans_type"] <- c("g1g1", "mixed", "g2g2")
  df[,"threshold_temporal"] <- c(3,4,5)
  df[,"proportion_sim_temporal"] <- NA
  expect_equal(set_cutoff("temporal", 3,4,5), df)
})


