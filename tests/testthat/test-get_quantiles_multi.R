test_that("output is as expected", {
  set.seed(123)
  res_df <- as.data.frame(matrix(ncol = 3, nrow = 3))
  colnames(res_df) <- c("trans_type", "threshold_temporal", "proportion_sim_temporal")
  res_df[,"trans_type"] <- c("g1g1","mixed","g2g2")
  res_df[,"threshold_temporal"] <- c(20.1466871766374, 20.137613324821, 20.1713057402521)
  res_df[,"proportion_sim_temporal"] <- c(0.39048076901592, 0.468867914637129, 0.140651316346951)

  set.seed(123)

  expect_equal(get_quantiles_multi(d_type = "temporal", distrib = "lognormal",
                                   obs = c(500, 300), rr = c(0.5, 0.8),
                                   params = c(5, 1, 5, 1, 5, 1),
                                   n = 10000000, q = 0.95, assort_mix = 1), res_df)
})



