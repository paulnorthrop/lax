# Check that when type = "none" (no adjustment) anova.lax and calculation
# of the LRT by hand give the same answer

got_evd <- requireNamespace("evd", quietly = TRUE)

if (got_evd) {
  library(evd)
  small <- fgev(ow$temp, nsloc = ow$loc)
  adj_small <- alogLik(small, cluster = ow$year)
  tiny <- fgev(ow$temp)
  adj_tiny <- alogLik(tiny, cluster = ow$year)

  from_anova <- anova(adj_small, adj_tiny, type = "none")$ALRTS[2]
  by_hand <- as.numeric(2 * (logLik(adj_small) - logLik(adj_tiny)))

  test_that("anova.lax vs by hand", {
    testthat::expect_equal(from_anova, by_hand)
  })
}
