context("anova")

# Check that when type = "none" (no adjustment) anova.lax and calculation
# of the LRT by hand give the same answer

got_evd <- requireNamespace("evd", quietly = TRUE)

if (got_evd) {
  library(evd)
  y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
  x <- rep(c(1, -1), each = length(y) / 2)
  owfit <- evd::fgev(y, nsloc = x)
  year <- rep(1:(length(y) / 2), 2)
  small <- alogLik(owfit, cluster = year)

  owfit <- evd::fgev(y)
  year <- rep(1:(length(y) / 2), 2)
  tiny <- alogLik(owfit, cluster = year)

  from_anova <- anova(small, tiny, type = "none")$ALRTS[2]
  by_hand <- as.numeric(2 * (logLik(small) - logLik(tiny)))

  test_that("anova.lax vs by hand", {
    testthat::expect_equal(from_anova, by_hand)
  })
}
