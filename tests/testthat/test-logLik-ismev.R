context("logLik, mev package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("mev", quietly = TRUE)) {
  library(mev)

  # mev::gev.fit

  # The example from the mev::gev.fit documentation
  gev_fit <- mev::fit.gev(revdbayes::portpirie, show = FALSE)
  temp <- gev_fit
  adj_gev_fit <- alogLik(gev_fit)
  class(temp) <- "mev_gev"
  test_that("mev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::gev_fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(adj_gev_fit))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("mev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}
