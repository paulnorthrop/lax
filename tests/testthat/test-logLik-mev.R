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
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(adj_gev_fit))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # mev::gpd.fit

  # An example from the mev::fit.gpd documentation
  data(eskrain)
  gpd_mev <- fit.gpd(eskrain, threshold = 35, method = 'Grimshaw')
  temp <- gpd_mev
  adj_gpd_mev <- alogLik(gpd_mev)
  class(temp) <- "mev_gpd"
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gpd_mev), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gpd_mev), logLik(adj_gpd_mev))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}
