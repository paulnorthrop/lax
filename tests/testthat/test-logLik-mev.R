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

  # mev::fit.egp

  # An example from the mev::fit.egp documentation
  if (requireNamespace("evd", quietly = TRUE)) {
    models <- c("egp1", "egp2", "egp3")
    set.seed(7102019)
    xdat <- evd::rgpd(n = 100, loc = 0, scale = 1, shape = 0.5)
    for (i in 1:3) {
      fitted <- fit.egp(xdat = xdat, thresh = 1, model = models[i],
                        show = FALSE)
      temp <- fitted
      adj_fitted <- alogLik(fitted)
      class(temp) <- "mev_egp"
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equivalent(logLik(fitted), logLik(logLikVec(temp)))
      })
      # Check that alogLik also returned the correct maximised log-likelihood
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equivalent(logLik(fitted), logLik(adj_fitted))
      })
      # Check logLik.evd_fgev: trivially correct (up to naming)
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equivalent(logLik(temp), logLik(logLikVec(temp)))
      })
    }
  }
}
