context("logLik")

# Check that logLik(object) and logLik(logLikVec(object)) agree

# evd::fgev

if (requireNamespace("evd", quietly = TRUE)) {
  library(evd)
  # An example from the evd::fgev documentation
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  temp <- M1
  class(temp) <- "evd_fgev"
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}

# evd::fpot

if (requireNamespace("evd", quietly = TRUE)) {
  library(evd)
  # An example from the evd::fpot documentation
  uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)

  # model = "gpd"
  M1 <- evd::fpot(uvdata, 1)
  temp <- M1
  class(temp) <- "evd_fpot"
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # model = "pp"
  M1 <- evd::fpot(uvdata, 1, model = "pp")
  temp <- M1
  class(temp) <- "evd_fpot"
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}

# ismev::gev.fit

if (requireNamespace("ismev", quietly = TRUE)) {
  library(ismev)
  # The example from the ismev::gev.fit documentation
  gev_fit <- ismev::gev.fit(revdbayes::portpirie, show = FALSE)
  temp <- gev_fit
  class(temp) <- "ismev_gev"
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # An example from page 113 of Coles (2001)
  data(fremantle)
  xdat <- fremantle[, "SeaLevel"]
  # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
  ydat <- cbind(fremantle[, 1] - 1896, fremantle[, 3])
  gev_fit <- oogev.fit(xdat, ydat, mul = 1:2, show = FALSE)
  temp <- gev_fit
  class(temp) <- "ismev_gev"
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}
