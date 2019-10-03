context("logLik, ismev package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("ismev", quietly = TRUE)) {
  library(ismev)

  # ismev::gev.fit

  # The example from the ismev::gev.fit documentation
  gev_fit <- ismev::gev.fit(revdbayes::portpirie, show = FALSE)
  temp <- gev_fit
  adj_gev_fit <- alogLik(gev_fit)
  class(temp) <- "ismev_gev"
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::gev_fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(adj_gev_fit))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("ismev::gev.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # An example from page 113 of Coles (2001)
  data(fremantle)
  xdat <- fremantle[, "SeaLevel"]
  # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
  ydat <- cbind(fremantle[, 1] - 1896, fremantle[, 3])
  gev_fit <- gev_refit(xdat, ydat, mul = 1:2, show = TRUE)
  temp <- gev_fit
  adj_gev_fit <- alogLik(gev_fit)
  class(temp) <- "ismev_gev"
  test_that("ismev::gev.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::gev_fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(gev_fit), logLik(adj_gev_fit))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("ismev::gev.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # ismev::gpd.fit

  # An example from the ismev::gpd.fit documentation
  data(rain)
  rain_fit <- gpd.fit(rain, 10, show = FALSE)
  temp <- rain_fit
  adj_rain_fit <- alogLik(rain_fit)
  class(temp) <- "ismev_gpd"
  test_that("ismev::gpd.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(rain_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::gpd.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(rain_fit), logLik(adj_rain_fit))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("ismev::gpd.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
  # Continuing to the regression example on page 119 of Coles (2001)
  ydat <- as.matrix((1:length(rain)) / length(rain))
  reg_rain_fit <- gpd_refit(rain, 30, ydat = ydat, sigl = 1, siglink = exp,
                            show = TRUE)
  temp <- reg_rain_fit
  adj_reg_rain_fit <- alogLik(reg_rain_fit)
  class(temp) <- "ismev_gpd"
  test_that("ismev::gpd.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(reg_rain_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::gpd.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(reg_rain_fit), logLik(adj_reg_rain_fit))
  })
  # Check logLik.gpd.fit: trivially correct
  test_that("ismev::gpd.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # ismev::pp.fit

  # An example from the ismev::pp.fit documentation
  data(rain)
  init <- c(40.55755732, 8.99195409, 0.05088103)
  muinit <- init[1]
  siginit <- init[2]
  shinit <- init[3]
  rain_fit <- pp_refit(rain, 10, muinit = muinit, siginit = siginit,
                       shinit = shinit, show = TRUE)
  temp <- rain_fit
  adj_rain_fit <- alogLik(rain_fit)
  class(temp) <- "ismev_pp"
  test_that("ismev::pp.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(rain_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::pp.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(reg_rain_fit), logLik(adj_reg_rain_fit))
  })
  # Check logLik.pp.fit: trivially correct
  test_that("ismev::pp.fit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
  # An example from chapter 7 of Coles (2001).
  # Code from demo ismev::wooster.temps
  data(wooster)
  x <- seq(along = wooster)
    usin <- function(x, a, b, d) {
    a + b * sin(((x - d) * 2 * pi) / 365.25)
  }
  wu <- usin(x, -30, 25, -75)
  ydat <- cbind(sin(2 * pi * x / 365.25), cos(2 * pi *x / 365.25))
  init <- c(-15.3454188, 9.6001844, 28.5493828, 0.5067104, 0.1023488,
            0.5129783, -0.3504231)
  muinit <- init[1:3]
  siginit <- init[4:6]
  shinit <- init[7]
  wooster.pp <- pp_refit(-wooster, threshold = wu, ydat = ydat, mul = 1:2,
                         sigl = 1:2, siglink = exp, method = "BFGS",
                         muinit = muinit, siginit = siginit, shinit = shinit,
                         show = FALSE)
  temp <- wooster.pp
  adj_wooster.pp <- alogLik(wooster.pp)
  class(temp) <- "ismev_pp"
  test_that("ismev::pp.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(wooster.pp), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("ismev::pp.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(wooster.pp), logLik(adj_wooster.pp))
  })
  # Check logLik.pp.fit: trivially correct
  test_that("ismev::pp.fit, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}
