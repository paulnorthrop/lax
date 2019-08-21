context("nobs, evd package")

# Check that nobs.evd bechaves correctly

if (requireNamespace("evd", quietly = TRUE)) {
  library(evd)

  # evd::fgev

  # An example from the evd::fgev documentation
  set.seed(4082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  test_that("evd::fgev, nobs.gev() vs. length(response)", {
    testthat::expect_equal(nobs(M1), length(uvdata))
  })

  # evd::fpot

  # An example from the evd::fpot documentation
  set.seed(4082019)
  uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
  u <- 1
  M2 <- evd::fpot(uvdata, u)
  test_that("evd::fpot, nobs.evd() vs. sum(response > u)", {
    testthat::expect_equivalent(nobs(M2), sum(uvdata > u))
  })

  # evd::fextreme

  # An example from the evd::fextreme documentation
  uvdata <- evd::rextreme(100, qnorm, mean = 0.56, mlen = 365)
  M3 <- evd::fextreme(uvdata, list(mean = 0, sd = 1), distn = "norm", mlen = 365)
  test_that("evd::fextreme, nobs.evd() vs. length(response)", {
    testthat::expect_equivalent(nobs(M3), length(uvdata))
  })
}
