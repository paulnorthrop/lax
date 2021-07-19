# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("evd", quietly = TRUE)) {
  library(evd)

  # evd::fgev

  # An example from the evd::fgev documentation
  set.seed(4082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  temp <- M1
  adj_fgev <- alogLik(M1)
  class(temp) <- "evd_fgev"
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(adj_fgev))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # No covariates

  M2 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  temp <- M2
  adj_fgev <- alogLik(M2)
  class(temp) <- "evd_fgev"
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("evd::fgev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M2), logLik(adj_fgev))
  })

  # evd::fpot

  # Check whether logLik.pot exists (from package POT)
  # If it does then reverse the class of M1 so that "evd" is first, not "pot"
  all_logLik_methods <- utils::methods(logLik)
  if ("logLik.pot" %in% all_logLik_methods) {
    pot_to_evd <- TRUE
  } else {
    pot_to_evd <- FALSE
  }

  # An example from the evd::fpot documentation
  set.seed(4082019)
  uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)

  # model = "gpd"
  M1 <- evd::fpot(uvdata, 1)
  adj_fpot <- alogLik(M1)
  temp <- M1
  class(temp) <- "evd_fpot"
  if (pot_to_evd) {
    class(M1) <- rev(class(M1))
  }
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fpot, gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("evd::fpot: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(adj_fpot))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fpot, gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # model = "pp"
  M1 <- evd::fpot(uvdata, 1, model = "pp")
  temp <- M1
  adj_fpot <- alogLik(M1)
  class(temp) <- "evd_fpot"
  if (pot_to_evd) {
    class(M1) <- rev(class(M1))
  }
  # Note: evd::logLik.evd returns non-standard attributes (no nobs)
  # Therefore, use expect_equivalent(), rather than expect_equal()
  test_that("evd::fpot, pp: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("evd::fpot: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(M1), logLik(adj_fpot))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("evd::fpot, pp: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}
