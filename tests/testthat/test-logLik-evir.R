context("logLik, evir package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

# evir::gev

if (requireNamespace("evir", quietly = TRUE)) {
  library(evir)

  # An example from the evir::gev documentation
  data(bmw)
  out <- evir::gev(bmw, "month")
  temp <- out
  adj_out <- alogLik(out)
  class(temp) <- "evir_gev"
  test_that("evir::gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returnd the correct maximised log-likelihood
  test_that("evd::gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(adj_out))
  })
  # Check logLik.evir_gev: trivially correct
  test_that("evir::gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # An example from the evir::gpd documentation
  data(danish)
  out <- evir::gpd(danish, 10)
  temp <- out
  adj_out <- alogLik(out)
  class(temp) <- "evir_gpd"
  test_that("evir::gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returnd the correct maximised log-likelihood
  test_that("evd::gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(adj_out))
  })
  # Check logLik.evir_gpd: trivially correct
  test_that("evir::gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # An example from the evir::pot documentation
  # We use lax::re_pot() to return the input data
  out <- re_pot(danish, 10)
  temp <- out
  adj_out <- alogLik(out)
  class(temp) <- "evir_pot"
  test_that("evir::pot logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returnd the correct maximised log-likelihood
  test_that("evd::pot: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(out), logLik(adj_out))
  })
  # Check logLik.evir_gpd: trivially correct
  test_that("evir::pot: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

}

