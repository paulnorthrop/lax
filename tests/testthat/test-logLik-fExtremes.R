# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("fExtremes", quietly = TRUE)) {
  library(fExtremes)

  # fExtremes::gevFit

  # An example from the POT::fitgpd documentation.
  set.seed(4082019)
  # An example from the fExtremes::gevFit documentation
  x <- fExtremes::gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
  # Fit GEV distribution by maximum likelihood estimation
  fit <- fExtremes::gevFit(x)
  temp <- fit
  adj_fit <- alogLik(fit)
  class(temp) <- "fExtremes_gev"

  test_that("fExtremes::gevFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(fit), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("fExtremes::gevFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(fit), logLik(adj_fit), ignore_attr = TRUE)
  })
  # Check logLik.POT_fitgpd: trivially correct
  test_that("fExtremes::gevFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # fExtremes::gpdFit

  # An example from the fExtremes::gpdFit documentation
  # Simulate GP data
  x <- fExtremes::gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
  # Fit GP distribution by maximum likelihood estimation
  fit <- fExtremes::gpdFit(x, u = min(x))
  temp <- fit
  adj_fit <- alogLik(fit)
  class(temp) <- "fExtremes_gpd"

  my_tol <- 1e-4

  test_that("fExtremes::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(fit), logLik(logLikVec(temp)),
                           tolerance = my_tol, ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("fExtremes::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(fit), logLik(adj_fit), tolerance = my_tol,
                           ignore_attr = TRUE)
  })
  # Check logLik.POT_fitgpd: trivially correct
  test_that("fExtremes::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

}
