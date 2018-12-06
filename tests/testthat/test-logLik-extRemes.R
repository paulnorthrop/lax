context("logLik, extRemes package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

# We need the extRemes and distillery packages
got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
got_distillery <- requireNamespace("distillery", quietly = TRUE)

# Examples from the extRemes::fevd documentation

if (got_extRemes & got_distillery) {
  library(extRemes)
  library(distillery)
  data(PORTw)

  # extRemes::fevd, GEV

  # An example from the extRemes::fgev documentation
  fit0 <- fevd(TMX1, PORTw, units = "deg C", use.phi = TRUE)
  temp <- fit0
  class(temp) <- "extRemes_gev"
  test_that("extRemes::fgev, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit0), logLik(logLikVec(temp)))
  })
  # Check logLik.extRemes_gev, GEV: trivially correct
  test_that("extRemes::fevd, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}

