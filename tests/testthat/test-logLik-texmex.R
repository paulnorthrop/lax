context("logLik, texmex package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

# We need the texmex package, and ismev for the fremantle and rain datasets
got_texmex <- requireNamespace("texmex", quietly = TRUE)
got_ismev <- requireNamespace("ismev", quietly = TRUE)

# Examples from the texmex::evm documentation

if (got_texmex) {
  library(texmex)

  # texmex::evm, GEV
  mod <- texmex::evm(SeaLevel, data = portpirie, family = gev)
  temp <- mod
  adj_mod <- alogLik(mod)
  class(temp) <- "texmex_evmOpt"

  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(mod), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(mod), logLik(adj_mod))
  })
  # Check logLik.extRemes_gev, GEV: trivially correct
  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

}

