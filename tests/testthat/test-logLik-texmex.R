# Check that logLik(object) and logLik(logLikVec(object)) agree

# We need the texmex package, and ismev for the fremantle and rain datasets
got_texmex <- requireNamespace("texmex", quietly = TRUE)
got_ismev <- requireNamespace("ismev", quietly = TRUE)

# Examples from the texmex::evm documentation

if (got_texmex) {
  library(texmex)

  # texmex::evm, GEV
  mod <- texmex::evm(SeaLevel, data = texmex::portpirie, family = gev)
  temp <- mod
  adj_mod <- alogLik(mod)
  class(temp) <- "texmex_evmOpt"

  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(adj_mod), ignore_attr = TRUE)
  })
  # Check logLik.texmex_evmOpt, GEV: trivially correct
  test_that("texmex::evm, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # texmex::evm, EGP3
  mod <- texmex::evm(rain, th = 30, family = egp3)
  temp <- mod
  adj_mod <- alogLik(mod)
  class(temp) <- "texmex_evmOpt"

  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(adj_mod), ignore_attr = TRUE)
  })
  # Check logLik.texmex_evmOpt, EGP3: trivially correct
  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # GP
  mod <- texmex::evm(rain, th = 30)
  temp <- mod
  adj_mod <- alogLik(mod)
  class(temp) <- "texmex_evmOpt"

  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(adj_mod), ignore_attr = TRUE)
  })
  # Check logLik.texmex_evmOpt, EGP3: trivially correct
  test_that("texmex::evm, EGP3: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # GEV regression
  # An example from page 113 of Coles (2001)
  if (got_ismev) {
    library(ismev)
    data(fremantle)
    new_fremantle <- fremantle
    # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
    new_fremantle[, "Year"] <- new_fremantle[, "Year"] - 1896
    mod <- texmex::evm(y = SeaLevel, data = new_fremantle, family = gev,
                       mu = ~ Year + SOI)
    temp <- mod
    adj_mod <- alogLik(mod)
    class(temp) <- "texmex_evmOpt"

    test_that("texmex::evm, GEV, reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equal(logLik(mod), logLik(logLikVec(temp)),
                             ignore_attr = TRUE)
    })
    # Check that alogLik also returned the correct maximised log-likelihood
    test_that("texmex::evm, GEV, reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equal(logLik(mod), logLik(adj_mod), ignore_attr = TRUE)
    })
    # Check logLik.texmex_evmOpt, EGP3: trivially correct
    test_that("texmex::evm, GEV, reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
    })
  }

  # GP regression
  # An example from page 119 of Coles (2001)
  n_rain <- length(rain)
  rain_df <- data.frame(rain = rain, time = 1:n_rain / n_rain)
  mod <- texmex::evm(y = rain, data = rain_df, family = gpd, th = 30,
                     phi = ~ time)
  temp <- mod
  adj_mod <- alogLik(mod)
  class(temp) <- "texmex_evmOpt"

  test_that("texmex::evm, GP, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("texmex::evm, GP, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mod), logLik(adj_mod), ignore_attr = TRUE)
  })
  # Check logLik.texmex_evmOpt, EGP3: trivially correct
  test_that("texmex::evm, GP, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })
}

