# Check that nobs.evd behaves correctly

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
  adj_M1 <- alogLik(M1)
  test_that("texmex::evm, nobs.evm_Opt vs. nobs.lax", {
    testthat::expect_equal(nobs(M1), nobs(adj_M1))
  })

  # evd::fpot
  # An example from the evd::fpot documentation
  set.seed(4082019)
  uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
  u <- 1
  M2 <- evd::fpot(uvdata, u)
  test_that("evd::fpot, nobs.evd() vs. sum(response > u)", {
    testthat::expect_equal(nobs(M2), sum(uvdata > u), ignore_attr = TRUE)
  })
  adj_M2 <- alogLik(M2)
  test_that("texmex::evm, nobs.evm_Opt vs. nobs.lax", {
    testthat::expect_equal(nobs(M2), nobs(adj_M2))
  })

  # evd::fextreme
  # An example from the evd::fextreme documentation
  uvdata <- evd::rextreme(100, qnorm, mean = 0.56, mlen = 365)
  M3 <- evd::fextreme(uvdata, list(mean = 0, sd = 1), distn = "norm",
                      mlen = 365)
  test_that("evd::fextreme, nobs.evd() vs. length(response)", {
    testthat::expect_equal(nobs(M3), length(uvdata), ignore_attr = TRUE)
  })
}

# Check that nobs.evmOpt behaves correctly

if (requireNamespace("texmex", quietly = TRUE)) {
  library(texmex)
  # texmex::evm, GEV
  mod <- texmex::evm(SeaLevel, texmex::portpirie, family = gev)
  test_that("texmex::evm, nobs.evmOpt vs. length(response)", {
    testthat::expect_equal(nobs(mod), length(texmex::portpirie$SeaLevel))
  })
  adj_mod <- alogLik(mod)
  test_that("texmex::evm, nobs.evmOpt vs. nobs.lax", {
    testthat::expect_equal(nobs(mod), nobs(adj_mod))
  })
}

# Check that nobs.pot behaves correctly

if (requireNamespace("POT", quietly = TRUE)) {
  library(POT)
  # An example from the POT::fitgpd documentation.
  set.seed(24082019)
  x <- POT::rgpd(200, 1, 2, 0.25)
  u <- 1.5
  mod <- POT::fitgpd(x, u, "mle")
  test_that("POT::fitgpd, nobs.pot vs. sum(x > u)", {
    testthat::expect_equal(nobs(mod), sum(x > u))
  })
  adj_mod <- alogLik(mod)
  test_that("POT::fitgpd, nobs.pot vs. nobs.lax", {
    testthat::expect_equal(nobs(mod), nobs(adj_mod))
  })
}
