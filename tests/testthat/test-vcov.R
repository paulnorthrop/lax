context("vcov")

# Check that vcov methods agree

# ---------------------------------- evir ------------------------------------#

if (requireNamespace("evir", quietly = TRUE)) {
  library(evir)

  # An example from the evir::gev documentation
  data(bmw)
  out <- evir::gev(bmw, "month")
  temp <- out
  class(temp) <- "evir_gev"
  test_that("evir::gev: vcov.gev vs vcov.evir_gev", {
    testthat::expect_equal(vcov(out), vcov(temp))
  })

  # An example from the evir::gpd documentation
  data(danish)
  out <- evir::gpd(danish, 10)
  temp <- out
  class(temp) <- "evir_gpd"
  test_that("evir::gpd: vcov.gpd vs vcov.evir_gpd", {
    testthat::expect_equal(vcov(out), vcov(temp))
  })

  # An example from the evir::pot documentation
  out <- evir::pot(danish, 10)
  temp <- out
  class(temp) <- "evir_pot"
  test_that("evir::pot: vcov.potd vs vcov.evir_pot", {
    testthat::expect_equal(vcov(out), vcov(temp))
  })
}

# ---------------------------------- evd ------------------------------------#

if (requireNamespace("evir", quietly = TRUE)) {
  library(evd)

  # An example from the evd::fgev documentation
  set.seed(3082019)
  uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
  M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
  temp <- M1
  class(temp) <- "evd_fgev"
  test_that("evd::fgev: vcov.gev vs vcov.evd_gev", {
    # column names differ
    testthat::expect_equivalent(vcov(M1), vcov(temp))
  })
}

# -------------------------------- extRemes --------------------------------- #

if (requireNamespace("evir", quietly = TRUE)) {
  library(extRemes)

  fitPORTstdmax <- extRemes::fevd(TMX1, PORTw, scale.fun = ~STDTMAX,
                                  use.phi = TRUE)
  temp <- fitPORTstdmax
  class(temp) <- "extRemes_gev"
  test_that("extRemes::fevd: vcov.fevd vs vcov.extRemes_gev", {
    # column names differ
    testthat::expect_equivalent(vcov(fitPORTstdmax), vcov(temp))
  })
}
