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

if (requireNamespace("evd", quietly = TRUE)) {
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

if (requireNamespace("extRemes", quietly = TRUE)) {
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

# -------------------------------- fExtremes --------------------------------- #

if (requireNamespace("fExtremes", quietly = TRUE)) {
  library(fExtremes)

  # An example from the fExtremes::gevFit documentation
  set.seed(4082019)
  x <- gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
  # Fit GEV distribution by maximum likelihood estimation
  fit <- gevFit(x)
  temp <- fit
  class(temp) <- "fExtremes_gev"

  test_that("fEextremes::gevFIT: vcov.fGEVFIT vs vcov.fExtremes_gev", {
    testthat::expect_equal(vcov(fit), vcov(temp))
  })

  # An example from the fExtremes::gpdFit documentation
  # Simulate GP data
  x <- gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
  # Fit GP distribution by maximum likelihood estimation
  fit <- gpdFit(x, u = min(x))
  temp <- fit
  class(temp) <- "fExtremes_gpd"

  test_that("fEextremes::gevFIT: vcov.fGPDFIT vs vcov.fExtremes_gpd", {
    testthat::expect_equal(vcov(fit), vcov(temp))
  })
}

# --------------------------------- texmex ---------------------------------- #

if (requireNamespace("texmex", quietly = TRUE)) {
  library(texmex)

  mod <- evm(rain, th = 30)
  temp <- mod
  class(temp) <- "texmex_evmOpt"
  test_that("texmex::evm: vcov.evmOpt vs vcov.texmex_evmOpt", {
    # column names differ
    testthat::expect_equivalent(vcov(mod), vcov(temp))
  })
}
