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

