context("vcov")

# Check that vcov methods agree

# evir

if (requireNamespace("evir", quietly = TRUE)) {
  library(evir)

  # An example from the evir::gev documentation
  data(bmw)
  out <- evir::gev(bmw, "month")
  out2 <- out
  class(out2) <- "evir_gev"
  test_that("evir::gev: vcov.gev vs vcov.evir_gev", {
    testthat::expect_equal(vcov(out), vcov(out2))
  })
}

