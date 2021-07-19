# Check that adjusted fits from ismev::rlarg.fit(), or lax::rlarg_refit()
# using only the largest observations in a period agree with those from
# ismev::gev.fit() and lax::gev_refit

if (requireNamespace("ismev", quietly = TRUE)) {
  library(ismev)

  # Based on the example from the ismev::rlarg.fit documentation
  # Use revdbayes::venice to avoid ambiguity
  # It's the same as ismev::venice but years are in row names not column 1
  gfit <- gev.fit(revdbayes::venice[, 1], muinit = 111.1, siginit = 17.18,
                  shinit = -0.0767, show = FALSE)
  rfit <- rlarg.fit(revdbayes::venice[, 1, drop = FALSE], muinit = 111.1,
                    siginit = 17.18, shinit = -0.0767, show = FALSE)
  adj_gfit <- alogLik(gfit)
  adj_rfit <- alogLik(rfit)
  s_gfit <- summary(adj_gfit)
  s_rfit <- summary(adj_rfit)
  test_that("ismev, gev vs rlarg", {
    testthat::expect_equal(s_gfit, s_rfit)
  })

  # An example from chapter 6 of Coles (2001)
  data(fremantle)
  xdat <- fremantle[, "SeaLevel"]
  # Set year 1897 to 1 for consistency with page 113 of Coles (2001)
  ydat <- cbind(fremantle[, "Year"] - 1896, fremantle[, "SOI"])
  #
  gev_fit <- gev_refit(xdat, ydat, mul = 1:2, show = FALSE)
  adj_gev_fit <- alogLik(gev_fit)
  s_gev_fit <- summary(adj_gev_fit)
  #
  xdat <- fremantle[, "SeaLevel", drop = FALSE]
  rlarg_fit <- rlarg_refit(xdat, ydat = ydat, mul = 1:2, show = FALSE)
  adj_rlarg_fit <- alogLik(gev_fit)
  s_rlarg_fit <- summary(adj_gev_fit)
  test_that("ismev, gev vs rlarg, regression", {
    testthat::expect_equal(s_gev_fit, s_rlarg_fit)
  })

}
