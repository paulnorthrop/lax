context("alogLik")

# Check that a user-supplied class would work
# Base this on gev_refit

large <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4, show = FALSE,
                   method = "BFGS")
adj_large <- alogLik(large, cluster = ow$year, cadjust = FALSE)
res <- summary(adj_large)

temp <- large
class(temp) <- "test"
logLikVec.test <- logLikVec.ismev_gev
nobs.test <- nobs.ismev_gev
vcov.test <- vcov.ismev_gev
coef.test <- coef.ismev_gev
adj_test <- alogLik(temp, cluster = ow$year, cadjust = FALSE)
res_test <- summary(adj_large)

test_that("alogLik.defulat works", {
  testthat::expect_equal(res, res_test)
})
