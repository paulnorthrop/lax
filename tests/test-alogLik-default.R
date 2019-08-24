context("alogLik")

# Check that a user-supplied class would work
# Base this on gev_refit

large <- gev_refit(ow$temp, ow, mul = 4, sigl = 4, shl = 4, show = FALSE,
                   method = "BFGS")

# Adjust using the in-built methods, including alogLik.gev.fit()
adj_large <- alogLik(large, cluster = ow$year, cadjust = FALSE)
res <- summary(adj_large)

# Set the class to "ismev_gev" so that all methods apart from alogLik work:
# logLikvec.ismev_gev(), nobs.ismev_gev(), vcov.ismev_gev(), coef.ismev_gev()
temp <- large
class(temp) <- "ismev_gev"

# Now use alogLik.default
adj_test <- alogLik(temp, cluster = ow$year, cadjust = FALSE)
res_test <- summary(adj_test)

test_that("alogLik.defualt works", {
  testthat::expect_equal(res, res_test)
})
