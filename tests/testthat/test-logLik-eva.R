# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("eva", quietly = TRUE)) {
  library(eva)

  # eva::gpdFit

  # An example from the eva::gpdFit documentation
  set.seed(7)
  x <- eva::rgpd(2000, loc = 0, scale = 2, shape = 0.2)
  mle_fit <- eva::gpdFit(x, threshold = 4, method = "mle")
  temp <- mle_fit
  adj_mle_fit <- alogLik(mle_fit)
  class(temp) <- "laxeva_gpd"
  # logLik.gpdFit() and logLik.laxeva_gpd() disagree on nobs:
  # logLik.gpdFit() gives the number of raw observations
  # loglik.laxeva_gpd() gives the number of threshold excesses
  # Therefore, use equivalent in the first two tests
  test_that("eva::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mle_fit), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("eva::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(mle_fit), logLik(adj_mle_fit),
                           ignore_attr = TRUE)
  })
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gpdFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # Check that the information argument to eva:gpdFit() makes no difference
  mle_fit_o <- eva::gpdFit(x, threshold = 4, method = "mle",
                           information = "observed")
  mle_fit_e <- eva::gpdFit(x, threshold = 4, method = "mle",
                           information = "expected")
  adj_mle_fit_o <- alogLik(mle_fit_o)
  adj_mle_fit_e <- alogLik(mle_fit_e)
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gpdFit: observed vs expected information", {
    testthat::expect_equal(summary(adj_mle_fit_o), summary(adj_mle_fit_e))
  })

  # Another example from the eva::gpdFit documentation
  # A linear trend in the scale parameter
  set.seed(7)
  n <- 300
  x2 <- eva::rgpd(n, loc = 0, scale = 1 + 1:n / 200, shape = 0)
  covs <- as.data.frame(seq(1, n, 1))
  names(covs) <- c("Trend1")
  result1 <- eva::gpdFit(x2, threshold = 0, scalevars = covs,
                         scaleform = ~ Trend1)
  adj_result1 <- alogLik(result1)
  temp <- result1
  adj_result1 <- alogLik(result1)
  class(temp) <- "laxeva_gpd"
  # logLik.gpdFit() and logLik.laxeva_gpd() disagree on nobs:
  # logLik.gpdFit() gives the number of raw observations
  # loglik.laxeva_gpd() gives the number of threshold excesses
  # Therefore, use equivalent in the first two tests
  test_that("eva::gpdFit inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(logLikVec(temp)),
                           ignore_attr = TRUE)
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("eva::gpdFit inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(adj_result1),
                           ignore_attr = TRUE)
  })
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gpdFit inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # eva::gevrFit

  # An example from the eva::gevr documentation
  set.seed(7)
  x1 <- eva::rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
  result1 <- eva::gevrFit(x1, method = "mle")
  temp <- result1
  adj_result1 <- alogLik(result1)
  class(temp) <- "laxeva_rlarg"
  # logLik.gpdFit() and logLik.laxeva_gpd() disagree on nobs:
  # logLik.gpdFit() gives the number of raw observations
  # loglik.laxeva_gpd() gives the number of threshold excesses
  # Therefore, use equivalent in the first two tests
  test_that("eva::gevrFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("eva::gevrFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(adj_result1))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gevrFit: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # Repeat for 2-largest order statistics
  set.seed(7)
  x1 <- eva::rgevr(500, 2, loc = 0.5, scale = 1, shape = 0.3)
  result1 <- eva::gevrFit(x1, method = "mle")
  temp <- result1
  adj_result1 <- alogLik(result1)
  class(temp) <- "laxeva_rlarg"
  test_that("eva::gevrFit, r = 2: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("eva::gevrFit, r = 2: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result1), logLik(adj_result1))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gevrFit, r = 2: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # Check that the information argument to eva:gevrFit() makes no difference
  result_o <- eva::gevrFit(x1, method = "mle", information = "observed")
  result_e <- eva::gevrFit(x1, method = "mle", information = "expected")
  adj_result_o <- alogLik(result_o)
  adj_result_e <- alogLik(result_e)
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gevrFit: observed vs expected information", {
    testthat::expect_equal(summary(adj_result_o), summary(adj_result_e))
  })

  # Another example from the eva::gevrFit documentation
  # A linear trend in the location and scale parameter
  n <- 100
  r <- 10
  x2 <- eva::rgevr(n, r, loc = 100 + 1:n / 50,  scale = 1 + 1:n / 300,
                   shape = 0)
  covs <- as.data.frame(seq(1, n, 1))
  names(covs) <- c("Trend1")
  # Create some unrelated covariates
  covs$Trend2 <- rnorm(n)
  covs$Trend3 <- 30 * runif(n)
  result2 <- eva::gevrFit(data = x2, method = "mle", locvars = covs,
                          locform = ~ Trend1 + Trend2*Trend3,
                          scalevars = covs, scaleform = ~ Trend1)
  temp <- result2
  adj_result2 <- alogLik(result2)
  class(temp) <- "laxeva_rlarg"
  test_that("eva::gevrFit, inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result2), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("eva::gevrFit, inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(result2), logLik(adj_result2))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("eva::gevrFit, inc. covariates: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

}
