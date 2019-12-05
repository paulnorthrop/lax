context("logLik, mev package")

# Check that logLik(object) and logLik(logLikVec(object)) agree

if (requireNamespace("mev", quietly = TRUE)) {
  library(mev)

  # mev::fit.gev

  # The example from the mev::gev.fit documentation
  gev_fit <- mev::fit.gev(revdbayes::portpirie, show = FALSE)
  temp <- gev_fit
  adj_gev_fit <- alogLik(gev_fit)
  class(temp) <- "laxmev_gev"
  # Use methods internal to lax, until new mev hits CRAN
  class(gev_fit) <- "laxmev_gev"
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gev_fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gev_fit), logLik(adj_gev_fit))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("mev::fit.gev: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # mev::fit.gpd

  # An example from the mev::fit.gpd documentation
  data(eskrain)
  gpd_mev <- fit.gpd(eskrain, threshold = 35, method = 'Grimshaw')
  temp <- gpd_mev
  adj_gpd_mev <- alogLik(gpd_mev)
  class(temp) <- "laxmev_gpd"
  # Use methods internal to lax, until new mev hits CRAN
  class(gpd_mev) <- "laxmev_gpd"
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gpd_mev), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(gpd_mev), logLik(adj_gpd_mev))
  })
  # Check logLik.evd_fgev: trivially correct
  test_that("mev::fit.gpd: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # mev::fit.pp

  # Use simulated data
  set.seed(1112019)
  x <- revdbayes::rgp(365 * 10, loc = 0, scale = 1, shape = 0.1)
  pfit <- mev::fit.pp(x, threshold = 1, npp = 365)
  # (To do: delete the next two lines after new mev hits CRAN)
  pfit$xdat <- x
  pfit$npp <- 365
  temp <- pfit
  adj_pfit <- alogLik(pfit)
  class(temp) <- "laxmev_pp"
  # Use methods internal to lax, until new mev hits CRAN
  class(pfit) <- "laxmev_pp"
  test_that("mev::fit.pp: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(pfit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.pp: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(pfit), logLik(adj_pfit))
  })
  # Check logLik.gev.fit: trivially correct
  test_that("mev::fit.pp: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # Check that if only threshold exceedances are supplied to mev::fit.pp()
  # then an error is thrown
  data(eskrain)
  pp_mle <- fit.pp(eskrain, threshold = 30, np = 6201)
  adj_pp_mle <- try(alogLik(pp_mle), silent = TRUE)
  test_that("mev::fit.pp: only exceedances throws an error", {
    testthat::expect_identical(class(adj_pp_mle), "try-error")
  })


  # ismev::fit.rlarg

  # An example from the mev::fit.rlarg documentation
  set.seed(31102019)
  xdat <- rrlarg(n = 10, loc = 0, scale = 1, shape = 0.1, r = 4)
  rfit <- fit.rlarg(xdat)
  temp <- rfit
  adj_rfit <- alogLik(rfit)
  class(temp) <- "laxmev_rlarg"
  # logLik.mev_rlarg() produces an object with an extra attribute names
  # equal to "scale".  logLik.mev_rlarg() and logLik.laxmev_rlarg() also
  # disagree on nobs. Therefore, use equivalent in the first two tests
  # Use methods internal to lax, until new mev hits CRAN
  class(rfit) <- "laxmev_rlarg"
  test_that("mev::fit.rlarg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(rfit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("mev::fit.rlarg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(rfit), logLik(adj_rfit))
  })
  # Check logLik.rlarg.fit: trivially correct
  test_that("mev::fit.rlarg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # mev::fit.egp

  # An example from the mev::fit.egp documentation
  if (requireNamespace("evd", quietly = TRUE)) {
    models <- c("egp1", "egp2", "egp3")
    set.seed(7102019)
    xdat <- revdbayes::rgp(n = 100, loc = 0, scale = 1, shape = 0.5)
    for (i in 1:3) {
      fitted <- fit.egp(xdat = xdat, thresh = 1, model = models[i],
                        show = FALSE)
      temp <- fitted
      adj_fitted <- alogLik(fitted)
      class(temp) <- "laxmev_egp"
      # Use methods internal to lax, until new mev hits CRAN
      class(fitted) <- "laxmev_egp"
      # logLik.mev_egp() produces an object with an extra attribute names
      # equal to "scale".  Therefore, use equivalent in the first two tests
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equivalent(logLik(fitted), logLik(logLikVec(temp)))
      })
      # Check that alogLik also returned the correct maximised log-likelihood
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equivalent(logLik(fitted), logLik(adj_fitted))
      })
      # Check logLik.evd_fgev: trivially correct (up to naming)
      test_that("mev::fit.egp: logLik() vs. logLik(logLikVec)", {
        testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
      })
    }
  }
}
