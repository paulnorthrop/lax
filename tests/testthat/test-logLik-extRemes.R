# Check that logLik(object) and logLik(logLikVec(object)) agree

# We need the extRemes and distillery packages
got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
got_distillery <- requireNamespace("distillery", quietly = TRUE)

# Examples from the extRemes::fevd documentation

if (got_extRemes & got_distillery) {
  library(extRemes)
  library(distillery)
  data(PORTw)

  ##### extRemes::fevd, GEV

  # An example from the extRemes::fgev documentation
  fit0 <- fevd(TMX1, PORTw, units = "deg C", use.phi = TRUE)
  temp <- fit0
  adj_fit0 <- alogLik(fit0)
  class(temp) <- "extRemes_gev"
  test_that("extRemes::fevd, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit0), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit0), logLik(adj_fit0))
  })
  # Check logLik.extRemes_gev, GEV: trivially correct
  test_that("extRemes::fevd, GEV: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, GEV regression in location

  fit1 <- fevd(TMX1, PORTw, location.fun = ~STDTMAX, use.phi = TRUE)
  temp <- fit1
  adj_fit1 <- alogLik(fit1)
  class(temp) <- "extRemes_gev"

  test_that("extRemes::fevd, reg, loc: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, reg, loc: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(adj_fit1))
  })
  # Check logLik.extRemes_gev, GEV: trivially correct
  test_that("extRemes::fevd, reg, loc: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, GEV regression in log(scale)

  fit1 <- fevd(TMX1, PORTw, scale.fun = ~STDTMAX, use.phi = TRUE)
  temp <- fit1
  adj_fit1 <- alogLik(fit1)
  class(temp) <- "extRemes_gev"

  test_that("extRemes::fevd, reg, phi: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, reg, phi: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(adj_fit1))
  })
  # Check logLik.extRemes_gev, GEV: trivially correct
  test_that("extRemes::fevd, reg, phi: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # Repeat for use.phi = FALSE, regression in scale

  fit1 <- fevd(TMX1, PORTw, scale.fun = ~STDTMAX, use.phi = FALSE)
  temp <- fit1
  adj_fit1 <- alogLik(fit1)
  class(temp) <- "extRemes_gev"

  test_that("extRemes::fgev, reg, sigma: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fgev, reg, sigma: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(adj_fit1))
  })
  # Check logLik.extRemes_gev, reg: trivially correct
  test_that("extRemes::fevd, reg, sigma: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  ##### extRemes::fevd, GP

  # An example from the extRemes::fevd documentation
  data(damage)
  fit1 <- fevd(Dam, damage, threshold = 6, type = "GP",
               time.units = "2.05/year")
  temp <- fit1
  adj_fit1 <- alogLik(fit1)
  class(temp) <- "extRemes_gp"

  test_that("extRemes::fevd, GP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, GP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit1), logLik(adj_fit1))
  })
  # Check logLik.extRemes_gp, GEV: trivially correct
  test_that("extRemes::fevd, GP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, GP regression, non-constant threshold

  data(Fort)
  fit <- fevd(Prec, Fort, threshold=0.475,
              threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
              type = "GP")
  temp <- fit
  adj_fit <- alogLik(fit)
  class(temp) <- "extRemes_gp"

  test_that("extRemes::fevd, GP reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, GP reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(adj_fit))
  })
  # Check logLik.extRemes_gp, GP reg: trivially correct
  test_that("extRemes::fevd, reg: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, PP model

  fit <- fevd(Prec, Fort, threshold = 0.475, type = "PP", units = "inches")
  adj_fit <- alogLik(fit)
  temp <- fit
  class(temp) <- "extRemes_pp"

  test_that("extRemes::fevd, PP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, PP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(adj_fit))
  })
  # Check logLik.extRemes_gp, PP: trivially correct
  test_that("extRemes::fevd, PP: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, PP model, non-constant threshold

  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
              type = "PP")
  adj_fit <- alogLik(fit)
  temp <- fit
  class(temp) <- "extRemes_pp"

  test_that("extRemes::fevd, PP, var u: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(logLikVec(temp)))
  })
  # Check that alogLik also returned the correct maximised log-likelihood
  test_that("extRemes::fevd, PP, var u: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equivalent(logLik(fit), logLik(adj_fit))
  })
  # Check logLik.extRemes_pp, PP: trivially correct
  test_that("extRemes::fevd, PP, var u: logLik() vs. logLik(logLikVec)", {
    testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
  })

  # extRemes::fevd, PP regression
  # Use example from the ismev package (from chapter 7 of Coles (2001))
  # Code from demo ismev::wooster.temps

  if (requireNamespace("ismev", quietly = TRUE)) {
    library(ismev)
    data(wooster)
    x <- seq(along = wooster)
    usin <- function(x, a, b, d) {
      a + b * sin(((x - d) * 2 * pi) / 365.25)
    }
    wu <- usin(x, -30, 25, -75)
    ydat <- cbind(sin(2 * pi * x / 365.25), cos(2 * pi *x / 365.25))
    df_wooster <- data.frame(y = -wooster, sin = ydat[, 1], cos = ydat[, 2])
    fitPP <- fevd(y, df_wooster, threshold = wu, location.fun = ~ sin + cos,
                  scale.fun = ~ sin + cos, type = "PP", use.phi = TRUE)
    adj_fit <- alogLik(fitPP)
    temp <- fitPP
    class(temp) <- "extRemes_pp"

    test_that("extRemes::fevd, PP reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equivalent(logLik(fitPP), logLik(logLikVec(temp)))
    })
    # Check that alogLik also returned the correct maximised log-likelihood
    test_that("extRemes::fevd, PP reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equivalent(logLik(fitPP), logLik(adj_fit))
    })
    # Check logLik.extRemes_gp, PP: trivially correct
    test_that("extRemes::fevd, PP reg: logLik() vs. logLik(logLikVec)", {
      testthat::expect_equal(logLik(temp), logLik(logLikVec(temp)))
    })
  }
}

