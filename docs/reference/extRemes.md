# Loglikelihood adjustment for extRemes fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the function `fevd` in the
`extRemes` package. The model must have been fitted using maximum
likelihood estimation.

## Usage

``` r
# S3 method for class 'fevd'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)
```

## Arguments

- x:

  A fitted model object with certain associated S3 methods. See
  **Details**.

- cluster:

  A vector or factor indicating from which cluster the respective
  log-likelihood contributions from `loglik` originate. The length of
  `cluster` must be consistent with the `estfun` method to be used in
  the estimation of the 'meat' `V` of the sandwich estimator of the
  covariance matrix of the parameters to be passed to
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).
  In most cases, `cluster` must have length equal to the number of
  observations in data. The exception is the GP (only) model
  (`binom = FALSE`), where the `cluster` may either contain a value for
  each observation in the raw data, or for each threshold exceedance in
  the data.

  If `cluster` is not supplied (is `NULL`) then it is assumed that each
  observation forms its own cluster. See **Details** for further
  details.

- use_vcov:

  A logical scalar. Should we use the `vcov` S3 method for `x` (if this
  exists) to estimate the Hessian of the independence loglikelihood to
  be passed as the argument `H` to
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)?
  Otherwise, `H` is estimated inside
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
  using [`optimHess`](https://rdrr.io/r/stats/optim.html).

- ...:

  Further arguments to be passed to the functions in the sandwich
  package
  [`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html)
  (if `cluster = NULL`), or
  [`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  (if `cluster` is not `NULL`).

## Value

An object inheriting from class `"chandwich"`. See
[`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html).
`class(x)` is a vector of length 5. The first 3 components are
`c("lax", "chandwich", "extRemes")`. The remaining 2 components depend
on the model that was fitted. The 4th component is: `"gev"` if
`x$type = "GEV"` or `x$type = "Gumbel"`; `"gp"` if `x$type = "GP"` or
`x$type = "Exponential"`; `"pp"` if `x$type = "PP"`. The 5th component
is `"stat"` if `is.fixedfevd = TRUE` and `"nonstat"` if
`is.fixedfevd = FALSE`.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
for details.

## References

Chandler, R. E. and Bate, S. (2007). Inference for clustered data using
the independence loglikelihood. *Biometrika*, **94**(1), 167-183.
[doi:10.1093/biomet/asm015](https://doi.org/10.1093/biomet/asm015)

Suveges, M. and Davison, A. C. (2010) Model misspecification in peaks
over threshold analysis, *The Annals of Applied Statistics*, **4**(1),
203-221. [doi:10.1214/09-AOAS292](https://doi.org/10.1214/09-AOAS292)

Zeileis (2006) Object-Oriented Computation and Sandwich Estimators.
*Journal of Statistical Software*, **16**, 1-16.
[doi:10.18637/jss.v016.i09](https://doi.org/10.18637/jss.v016.i09)

## See also

[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md):
loglikelihood adjustment for model fits.

## Examples

``` r
# We need the extRemes and distillery packages
got_extRemes <- requireNamespace("extRemes", quietly = TRUE)
got_distillery <- requireNamespace("distillery", quietly = TRUE)

if (got_extRemes & got_distillery) {
  library(extRemes)
  library(distillery)
  # Examples from the extRemes::fevd documentation
  data(PORTw)

  # GEV
  fit0 <- fevd(TMX1, PORTw, units = "deg C", use.phi = TRUE)
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GEV regression
  fitPORTstdmax <- fevd(TMX1, PORTw, scale.fun = ~STDTMAX, use.phi = TRUE)
  adj_fit1 <- alogLik(fitPORTstdmax)
  summary(adj_fit1)
  fitPORTstdmax2 <- fevd(TMX1, PORTw, location.fun = ~STDTMAX,
                         scale.fun = ~STDTMAX, use.phi = TRUE)
  adj_fit2 <- alogLik(fitPORTstdmax2)
  summary(adj_fit2)
  anova(adj_fit0, adj_fit1)
  anova(adj_fit1, adj_fit2)
  anova(adj_fit0, adj_fit2)
  anova(adj_fit0, adj_fit1, adj_fit2)

  # Gumbel
  fit0 <- fevd(TMX1, PORTw, type = "Gumbel", units = "deg C")
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GP
  data(damage)
  fit1 <- fevd(Dam, damage, threshold = 6, type = "GP",
               time.units = "2.05/year")
  adj_fit1 <- alogLik(fit1)
  summary(adj_fit1)

  # Exponential
  fit0 <- fevd(Dam, damage, threshold = 6, type="Exponential",
               time.units = "2.05/year")
  adj_fit0 <- alogLik(fit0)
  summary(adj_fit0)

  # GP non-constant threshold
  data(Fort)
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
              type = "GP")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # Exponential non-constant threshold
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun = ~I(-0.15 * cos(2 * pi * month / 12)),
              type = "Exponential")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # PP model
  fit <- fevd(Prec, Fort, threshold = 0.475, type = "PP", units = "inches")
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # PP non-constant threshold
  fit <- fevd(Prec, Fort, threshold = 0.475,
              threshold.fun=~I(-0.15 * cos(2 * pi * month / 12)),
              type = "PP")
  adj_fit <- alogLik(fit)
  summary(adj_fit)
}
```
