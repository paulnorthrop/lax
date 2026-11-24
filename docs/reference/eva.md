# Loglikelihood adjustment for eva fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions `gevrFit` and
`gpdFit` in the eva package.

## Usage

``` r
# S3 method for class 'gevrFit'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'gpdFit'
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
`c("lax", "chandwich", "eva")`. The 4th component depends on which model
was fitted. `"rlarg"` if `gevrFit` was used; `"gpd"` if `gpdFit` was
used. The 5th component is `"stat"` if there are no covariates in the
mode and `"nonstat"` otherwise.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
for details.

In the stationary case (no covariates) the function `gevrFit` and
`gpdFit` in the eva package offer standard errors based on the expected
information or on the observed information, via the argument
`information`. In contrast,
[`alogLik()`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
always bases calculations on the observed information matrix. Therefore,
unadjusted standard errors resulting from
[`alogLik()`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
may be different the corresponding standard errors from `gevrFit` or
`gpdFit`.

For `gevrFit` only GEV fits (`gumbel = FALSE`) are supported.

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
# We need the eva package
got_eva <- requireNamespace("eva", quietly = TRUE)

if (got_eva) {
  library(eva)
  # An example from the eva::gpdFit documentation
  set.seed(7)
  x <- eva::rgpd(2000, loc = 0, scale = 2, shape = 0.2)
  mle_fit <- eva::gpdFit(x, threshold = 4, method = "mle")
  adj_mle_fit <- alogLik(mle_fit)
  summary(adj_mle_fit)

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
  summary(adj_result1)

  # An example from the eva::gevrFit documentation
  set.seed(7)
  x1 <- eva::rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
  result1 <- eva::gevrFit(x1, method = "mle")
  adj_result1 <- alogLik(result1)
  summary(adj_result1)

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
  adj_result2 <- alogLik(result2)
  summary(adj_result2)
}
```
