# Loglikelihood adjustment for fExtremes fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions `gevFit`,
`gumbelFit` and `gpdFit` in the `fExtremes` package. The model must have
been fitted using maximum likelihood estimation.

## Usage

``` r
# S3 method for class 'fGEVFIT'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'fGPDFIT'
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
`c("lax", "chandwich", "fExtremes")`. The remaining 2 components depend
on the model that was fitted. If `gevFit` or `gumbelFit` was used then
these components are `c("gev", "stat")`. If `gpdFit` was used then these
components are `c("gpd", "stat")`.

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
# We need the fExtremes package
got_fExtremes <- requireNamespace("fExtremes", quietly = TRUE)
if (got_fExtremes) {
  library(fExtremes)

  # GEV
  # An example from the fExtremes::gevFit documentation
  set.seed(4082019)
  x <- fExtremes::gevSim(model = list(xi=0.25, mu=0, beta=1), n = 1000)
  # Fit GEV distribution by maximum likelihood estimation
  fit <- fExtremes::gevFit(x)
  adj_fit <- alogLik(fit)
  summary(adj_fit)

  # GP
  # An example from the fExtremes::gpdFit documentation
  # Simulate GP data
  x <- fExtremes::gpdSim(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000)
  # Fit GP distribution by maximum likelihood estimation
  fit <- fExtremes::gpdFit(x, u = min(x))
  adj_fit <- alogLik(fit)
  summary(adj_fit)
}
```
