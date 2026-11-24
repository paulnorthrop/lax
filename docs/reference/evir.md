# Loglikelihood adjustment for evir fits

S3 `alogLik` method to perform loglikelihood adjustment for fitted
extreme value model objects returned from the functions `gev`, `gpd` and
`pot` in the evir package. If `x` was returned from `pot` then the model
will need to be re-fitted using
[`pot_refit`](https://github.com/paulnorthrop/lax/reference/pot_refit.md).

## Usage

``` r
# S3 method for class 'gev'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'gpd'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)

# S3 method for class 'potd'
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
`c("lax", "chandwich", "evir")`. The remaining 2 components depend on
the model that was fitted. If `gev` was used then these components are
`c("gev", "stat")`. If `gpd` was used then these components are
`c("gpd", "stat")`. If
[`pot_refit`](https://github.com/paulnorthrop/lax/reference/pot_refit.md)
was used then these components are `c("potd", "stat")`.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md)
for details.

If `pot` was used then `x` does not contain the raw data that `alogLik`
needs. The model will need to be re-fitted using
[`pot_refit`](https://github.com/paulnorthrop/lax/reference/pot_refit.md)
and the user will be prompted to do this by an error message produced by
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md).

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
# We need the evir package
got_evir <- requireNamespace("evir", quietly = TRUE)
if (got_evir) {
  library(evir)
  # An example from the evir::gev documentation
  data(bmw)
  out <- gev(bmw, "month")
  adj_out <- alogLik(out)
  summary(adj_out)

  # An example from the evir::gpd documentation
  data(danish)
  out <- gpd(danish, 10)
  adj_out <- alogLik(out)
  summary(adj_out)

  # An example from the evir::pot documentation
  # We use lax::pot_refit() to return the input data
  out <- pot_refit(danish, 10)
  adj_out <- alogLik(out)
  summary(adj_out)
}
```
