# Inference for the Bernoulli distribution

Functions involved in making inferences about the probability of success
in a Bernoulli distribution.

## Usage

``` r
fit_bernoulli(data)

# S3 method for class 'bernoulli'
logLikVec(object, pars = NULL, ...)

# S3 method for class 'bernoulli'
nobs(object, ...)

# S3 method for class 'bernoulli'
coef(object, ...)

# S3 method for class 'bernoulli'
vcov(object, ...)

# S3 method for class 'bernoulli'
logLik(object, ...)

# S3 method for class 'bernoulli'
alogLik(x, cluster = NULL, use_vcov = TRUE, ...)
```

## Arguments

- data:

  A numeric vector of outcomes from Bernoulli trials: 0 for a failure, 1
  for a success. Alternatively, a logical vector with FALSE for a
  failure and TRUE for a success.

- pars:

  A numeric parameter vector of length 1 containing the value of the
  Bernoulli success probability.

- ...:

  Further arguments to be passed to the functions in the sandwich
  package
  [`meat`](https://zeileis.codeberg.page/sandwich/reference/meat.html)
  (if `cluster = NULL`), or
  [`meatCL`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  (if `cluster` is not `NULL`).

- x, object:

  A fitted model object returned from `fit_bernoulli()`.

- cluster:

  A vector or factor indicating from which cluster each observation in
  `data` originates.

- use_vcov:

  A logical scalar. Should we use the `vcov` S3 method for `x` (if this
  exists) to estimate the Hessian of the independence loglikelihood to
  be passed as the argument `H` to
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)?
  Otherwise, `H` is estimated inside
  [`adjust_loglik`](https://github.com/paulnorthrop/chandwich/reference/adjust_loglik.html)
  using [`optimHess`](https://rdrr.io/r/stats/optim.html).

## Value

`fit_bernoulli` returns an object of class `"bernoulli"`, a list with
components: `logLik, mle, nobs, vcov, data, obs_data`, where `data` are
the input data and `obs_data` are the input data after any missing
values have been removed, using
[`na.omit`](https://rdrr.io/r/stats/na.fail.html).

`logLikVec.bernoulli` returns an object of class `"logLikVec"`, a vector
of length `length(data)` containing the loglikelihood contributions from
the individual observations in `data`.

## Details

`fit_bernoulli`: fit a Bernoulli distribution

`logLikVec.bernoulli`: calculates contributions to a loglikelihood based
on the Bernoulli distribution. The loglikelihood is calculated up to an
additive constant.

`nobs, coef, vcov` and `logLik` methods are provided.

## See also

[`Binomial`](https://rdrr.io/r/stats/Binomial.html). The Bernoulli
distribution is the special case where `size = 1`.

## Examples

``` r
# Set up data
x <- exdex::newlyn
u <- quantile(x, probs = 0.9)
exc <- x > u

# Fit a Bernoulli distribution
fit <- fit_bernoulli(exc)

# Calculate the loglikelihood at the MLE
res <- logLikVec(fit)

# The logLik method sums the individual loglikelihood contributions.
logLik(res)
#> 'log Lik.' -939.9109 (df=1)

# nobs, coef, vcov, logLik methods for objects returned from fit_bernoulli()
nobs(fit)
#> [1] 2894
coef(fit)
#>       prob 
#> 0.09986178 
vcov(fit)
#>              prob
#> prob 3.106061e-05
logLik(fit)
#> 'log Lik.' -939.9109 (df=1)

# Adjusted loglikelihood
# Create 5 clusters each corresponding approximately to 1 year of data
cluster <- rep(1:5, each = 579)[-1]
afit <- alogLik(fit, cluster = cluster, cadjust = FALSE)
summary(afit)
#>          MLE       SE adj. SE
#> prob 0.09986 0.005573 0.01831
```
