# Sum loglikelihood contributions from individual observations

S3 logLik method for logLikVec objects.

## Usage

``` r
# S3 method for class 'logLikVec'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"logLikVec"` return from a `logLikVec` method.

- ...:

  Further arguments.

## Value

An object of class `"logLik"`. This is the value of the loglikelihood,
with attributes `"df"` (degrees of freedom) giving the number of free
parameters, and `"nobs"` giving the number of observations.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md):
loglikelihood adjustment for model fits.

## Examples

See the example in
[`bernoulli`](https://github.com/paulnorthrop/lax/reference/bernoulli.md).
