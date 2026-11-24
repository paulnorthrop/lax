# Evaluate loglikelihood contributions from specific observations

Generic function for calculating loglikelihood contributions from
individual observations for a fitted model.

## Usage

``` r
logLikVec(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Further arguments.

## Value

An object of class `"logLikVec"`, a vector containing individual
loglikelihood contributions.

## Details

See
[`alogLik`](https://github.com/paulnorthrop/lax/reference/alogLik.md):
loglikelihood adjustment for model fits.

## Examples

See the example in
[`bernoulli`](https://github.com/paulnorthrop/lax/reference/bernoulli.md).
