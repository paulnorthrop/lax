# Plot diagnostics for a retlev object

`plot` method for an objects of class `c("retlev", "lax")`.

## Usage

``` r
# S3 method for class 'retlev'
plot(x, y = NULL, level = NULL, legend = TRUE, digits = 3, plot = TRUE, ...)
```

## Arguments

- x:

  an object of class `c("retlev", "lax")`, a result of a call to
  [`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md),
  using `prof = TRUE`.

- y:

  Not used.

- level:

  A numeric scalar in (0, 1). The confidence level required for the
  confidence interval for the `m`-year return level. If `level` is not
  supplied then `x$level` is used. `level` must be no larger than
  `x$level`.

- legend:

  A logical scalar. Should we add a legend (in the top right of the
  plot) that gives the approximate values of the MLE and 100`level`%
  confidence limits?

- digits:

  An integer. Passed to `signif` to round the values in the legend.

- plot:

  A logical scalar. If `TRUE` then the plot is produced. Otherwise, it
  is not, but the MLE and confidence limits are returned.

- ...:

  Further arguments to be passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A numeric vector of length 3 containing the lower 100`level`% confidence
limit, the MLE and the upper 100`level`% confidence limit.

## Details

Plots the profile loglikelihood for a return level, provided that `x`
returned by a call to
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md)
using `prof = TRUE`. Horizontal lines indicate the values of the
maximised loglikelihood and the critical level used to calculate the
confidence limits. If `level` is smaller than `x$level` then approximate
100`level`% confidence limits are recalculated based on the information
contained in `x$for_plot`.

## Examples

See the examples in
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

## See also

[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md)
to perform inferences about return levels.
