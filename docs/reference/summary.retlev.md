# Summary method for a `"retlev"` object

`summary` method for an objects of class `c("retlev", "lax")`.

## Usage

``` r
# S3 method for class 'retlev'
summary(object, digits, ...)
```

## Arguments

- object:

  an object of class `c("retlev", "lax")`, a result of a call to
  [`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

- digits:

  An integer. Used for number formatting with `signif`. If `digits` is
  not specified (i.e. [`missing`](https://rdrr.io/r/base/missing.html))
  then [`signif()`](https://rdrr.io/r/base/Round.html) will not be
  called (i.e. no rounding will be performed).

- ...:

  Additional arguments. None are used in this function.

## Value

Returns a list containing the list element `object$call` and a numeric
matrix `matrix` containing the MLE and estimated SE of the return level.

## Examples

See the examples in
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

## See also

[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).
