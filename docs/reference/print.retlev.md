# Print method for retlev object

`print` method for an objects of class `c("retlev", "lax")`.

## Usage

``` r
# S3 method for class 'retlev'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an object of class `c("retlev", "lax")`, a result of a call to
  [`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

- digits:

  The argument `digits` to
  [`print.default`](https://rdrr.io/r/base/print.default.html).

- ...:

  Additional arguments. None are used in this function.

## Value

The argument `x`, invisibly, as for all
[`print`](https://rdrr.io/r/base/print.html) methods.

## Details

Prints the call to
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md)
and the estimates and 100`x$level`% confidence limits for the `x$m`-year
return level.

## Examples

See the examples in
[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).

## See also

[`return_level`](https://github.com/paulnorthrop/lax/reference/return_level.md).
