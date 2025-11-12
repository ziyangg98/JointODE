# Print Summary of MarginalODE Fit

Print Summary of MarginalODE Fit

## Usage

``` r
# S3 method for class 'summary.MarginalODE'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  An object of class `summary.MarginalODE`

- digits:

  Number of digits to display (default: max(3L, getOption("digits") -
  3L))

- signif.stars:

  Logical; show significance stars (default:
  getOption("show.signif.stars"))

- ...:

  Additional arguments passed to `printCoefmat`
