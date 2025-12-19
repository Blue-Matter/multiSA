# Additional methods for AD types

Methods for RTMB AD class

## Usage

``` r
# S4 method for class 'ad,matrix'
x %*% y
```

## Arguments

- x:

  AD object

- y:

  Non-AD matrix

## Functions

- `x %*% y`: Matrix product function implemented for mixed AD and non-AD
  objects with `colSums(x * y)`. See
  [ADmatrix](https://rdrr.io/pkg/RTMB/man/ADmatrix.html).
