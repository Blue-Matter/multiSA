# If statements compatible with RTMB

Convenience functions that allow taping of gradients in RTMB with if
expressions, following the corresponding `CppAD` functions.

## Usage

``` r
CondExpLt(left, right, if_true, if_false)

CondExpLe(left, right, if_true, if_false)

CondExpGt(left, right, if_true, if_false)

CondExpGe(left, right, if_true, if_false)

CondExpEq(left, right, if_true, if_false)
```

## Arguments

- left:

  Numeric on left hand side of the evaluation

- right:

  Numeric on right hand side of the evaluation

- if_true:

  Numeric if expression is true

- if_false:

  Numeric if expression is false

## Value

Numeric

## Details

Functions should be vectorized.

`CondExpLt` evaluates whether `left < right`

`CondExpLe` evaluates whether `left <= right`

`CondExpGt` evaluates whether `left > right`

`CondExpGe` evaluates whether `left >= right`

`CondExpEq` evaluates whether `left == right`

## Examples

``` r
library(RTMB)
TapeConfig(comparison = "tape")
f <- function(x) CondExpLt(x, 3, 0, x^2)
g <- MakeTape(f, numeric(1))
x <- seq(0, 5)

# Does not work!
f2 <- function(x) if (x < 3) 0 else x^2
g2 <- MakeTape(f2, numeric(1))

# Compare the real answer (deriv) with various values returned by RTMB
data.frame(
  x = x,
  deriv = ifelse(x < 3, 0, 2 * x),
  deriv_f = sapply(x, g$jacobian),
  deriv_f2 = sapply(x, g2$jacobian)
)
#>   x deriv deriv_f deriv_f2
#> 1 0     0       0        0
#> 2 1     0       0        0
#> 3 2     0       0        0
#> 4 3     6       6        0
#> 5 4     8       8        0
#> 6 5    10      10        0
```
