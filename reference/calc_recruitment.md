# Calculate recruitment from stock-recruit function

Calculate recruitment from stock-recruit function

## Usage

``` r
calc_recruitment(x, SRR = c("BH", "Ricker"), eq = FALSE, ...)
```

## Arguments

- x:

  Numeric, either the spawning output or the equilibrium spawners per
  recruit, from which the recruitment will be calculated. See argument
  `eq`.

- SRR:

  Character to indicate the functional form of the stock recruit
  function

- eq:

  Logical, indicates whether `x` is the spawning output (`FALSE`) or
  equilibrium spawners per recruit (`TRUE`)

- ...:

  Parameters of the SRR function. Provide one of two sets of variables:

  1.  `h`, `R0` and `phi0`, or

  2.  `a` and `b` (alpha, beta values)

## Value

Numeric of length `x`

## Examples

``` r
calc_recruitment(10, SRR = "Ricker", a = 2, b = 0.5)
#> [1] 0.1347589
calc_recruitment(10, SRR = "Ricker", h = 0.9, R0 = 1, phi0 = 1)
#> [1] 4.480838e-07
```
