# Jitter starting values from fitted model

Run additional model fits with jittered starting values.

## Usage

``` r
do_jitter(x, n = 1, use_fitted = TRUE, amount = NULL, cores = 1, seed, ...)
```

## Arguments

- x:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- n:

  Integer, number of jittered model runs

- use_fitted:

  Logical, whether to jitter from estimated parameters (`TRUE`) or
  original starting values (`FALSE`)

- amount:

  Numeric or NULL, passed to
  [`base::jitter()`](https://rdrr.io/r/base/jitter.html)

- cores:

  Integer, number of CPUs for parallel processing

- seed:

  Integer, for replicating the sampling function. Optional.

- ...:

  Other arguments to pass to
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

## Value

List (length `n`) containing
[MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
objects

## Details

The new starting parameters are: `pars + amount * runif(n, -1, 1)` where
`pars` are either the fitted values or original starting values
depending on `use_fitted`.
