# Simulate data

Simulate data observations from fitted MSA model.

## Usage

``` r
# S4 method for class 'MSAassess'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- nsim:

  Integer, number of simulations

- seed:

  Random number generator seed

- ...:

  Not used

## Value

A list of `nsim` length with data observations
