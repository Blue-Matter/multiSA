# Softmax function

Takes a vector of real numbers and returns the corresponding vector of
probabilities

## Usage

``` r
softmax(eta, log = FALSE)
```

## Arguments

- eta:

  Vector

- log:

  Logical, whether to return the value of the logarithm

## Value

Numeric, vector length of `eta`: \\\exp(\eta)/\sum\exp(\eta)\\

## Details

Uses `multiSA:::logspace.add` for numerical stability
