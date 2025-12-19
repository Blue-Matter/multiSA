# Calculate model residuals

Extract residuals from fitted model

## Usage

``` r
# S3 method for class 'MSAassess'
residuals(object, vars, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- vars:

  Character vector to indicate which residuals will be calculated.
  Available choices from
  [MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
  object are: "Cinit_mfr", "Cobs_ymfr", "CAAobs_ymafr", "CALobs_ymlfr",
  "Iobs_ymi", "IAAobs_ymai", "IALobs_ymli", "SC_ymafrs"

- type:

  Character, 'response' for the `log(observed/predicted)` values or
  'pearson' for calculating Z-scores. Composition data always use
  'pearson'.

- ...:

  Not used

## Value

A named list based on `vars` argument.
