# Calculate standard errors

A wrapper function to return standard errors. Various numerical
techniques are employed to obtain a positive-definite covariance matrix.

## Usage

``` r
get_sdreport(obj, getReportCovariance = FALSE, silent = FALSE, ...)
```

## Arguments

- obj:

  The list returned by
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- getReportCovariance:

  Logical, passed to
  [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- silent:

  Logical, whether to report progress to console. See details.

- ...:

  Other arguments to
  [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
  besides `par.fixed, hessian.fixed, getReportCovariance`

## Value

Object returned by
[`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html). A
correlation matrix is generated and stored in:
`get_sdreport(obj)$env$corr.fixed`

## Details

In marginal cases where the determinant of the Hessian matrix is less
than 0.1, several steps are utilized to obtain a positive-definite
covariance matrix, in the following order:

1.  Invert hessian returned by `h <- obj@he(obj$env.last.par.best)`
    (skipped in models with random effects)

2.  Invert hessian returned by
    `h <- stats::optimHess(obj$env.last.par.best, obj$fn, obj$gr)`

3.  Invert hessian returned by
    `h <- numDeriv::jacobian(obj$gr, obj$env.last.par.best)`

4.  Calculate covariance matrix from `chol2inv(chol(h))`
