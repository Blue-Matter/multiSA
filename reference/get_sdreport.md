# Calculate standard errors

A wrapper function to return standard errors. Various numerical
techniques are employed to obtain a positive-definite covariance matrix
in marginal cases.

## Usage

``` r
get_sdreport(
  obj,
  par.fixed,
  exact = FALSE,
  getReportCovariance = FALSE,
  silent = FALSE,
  ...
)
```

## Arguments

- obj:

  The list returned by
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- par.fixed:

  Numeric vector of parameters from which to calculate covariance
  matrix. Optional

- exact:

  Logical, whether to use autodiff or finite-difference approximation
  for the hessian. See details.

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
[`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

A correlation matrix is generated and stored in:
`get_sdreport(obj)$env$corr.fixed`

The hessian is stored in `get_sdreport(obj)$env$hessian`

## Details

Uses [`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) if
`exact = FALSE`. Autodiff with `exact = TRUE` is only available for TMB
models without random effects, but is also memory-intensive.

In numerically marginal cases where the determinant of the Hessian
matrix is less than 0.1, the function will attempt to calculate the
hessian with
[`numDeriv::jacobian()`](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
and the gradient from TMB.

Finally, in other marginal cases where
[`chol()`](https://rdrr.io/r/base/chol.html) identifies a
positive-definite Hessian but
[`solve()`](https://rdrr.io/r/base/solve.html) fails to invert the
matrix, the covariance matrix will be updated with `chol2inv(chol(h))`
