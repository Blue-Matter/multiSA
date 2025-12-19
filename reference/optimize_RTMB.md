# Optimize RTMB model

A convenient function that fits a RTMB model and calculates standard
errors.

## Usage

``` r
optimize_RTMB(
  obj,
  hessian = FALSE,
  restart = 0,
  do_sd = TRUE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  lower = -Inf,
  upper = Inf,
  silent = FALSE
)
```

## Arguments

- obj:

  The list returned by
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- hessian:

  Logical, whether to pass the Hessian function `obj$he` to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). Only used if
  there are no random effects in the model.

- restart:

  Integer, the maximum number of additional attempts to fit the model.
  See details.

- do_sd:

  Logical, whether to calculate standard errors through
  [`get_sdreport()`](https://blue-matter.github.io/multiSA/reference/get_sdreport.md)

- control:

  List of options passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

- lower:

  Lower bounds of parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

- upper:

  Upper bounds of parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

- silent:

  Logical, whether to report progress to console

## Value

A named list: "opt" is the output of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) and "SD" is the
output of
[`get_sdreport()`](https://blue-matter.github.io/multiSA/reference/get_sdreport.md)

## Details

Argument `restart` allows for recursive model fitting to obtain
convergence, through the following procedure:

1.  Optimize model with
    [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

2.  Determine convergence, defined by
    [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
    by whether the Cholesky decomposition of the covariance matrix is
    possible.

3.  If convergence is not achieved, jitter parameter estimates with
    multiplicative factor `rlnorm(mean = 0, sd = 1e-3)` and return to
    step 1.

## See also

[`get_sdreport()`](https://blue-matter.github.io/multiSA/reference/get_sdreport.md)
