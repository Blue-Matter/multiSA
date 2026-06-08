# Optimize RTMB model

A convenient function to fit a RTMB model with
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

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

  Deprecated.

- do_sd:

  Deprecated.

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

A named list, output of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

## See also

[`get_sdreport()`](https://blue-matter.github.io/multiSA/reference/get_sdreport.md)
