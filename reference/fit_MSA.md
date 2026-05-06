# Fit MSA model

Wrapper function that calls RTMB to create the model and perform the
numerical optimization

## Usage

``` r
fit_MSA(
  x,
  parameters,
  map = list(),
  random = NULL,
  run_model = TRUE,
  do_sd = TRUE,
  report = TRUE,
  silent = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  ...
)
```

## Arguments

- x:

  Data object. Class
  [MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md),
  validated by
  [`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md).
  Alternatively,
  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  that will be fitted again.

- parameters:

  List of parameters, e.g., returned by
  [`make_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md)
  and validated by
  [`check_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md).

- map:

  List of parameters indicated whether they are fixed and how they are
  shared, e.g., returned by
  [`make_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md).
  See
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

- random:

  Character vector indicating the parameters that are random effects,
  e.g., returned by
  [`make_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md).

- run_model:

  Logical, whether to fit the model through
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- do_sd:

  Logical, whether to calculate the standard errors with
  [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

- report:

  Logical, whether to return the report list with
  `obj$report(obj$env$last.par.best)`.

- silent:

  Logical, whether to report progress to console. **Not passed to
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).**
  Recommend to set to `TRUE` to speed up run time, e.g., when running
  simulations, multiple fits, etc.

- control:

  Passed to [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

- ...:

  Other arguments to
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

## Value

A
[MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
object.

## See also

[`report()`](https://blue-matter.github.io/multiSA/reference/report.md)
[`retrospective()`](https://blue-matter.github.io/multiSA/reference/retrospective.md)
