# MSAassess S4 object

S4 object that returns output from MSA model

## Slots

- `obj`:

  RTMB object returned by
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- `opt`:

  List returned by
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

- `SD`:

  List returned by
  [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)

- `report`:

  List of model output at the parameter estimates, returned by
  `obj$report(obj$env$last.par.best)`

- `Misc`:

  List, miscellaneous items
