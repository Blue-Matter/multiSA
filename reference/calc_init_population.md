# Initial population projection

Project a population forward in time using
[`calc_population()`](https://blue-matter.github.io/multiSA/reference/calc_population.md),
an alternative to
[`calc_phi_project()`](https://blue-matter.github.io/multiSA/reference/calc_phi_project.md)
to establish initial age structure.

## Usage

``` r
calc_init_population(
  ny = 10,
  nm = 4,
  na = 20,
  nf = 1,
  nr = 4,
  ns = 1,
  initN_ars = array(1, c(na, nr, ns)),
  F_mfr = array(0, c(nm, nf, nr)),
  sel_mafs = array(1, c(nm, na, nf, ns)),
  fwt_mafs = array(1, c(nm, na, nf, ns)),
  q_fs = matrix(1, nf, ns),
  M_as,
  mov_marrs,
  mat_as,
  fec_as,
  SRR_s,
  sralpha_s,
  srbeta_s,
  m_spawn = 1,
  m_advanceage = 1,
  delta_s = rep(0, ns),
  natal_rs = matrix(1, nr, ns),
  recdist_rs = matrix(1/nr, nr, ns)
)
```

## Arguments

- ny:

  Integer, number of years for the projection

- nm:

  Integer, number of seasons

- na:

  Integer, number of age classes

- nf:

  Integer, number of fleets

- nr:

  Integer, number of regions

- ns:

  Integer, number of stocks

- initN_ars:

  Abundance in the first year, first season. Array `[a, r, s]`

- F_mfr:

  Equilibrium fishing mortality (per season). Matrix `[m, f, r]`

- sel_mafs:

  Selectivity by season, age, fleet, stock. Array `[m, a, f, s]`

- fwt_mafs:

  Fishery weight array by season, age, fleet, stock. Array
  `[m, a, r, r]`. Can be used calculate yield per recruit.

- q_fs:

  Relative catchability of stock `s` for fleet `f`. Defaults to 1 if
  missing. Matrix `[f, s]`

- M_as:

  Natural mortality. Matrix `[a, s]`

- mov_marrs:

  Movement array `[m, a, r, r, s]`. If missing, uses a diagonal matrix
  (no movement among areas).

- mat_as:

  Maturity at age. Matrix `[a, s]`

- fec_as:

  Fecundity at age. Matrix `[a, s]`

- SRR_s:

  Character vector by `s` for the stock recruit relationship. See
  [`calc_recruitment()`](https://blue-matter.github.io/multiSA/reference/calc_recruitment.md)
  for options

- sralpha_s:

  Numeric vector by `s` for the stock recruit alpha parameter

- srbeta_s:

  Numeric vector by `s` for the stock recruit beta parameter

- m_spawn:

  Integer, season of spawning

- m_advanceage:

  Integer, season at which to advance integer year age classes

- delta_s:

  Numeric vector by `s`. Fraction of season that elapses when spawning
  occurs, e.g., midseason spawning when `delta_s = 0.5`.

- natal_rs:

  Matrix `[r, s]`. The fraction of the mature stock `s` in region `r`
  that spawns at time of spawning. See example in
  [Dstock](https://blue-matter.github.io/multiSA/reference/Dstock-class.md).

- recdist_rs:

  Matrix `[r, s]`. The fraction of the incoming recruitment of stock `s`
  that settles in region `r`.

## Value

A named list returned by
[`calc_population()`](https://blue-matter.github.io/multiSA/reference/calc_population.md).

## Details

The initial population vector will be the survival at age evenly divided
by the number of regions `nr`.
