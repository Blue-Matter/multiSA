# Multi-fleet, multi-area, multi-stock population dynamics model

Project age-structured populations forward in time. Also used by
[`calc_phi_project()`](https://blue-matter.github.io/multiSA/reference/calc_phi_project.md)
to calculate equilibrium abundance and biomass for which there is no
analytic solution due to seasonal movement.

## Usage

``` r
calc_population(
  ny = 10,
  nm = 4,
  na = 20,
  nf = 1,
  nr = 4,
  ns = 2,
  initN_ars = array(1, c(na, nr, ns)),
  mov_ymarrs,
  M_yas = array(0.3, c(ny, na, ns)),
  SRR_s = rep("BH", ns),
  sralpha_s = rep(1e+16, ns),
  srbeta_s = rep(1e+16, ns),
  mat_yas = array(1, c(ny, na, ns)),
  fec_yas = array(1, c(ny, na, ns)),
  Rdev_ys = matrix(1, ny, ns),
  m_spawn = 1,
  m_advanceage = 1,
  delta_s = rep(0, ns),
  natal_rs = matrix(1, nr, ns),
  recdist_rs = matrix(1/nr, nr, ns),
  fwt_ymafs = array(1, c(ny, nm, na, nf, ns)),
  q_fs = matrix(1, nf, ns),
  sel_ymafs = array(1, c(ny, nm, na, nf, ns)),
  condition = c("F", "catch"),
  F_ymfr = array(0, c(ny, nm, nf, nr)),
  Cobs_ymfr = array(1e-08, c(ny, nm, nf, nr)),
  Fmax = 2,
  nitF = 5L
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

- mov_ymarrs:

  Movement array `[y, m, a, r, r, s]`. If missing, uses a diagonal
  matrix (no movement among areas).

- M_yas:

  Natural mortality (per year). Array `[y, a, s]`

- SRR_s:

  Character vector by `s` for the stock recruit relationship. See
  [`calc_recruitment()`](https://blue-matter.github.io/multiSA/reference/calc_recruitment.md)
  for options

- sralpha_s:

  Numeric vector by `s` for the stock recruit alpha parameter

- srbeta_s:

  Numeric vector by `s` for the stock recruit beta parameter

- mat_yas:

  Maturity ogive. Array `[y, a, s]`

- fec_yas:

  Fecundity schedule (spawning output of mature individuals). Array
  `[y, a, s]`

- Rdev_ys:

  Recruitment deviations. Matrix `[y, s]`

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

- fwt_ymafs:

  Fishery weight at age. Array `[y, m, a, f, s]`

- q_fs:

  Relative catchability of stock `s` for fleet `f`. Defaults to 1 if
  missing. Matrix `[f, s]`

- sel_ymafs:

  Fishery selectivity. Array `[y, m, a, f, s]`

- condition:

  Whether the fishing mortality is conditioned on the catch or specified
  F argument.

- F_ymfr:

  Fishing mortality (per season). Array `[y, m, f, r]`. Only used if
  `condition = "F"`.

- Cobs_ymfr:

  Fishery catch (weight). Array `[y, m, f, r]`. Only used if
  `condition = "catch"` to solve for F (see
  [`calc_F()`](https://blue-matter.github.io/multiSA/reference/calc_F.md)).

- Fmax:

  Numeric, the maximum Findex value

- nitF:

  Integer, number of iterations for the Newton-Raphson routine

## Value

A named list containing:

- `N_ymars` Stock abundance

- `F_ymars` Fishing mortality (summed across fleets)

- `F_ymfr` Fishing mortality (by fleet and region)

- `Z_ymars` Total mortality

- `F_ymafrs` Fishing mortality (disaggregated by fleet)

- `CN_ymafrs` Catch at age (abundance)

- `CB_ymfrs` Fishery catch (weight)

- `VB_ymfrs` Vulnerable biomass available to the fishing fleets

- `Nsp_yars` Spawning abundance (in the spawning season)

- `Npsp_yars` Potentail spawners, mature animals in the spawning season
  that do not spawn if outside natal regions

- `S_yrs` Spawning output

- `R_ys` Recruitment

- `penalty` Numeric quadratic penalty if apical fishing mortality (by
  fleet) exceeds `Fmax`. See
  [`calc_F()`](https://blue-matter.github.io/multiSA/reference/calc_F.md).

## Examples

``` r
unfished_pop <- calc_population()
```
