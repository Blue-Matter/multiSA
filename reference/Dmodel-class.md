# Dmodel S4 object

S4 class that organizes the various data inputs for the MSA model.
`MSAdata` simply inherits the slots from 6 component classes: `Dmodel`,
`Dstock`, `Dfishery`, `Dsurvey` `DCKMR`, and `Dtag`, where the `D`-
prefix denotes an object for data inputs (or model configuration).

## Details

For convenience, most arrays and matrices have the associated dimensions
in the variable name. For example, `Cobs_ymfr` represents observed catch
with the dimension following the underscore, following this template:

|     |         |
|-----|---------|
| `y` | Year    |
| `m` | Season  |
| `a` | Age     |
| `r` | Region  |
| `f` | Fishery |
| `i` | Index   |
| `s` | Stock   |

## Slots inherited from Dmodel

- `ny`:

  Integer, number of years

- `nm`:

  Integer, number of seasons

- `na`:

  Integer, number of ages. The first age class is zero and the last age
  class (plus group is age `na - 1`).

- `nl`:

  Integer, number of length bins. Set to zero if lengths are not
  modeled.

- `nr`:

  Integer, number of spatial regions

- `ns`:

  Integer, number of stocks

- `lbin`:

  Vector of lower boundary of length bins. Length `nl + 1`

- `lmid`:

  Vector of midpoint of length bins. Length `nl`

- `Fmax`:

  Numeric, maximum allowable instantaneous fishing mortality rate (units
  of per season). Defaults to 3.

- `y_phi`:

  Integer, the year from which to obtain values of natural mortality and
  fecundity for the unfished stock-recruit replacement line (`phi`).
  Relevant if natural mortality or fecundity are time-varying. Defaults
  to 1.

- `scale_s`:

  Vector, length `ns`. Multiplicative scaling factor that informs
  relative stock size to aid parameter estimation. Larger values implies
  larger stocks. Default set to 1. See
  [`make_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md).

- `nyinit`:

  Integer, number of years of spool-up to calculate equilibrium unfished
  and starting conditions for the population model to account for
  seasonal and spatial dynamics. The numerical spool-up is not needed
  when both `nm = 1` and `nr = 1`, i.e., `nyinit = 1`. Otherwise, set to
  `2 * na` by default.

- `condition`:

  Character, either to specify the model estimates fishing mortality as
  a parameter (`"F"`, default) or equal to the catch (`"catch"`).

- `nitF`:

  Integer, number of iterations to solve Baranov catch equation from
  observed catch if `condition = "catch"`. Defaults to 5.

- `y_Fmult_f`:

  Integer vector by fleet, the year in which to directly estimate F.
  Choose a year/season/region combination when the catch is average
  relative to the time series. Only used if `condition = "F"`.

- `m_Fmult_f`:

  Integer vector by fleet, the season in which to directly estimate F.
  Choose a year/season/region combination when the catch is average
  relative to the time series. Only used if `condition = "F"`.

- `r_Fmult_f`:

  Integer vector by fleet, the region in which to directly estimate F.
  Choose a year/season/region combination when the catch is average
  relative to the time series. Only used if `condition = "F"`.

- `pbc_rdev_ys`:

  Numeric matrix, for the fraction of lognormal bias correction
  (`-0.5 * sd_r^2`) applied to the recruitment estimates in the model.
  Typically between 0-1, with default of 1.

- `pbc_initrdev_as`:

  Numeric matrix, for the fraction of lognormal bias correction
  (`-0.5 * sd_r^2`) applied to the initial abundance vector in the
  model. Typically between 0-1, with default of 1.

- `prior`:

  Character vector to be evaluated in the model to return the log prior
  for a parameter. See example in documentation for
  [prior](https://blue-matter.github.io/multiSA/reference/prior.md).

- `nyret`:

  Integer, number of recent years of data to remove from the likelihood
  for retrospective analysis (positive numbers). Default is zero.

## See also

[MSAdata-class](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
Dmodel-class
[Dstock-class](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
[Dfishery-class](https://blue-matter.github.io/multiSA/reference/Dfishery-class.md)
[Dsurvey-class](https://blue-matter.github.io/multiSA/reference/Dsurvey-class.md)
[DCKMR-class](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
[Dtag-class](https://blue-matter.github.io/multiSA/reference/Dtag-class.md)
