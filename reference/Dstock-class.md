# Dstock S4 object

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

## Slots inherited from Dstock

- `m_spawn`:

  Integer, season of spawning. Defaults to 1. The progeny will enter the
  age structure in the same season.

- `m_advanceage`:

  Integer, season in which to advance age classes (corresponding to
  calendar year) to the next age class. Defaults to 1.

- `len_ymas`:

  Length-at-age. Only needed if `Dmodel@nl > 0` to fit length
  composition (you may want to calculate the growth at the middle of the
  time step).
  [`calc_growth()`](https://blue-matter.github.io/multiSA/reference/calc_growth.md)
  may be a helpful function.

- `sdlen_ymas`:

  Standard deviation in length-at-age

- `LAK_ymals`:

  Length-at-age probability array. If empty, values will be calculated
  by
  [`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
  with
  [`calc_LAK()`](https://blue-matter.github.io/multiSA/reference/calc_LAK.md).

- `mat_yas`:

  Proportion mature by age class. Ignored if maturity ogive is
  estimated, e.g., when fitting to close-kin genetic data.

- `swt_ymas`:

  Stock weight-at-age (at beginning of time step). See
  [`calc_growth()`](https://blue-matter.github.io/multiSA/reference/calc_growth.md)
  example.

- `fec_yas`:

  Fecundity, i.e., spawning output, of mature animals. Default uses
  stock weight at age.

- `M_yas`:

  Natural mortality. Instantaneous units per year. Ignored if M is
  estimated.

- `SRR_s`:

  Character vector of stock-recruit relationship by stock. See `SRR`
  argument in
  [`calc_recruitment()`](https://blue-matter.github.io/multiSA/reference/calc_recruitment.md)
  for options.

- `delta_s`:

  Fraction of season that elapses when spawning occurs, e.g., midseason
  spawning occurs when `delta_s = 0.5`. Default is zero.

- `presence_rs`:

  Logical matrix indicating presence/absence of stock `s` in region `r`.
  Used to constrain movement matrix. Default is TRUE for all stocks and
  regions.

- `natal_rs`:

  The fraction of the mature stock `s` in region `r` that spawns at time
  of spawning. See example. Default is 1 for all stocks and regions.

## See also

[MSAdata-class](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
[Dmodel-class](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
Dstock-class
[Dfishery-class](https://blue-matter.github.io/multiSA/reference/Dfishery-class.md)
[Dsurvey-class](https://blue-matter.github.io/multiSA/reference/Dsurvey-class.md)
[DCKMR-class](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
[Dtag-class](https://blue-matter.github.io/multiSA/reference/Dtag-class.md)

## Examples

``` r
# Set natal_rs matrix so that the spawning output of stock 1 is
# calculated from mature animals present in regions 1, 2.
# Similarly for stock 2, spawning output from areas 2 and 3.
nr <- 4
ns <- 2
natal_rs <- matrix(0, nr, ns)
natal_rs[1:2, 1] <- natal_rs[2:3, 2] <- 1
```
