# Dfishery S4 object

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

## Slots inherited from Dfishery

- `nf`:

  Integer, number of fleets

- `Cobs_ymfr`:

  Total fishery catch

- `Csd_ymfr`:

  Lognormal standard deviation of the fishery catch. Only used if
  `Dmodel@condition = "F"`. Default of 0.01.

- `fwt_ymafs`:

  Fishery weight at age. Set to 1 when fleet catch is in units of
  abundance. Set to stock weight at age by default (values at beginning
  of time step).

- `CAAobs_ymafr`:

  Fishery catch at age composition

- `CALobs_ymlfr`:

  Fishery catch at length composition

- `fcomp_like`:

  Character, likelihood for the fishery composition data. See `type`
  argument of
  [`like_comp()`](https://blue-matter.github.io/multiSA/reference/like_comp.md)
  for options

- `CAAN_ymfr`:

  Sample size of the catch at age vector by season if using the
  multinomial or Dirichlet-multinomial likelihoods

- `CALN_ymfr`:

  Sample size of the catch at length vector by season if using the
  multinomial or Dirichlet-multinomial likelihoods

- `CAAtheta_f`:

  Catch at age dispersion parameter if using the Dirichlet-multinomial
  likelihood. Default set to 1.

- `CALtheta_f`:

  Catch at length dispersion parameter if using the
  Dirichlet-multinomial likelihood. Default set to 1.

- `sel_block_yf`:

  Index of dummy fleets to model time blocks of selectivity

- `sel_f`:

  Character vector of the functional form for selectivity. Choose
  between:
  `"logistic_length", "dome_length", "logistic_age", "dome_age", "SB", "B"`

- `Cinit_mfr`:

  Equilibrium seasonal catch prior to the first year. One way to
  initialize the abundance at the start of the first year in the model.
  Default of zero.

- `SC_ymafrs`:

  Stock composition data.

- `SC_aa`:

  Boolean matrix that aggregates age classes for the stock composition
  data. See example.

- `SC_ff`:

  Boolean matrix that aggregates fleets for the stock composition data.
  See example.

- `SC_like`:

  Character, likelihood for the stock composition data. See `type`
  argument of
  [`like_comp()`](https://blue-matter.github.io/multiSA/reference/like_comp.md)
  for options

- `SCN_ymafr`:

  Sample size of the stock composition vector if using the multinomial
  or Dirichlet-multinomial likelihoods

- `SCtheta_f`:

  Stock composition dispersion parameter if using the
  Dirichlet-multinomial likelihood. Default set to 1.

- `SCstdev_ymafrs`:

  Stock composition standard deviation if using the lognormal
  likelihood. Default set to 0.1.

## See also

[MSAdata-class](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
[Dmodel-class](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
[Dstock-class](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
Dfishery-class
[Dsurvey-class](https://blue-matter.github.io/multiSA/reference/Dsurvey-class.md)
[DCKMR-class](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
[Dtag-class](https://blue-matter.github.io/multiSA/reference/Dtag-class.md)

## Examples

``` r
# Aggregate stock composition for ages 1-4 and 5-10 across all fleets
na <- 10
na_SC <- 2
SC_aa <- matrix(0, na_SC, na) # Assumes dim(SC_ymafrs)[3] = na_SC
SC_aa[1, 1:4] <- SC_aa[2, 5:10] <- 1

nf <- 3
nf_SC <- 1
SC_ff <- matrix(1, nf_SC, nf) # Assumes dim(SC_ymafrs)[4] <- nf_SC
```
