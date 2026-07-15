# MSAdata S4 object

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

## Slots

- `Dmodel`:

  Class
  [Dmodel](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
  containing parameters for model structure (number of years, ages,
  etc.)

- `Dstock`:

  Class
  [Dstock](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
  containing stock parameters (growth, natural mortality, etc.)

- `Dfishery`:

  Class
  [Dfishery](https://blue-matter.github.io/multiSA/reference/Dfishery-class.md)
  containing fishery data (catch, size and stock composition, etc.)

- `Dsurvey`:

  Class
  [Dsurvey](https://blue-matter.github.io/multiSA/reference/Dsurvey-class.md)
  containing survey data (indices of abundance)

- `DCKMR`:

  Class
  [DCKMR](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
  containing genetic close-kin data

- `Dtag`:

  Class
  [Dtag](https://blue-matter.github.io/multiSA/reference/Dtag-class.md)
  containing tagging data

- `Dlabel`:

  Class
  [Dlabel](https://blue-matter.github.io/multiSA/reference/Dlabel-class.md)
  containing names for various dimensions. Used for plotting.

- `Misc`:

  List for miscellaneous inputs as needed

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

  Character vector of the functional form for selectivity for each
  fleet. See `type` argument in
  [`conv_selpar()`](https://blue-matter.github.io/multiSA/reference/conv_selpar.md)
  for all options, for example,
  `"logistic_length", "dome_length", "logistic_age", "dome_age", "total", "mature"`

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

- `lambdaCobs_f`:

  Likelihood weight (by fleet) for the catch. Default is 1.

- `lambdaCAA_f`:

  Likelihood weight (by fleet) for the catch at age. Default is 1.

- `lambdaCAL_f`:

  Likelihood weight (by fleet) for the catch at length. Default is 1.

- `lambdaSC_f`:

  Likelihood weight (by fleet) for the stock composition. Default is 1.

## Slots inherited from Dsurvey

- `ni`:

  Integer, number of indices of abundance. Zero is possible.

- `Iobs_ymi`:

  Observed indices

- `Isd_ymi`:

  Lognormal standard deviation of the observed indices

- `unit_i`:

  Character vector, units of the index. Set to `"B"` to use stock weight
  at age (default) or `"N"` for abundance (numbers).

- `IAAobs_ymai`:

  Array, survey age composition

- `IALobs_ymli`:

  Array, survey length composition

- `icomp_like`:

  Character, likelihood for the composition data. See
  [`like_comp()`](https://blue-matter.github.io/multiSA/reference/like_comp.md)
  for options

- `IAAN_ymi`:

  Array, sample size of the index age composition by season if using the
  multinomial or Dirichlet-multinomial likelihoods

- `IALN_ymi`:

  Array, sample size of the index length composition by season if using
  the multinomial or Dirichlet-multinomial likelihoods

- `IAAtheta_i`:

  Numeric vector, survey age composition dispersion parameter if using
  the Dirichlet-multinomial likelihood

- `IALtheta_i`:

  Numeric vector, index length composition dispersion parameter if using
  the Dirichlet-multinomial likelihood

- `samp_irs`:

  Boolean array that specifies the regions and stocks sampled by the
  index. `samp[i, r, s]` indicates whether index `i` operates in region
  `r` and catches stock `s`.

- `sel_i`:

  Character vector, functional forms for selectivity. See `"type"`
  argument in
  [`conv_selpar()`](https://blue-matter.github.io/multiSA/reference/conv_selpar.md)
  for options.

- `delta_i`:

  Numeric vector, The elapsed fraction of time in the seasonal time step
  (between 0 - 1) when the index samples the population. Set to a
  negative number (-1) to sample over the duration of the timestep,
  i.e., `(1 - exp(-Z))/Z`.

- `lambdaI_i`:

  Likelihood weight (by survey) for the indices. Default is 1.

- `lambdaIAA_i`:

  Likelihood weight (by fleet) for the index age composition. Default is
  1.

- `lambdaIAL_i`:

  Likelihood weight (by fleet) for the index length composition. Default
  is 1.

## Slots inherited from DCKMR

- `POP_s`:

  A list by stock of data frames for parent-offspring pairs. Each row in
  the data frame corresponds to a "sampling unit" defined by the
  columns:

  |     |                                |
  |-----|--------------------------------|
  | `a` | Capture year of parent         |
  | `t` | Age at capture of parent       |
  | `y` | Birth year of offspring        |
  | `n` | Number of pairwise comparisons |
  | `m` | Number of POPs                 |

- `HSP_s`:

  A list by stock of data frames for half-sibling pairs. Each row in the
  data frame corresponds to a "sampling unit" defined by the columns:

  |      |                                |
  |------|--------------------------------|
  | `yi` | Birth year of older sibling    |
  | `yj` | Birth year of younger sibling  |
  | `n`  | Number of pairwise comparisons |
  | `m`  | Number of HSPs                 |

- `CKMR_like`:

  Character, likelihood for the POP and HSP sampling units. See `type`
  argument in
  [`like_CKMR()`](https://blue-matter.github.io/multiSA/reference/like_CKMR.md)
  for options.

## Slots inherited from Dtag

- `tag_ymarrs`:

  Array. Number of tags that move between regions. Informs movement
  matrices of stocks between time steps.

- `tag_ymars`:

  Array. Number of tags distributed among regions. Informs stock
  distribution (within time step).

- `tag_yy`:

  Boolean matrix that aggregates years for the tag data. Only used for
  the tag movement array `tag_ymarrs`.

- `tag_aa`:

  Boolean matrix that aggregates ages for the tag data.

- `tag_like`:

  Character. Likelihood for the tagging data, either the vector of
  proportions by region of origin for `tag_ymarrs`, or by region of
  stock distribution for `tag_ymars`. See `type` argument of
  [`like_comp()`](https://blue-matter.github.io/multiSA/reference/like_comp.md)
  for options

- `tagN_ymars`:

  Array. Sample size of the tag movement vectors if using the
  multinomial or Dirichlet-multinomial likelihoods.

- `tagN_ymas`:

  Array. Sample size of the tag distribution vectors if using the
  multinomial or Dirichlet-multinomial likelihoods.

- `tagtheta_s`:

  Array. Tag dispersion parameter (by stock) if using the
  Dirichlet-multinomial likelihoods. Default set to 1.

- `tagstdev_s`:

  Array. Tag standard deviation (by stock) if using the lognormal
  likelihood. Default set to 0.1.

- `lambdaTag_s`:

  Array. Likelihood weight for the tag data (by stock). Default is 1.

## Slots inherited from Dlabel

- `year`:

  Vector of years. Length `Dmodel@ny`

- `season`:

  Vector of season names. Length `Dmodel@nm`

- `age`:

  Vector of ages. Length `Dmodel@na`

- `region`:

  Vector of region names. Length `Dmodel@nr`

- `stock`:

  Vector of stock names. Length `Dmodel@ns`

- `fleet`:

  Vector of fleet names. Length `Dfishery@nf`

- `index`:

  Vector of index of abundance names. Length `Dsurvey@ni`

## See also

MSAdata-class
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
[Dmodel-class](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
[Dstock-class](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
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

# Aggregate stock composition for ages 1-4 and 5-10 across all fleets
na <- 10
na_SC <- 2
SC_aa <- matrix(0, na_SC, na) # Assumes dim(SC_ymafrs)[3] = na_SC
SC_aa[1, 1:4] <- SC_aa[2, 5:10] <- 1

nf <- 3
nf_SC <- 1
SC_ff <- matrix(1, nf_SC, nf) # Assumes dim(SC_ymafrs)[4] <- nf_SC
```
