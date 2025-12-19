# Dsurvey S4 object

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

## See also

[MSAdata-class](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
[Dmodel-class](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
[Dstock-class](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
[Dfishery-class](https://blue-matter.github.io/multiSA/reference/Dfishery-class.md)
Dsurvey-class
[DCKMR-class](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
[Dtag-class](https://blue-matter.github.io/multiSA/reference/Dtag-class.md)
