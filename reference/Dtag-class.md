# Dtag S4 object

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

## See also

[MSAdata-class](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
[`check_data()`](https://blue-matter.github.io/multiSA/reference/check_data.md)
[Dmodel-class](https://blue-matter.github.io/multiSA/reference/Dmodel-class.md)
[Dstock-class](https://blue-matter.github.io/multiSA/reference/Dstock-class.md)
[Dfishery-class](https://blue-matter.github.io/multiSA/reference/Dfishery-class.md)
[Dsurvey-class](https://blue-matter.github.io/multiSA/reference/Dsurvey-class.md)
[DCKMR-class](https://blue-matter.github.io/multiSA/reference/DCKMR-class.md)
Dtag-class
