# Selectivity at age and length

Calculate selectivity at age and length from a matrix of parameters.

- `conv_selpar()` converts parameters from log or logit space to real
  units.

- `calc_sel_len()` calculates selectivity at length.

- `calc_fsel_age()` calculates selectivity at age for fisheries, and
  coordinates dummy fleets.

- `calc_isel_age()` calculates selectivity at age for indices, and can
  map selectivity from fisheries or population parameters (e.g, mature
  or total biomass).

## Usage

``` r
conv_selpar(x, type, maxage, maxL)

calc_sel_len(sel_par, lmid, type, fsel_type, fsel_len)

calc_fsel_age(
  sel_len,
  LAK,
  type,
  sel_par,
  sel_block = seq(1, length(type)),
  mat,
  a = seq(1, nrow(LAK)) - 1
)

calc_isel_age(
  sel_len,
  LAK,
  type,
  sel_par,
  fsel_age,
  maxage,
  mat,
  a = seq(1, nrow(LAK)) - 1,
  fsel_type,
  fsel_len
)
```

## Arguments

- x:

  Estimated parameters. Matrix `[3, f]`

- type, fsel_type:

  Character string to indicate the functional form of selectivity. See
  details below.

- maxage:

  Maximum value of the age of full selectivity

- maxL:

  Maximum value of the length of full selectivity

- sel_par:

  Matrix of parameters returned by `conv_selpar()`

- lmid:

  Midpoints of length bins for calculating selectivity at length

- fsel_len:

  Selectivity at length matrix for fleets, returned by previous call to
  `calc_sel_len()`

- sel_len:

  Selectivity at length matrix returned by `calc_sel_len()`

- LAK:

  Length-at-age probability matrix. Matrix `[a, length(lmid)]`

- sel_block:

  Integer vector. Length `length(type)`. See details below.

- mat:

  Maturity at age vector

- a:

  Integer vector of ages corresponding to the rows of `LAK` (as well as
  `mat`)

- fsel_age:

  Matrix returned by `calc_fsel_age()`

## Value

`conv_selpar()` returns a matrix of the same dimensions as `x`.

`calc_sel_len()` returns a matrix `[l, f]`, i.e.,
`[length(lmid), length(type)]`.

`calc_fsel_age()` returns a matrix `[a, f]`, i.e.,
`[a, length(sel_block)]`

`calc_isel_age()` returns a matrix `[a, i]`, i.e., `[a, length(type)]`

## Details

Options for argument `type` include:

- functional forms with respect to length:
  `"logistic_length", "dome_length"`

- functional forms with respect to age: `"logistic_age", "dome_age"`

- for surveys, an integer (`f`) to map index selectivity at age to fleet
  `f` (will be coerced to integer)

- `"SB"` to fix to maturity at age schedule

- `"B"` to fix selectivity to 1 for all ages

- `"length_x_y"` to specify selectivity to 1 between lengths `x` and `y`
  (example: `"length_20_50"`)

- `"age_x_y"` to specify selectivity to 1 between age `x` and `y`
  (example: `"age_2_4"`)

- for surveys, `"f_x_y"` uses selectivity values from fleet `f` between
  bins `x` and `y` (either length or age depending on definition of
  selectivity for `f`) (example: `"3_2_4"`)

## Converting selectivity parameters (conv_selpar)

The first row of `x` corresponds to the length or age of full
selectivity: \\\mu_f = L\_{max}/(1 + \exp(-x\_{1,f}))\\

The second row of `x` corresponds to the width of the ascending limb for
selectivity: \\\sigma_f^{asc} = \exp(x\_{2,f})\\

The third row of `x` corresponds to the width of the descending limb for
selectivity (if dome-shaped): \\\sigma_f^{des} = \exp(x\_{3,f})\\

## Length selectivity (calc_sel_len)

The selectivity at length is \$\$ s\_{\ell} = \begin{cases} \exp(\alpha)
& L\_{\ell} \< \mu_f\\ \exp(\beta) & L\_{\ell} \ge \mu_f\\ \end{cases}
\$\$ where \\ \alpha = -0.5(L\_\ell - \mu_f)^2/(\sigma_f^{asc})^2 \\ and
\\ \beta = -0.5(L\_\ell - \mu_f)^2/(\sigma_f^{des})^2 \\

## Age selectivity (calc_fsel_age)

The equivalent selectivity at age is converted from the length values
(`sel_len`) as \$\$ s_a = \sum\_\ell s\_\ell \times
\textrm{Prob}(L\_{\ell}\|a) \$\$

If selectivity is explicitly in age units, then the values are directly
calculated from parameters in `sel_par`.

Vector `sel_block` assigns the output selectivity from a different
column of the input matrix and facilitates time-varying selectivity in
time blocks. For example, `sel_block[1] <- 2` means that selectivity
values in the first column of the output is based on the second column
of the input matrices (`sel_len[, 2]` or `sel_par[, 2]`).
