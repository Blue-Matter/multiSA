# Make list of parameters for RTMB

Sets up the list of parameters, map of parameters (see `map` argument in
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)), and
identifies some random effects parameters based on the input data and
some user choices on model configuration.

These functions provide a template for the parameter and map setup that
can be adjusted for alternative configurations. `check_parameters()`
checks whether custom made parameter lists are of the correct dimension.

## Usage

``` r
make_parameters(
  MSAdata,
  start = list(),
  map = list(),
  est_mov = c("none", "dist_random", "gravity_fixed"),
  silent = FALSE,
  ...
)

make_map(
  p,
  MSAdata,
  map = list(),
  est_M = FALSE,
  est_h = FALSE,
  est_mat = FALSE,
  est_sdr = FALSE,
  est_mov = c("none", "dist_random", "gravity_fixed"),
  est_qfs = FALSE,
  silent = FALSE
)

check_parameters(p = list(), map, MSAdata, silent = FALSE)
```

## Arguments

- MSAdata:

  S4 data object

- start:

  An optional list of parameters. Named list of parameters with the
  associated dimensions and transformations below. Overrides default
  values created by `make_parameters()`.

- map:

  List of mapped parameters. Used by `check_parameters()` only to count
  parameters.

- est_mov:

  Character describing structure of stock movement parameters. Only used
  if map arguments for the movement parameters is NULL. See details
  below.

- silent:

  Logical, whether `make_map()` reports messages to the console

- ...:

  Various arguments for `make_map()` (could be important!)

- p:

  List of parameters, e.g., returned by `make_parameters()`

- est_M:

  Logical, estimate natural mortality? Only used if `map$log_M_s` is
  NULL

- est_h:

  Logical, estimate steepness? Only used if `map$t_h_s` is NULL

- est_mat:

  Logical, estimate maturity? Only used if `map$mat_ps` is NULL

- est_sdr:

  Logical, estimate standard deviation of recruitment deviates? Only
  used if `map$log_sdr_s` is NULL

- est_qfs:

  Logical, estimate relative catchability of stocks by each fleet? Fix
  `log_q_fs` for the first stock if `TRUE`. Only used if `map$log_q_fs`
  is NULL

## Value

`make_parameters()` returns a list of parameters (`"p"`) concatenated
with the output of `make_map()`.

`make_map()` returns a named list containing parameter mappings
(`"map"`) and a character vector of random effects (`"random"`).

`check_parameters()` invisibly returns the parameter list if no problems
are encountered.

## Parameters

Generally parameter names will have up to three components, separated by
underscores. For example, `log_M_s` represents the natural logarithm of
natural mortality, and is a vector by stock `s`.

The first component describes the transformation from the estimated
parameter space to the normal parameter space, frequently, `log` or
`logit`. Prefix `t` indicates some other custom transformation that is
described below.

Second is the parameter name, e.g., `M` for natural mortality, `rdev`
for recruitment deviates, etc.

Third is the dimension of the parameter variable and the indexing for
the vectors, matrices, and arrays, e.g., `y` for year, `s` for stock.
See
[MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md).
Here, an additional index `p` represents some other number of parameters
that is described below.

- `t_R0_s`:

  Vector by `s`. Unfished recruitment, i.e., intersection of unfished
  replacement line and stock recruit function, is represented as:
  `R0_s <- exp(t_R0_s) * MSAdata@Dmodel@scale_s`. By default,
  `t_R0_s = 3`

- `t_h_s`:

  Vector by `s`. Steepness of the stock-recruit function. Logit space
  for Beverton-Holt (`h = 0.8 * plogis(t_h_s) + 0.2`) and log space for
  Ricker function (`h = exp(t_h_s) + 0.2`). Default steepness value of
  0.8

- `mat_ps`:

  Matrix `[2, s]`. Maturity parameters (can be estimated or specified in
  data object). Logistic functional form. The parameter in the first row
  is the age of 50 percent maturity in logit space:
  `a50_s <- plogis(mat_ps[1, ] * na)`. In the second row is the age of
  95 percent maturity as a logarithmic offset:
  `a95_s <- a50_s + exp(mat_ps[2, ])`. Default `a50_s <- 0.5 * na` and
  `a95_s <- a50_s + 1`

- `log_M_s`:

  Vector by `s`. Natural logarithm of natural mortality (can be
  estimated or specified in data object). Default parameter value for
  all stocks: `M <- -log(0.05)/MSAdata@Dmodel@na`

- `log_rdev_ys`:

  Matrix `[y, s]`. Log recruitment deviations. By default, all start
  values are at zero.

- `log_sdr_s`:

  Vector by `s`. log-Standard deviation of the log recruitment
  deviations. Default SD = 0.4

- `log_q_fs`:

  Matrix `[f, s]`. The natural logarithm of `q_fs`, the relative fishing
  efficiency of `f` for stock `s`. Equal values imply equal catchability
  of all stocks. See equations in
  [`calc_F()`](https://blue-matter.github.io/multiSA/reference/calc_F.md).
  Default sets all values to zero.

- `log_Fdev_ymfr`:

  Array `[y, m, f, r]`. Fishing mortality parameters. For each fleet,
  the log of F is estimated directly for the reference year, season,
  region. For other strata, F is an offset from this value: \$\$
  F\_{y,m,f,r} = \begin{cases} \exp(x^{\textrm{Fmult}}\_f) \quad & y =
  y\_{\textrm{ref}}, m = m\_{\textrm{ref}}, r = r\_{\textrm{ref}}\\
  \exp(x^{\textrm{Fmult}}\_f + x^{\textrm{Fdev}}\_{y,m,r}) \quad &
  \textrm{otherwise} \end{cases} \$\$

- `sel_pf`:

  Matrix `[3, f]`. Fishery selectivity parameters in logit or log space.
  See equations in
  [`conv_selpar()`](https://blue-matter.github.io/multiSA/reference/conv_selpar.md),
  where `sel_pf` is the `x` matrix.

- `sel_pi`:

  Matrix `[3, i]`. Index selectivity parameters in logit or log space.
  See equations in
  [`conv_selpar()`](https://blue-matter.github.io/multiSA/reference/conv_selpar.md),
  where `sel_pi` is the `x` matrix.

- `mov_x_marrs`:

  Array `[m, a, r, r, s]`. Base movement matrix. Set to -1000 to
  effectively exclude movements from region pairs. See equations in
  [`conv_mov()`](https://blue-matter.github.io/multiSA/reference/conv_mov.md)

- `mov_g_ymars`:

  Array `[y, m, a, r, s]`. Attractivity term in gravity model for
  movement. If `x` and `v` are zero, this matrix specifies the
  distribution of total stock abundance into the various regions. See
  equations in
  [`conv_mov()`](https://blue-matter.github.io/multiSA/reference/conv_mov.md)

- `mov_v_ymas`:

  Array `[y, m, a, s]`. Viscosity term in gravity model for movement.
  See equations in
  [`conv_mov()`](https://blue-matter.github.io/multiSA/reference/conv_mov.md)

- `log_sdg_rs`:

  Array `[r, s]`. Marginal log standard deviation in the stock
  distribution (`mov_g_ymars`) among regions for stock `s`. Only used
  when `est_mov = "dist_random"`. Default SD of 0.1.

- `t_corg_ps`:

  Array `[sum(1:(nr - 1)), s]`. Lower triangle of the correlation matrix
  for `mov_g_ymars`, to be obtained with the Cholesky factorization.
  Only used when `est_mov = dist_random`. Default values of zero.

- `log_initF_mfr`:

  Array `[m, f, r]`. Initial F corresponding to the equilibrium catch.
  Default is `log(0.1)` for seasons, fleets, and regions where
  `Dfishery@Cinit_mfr > 1e-8`, otherwise, not estimated (Finit = 0).

- `log_initrdev_as`:

  Array `[na - 1, s]`. Recruitment deviations for the initial
  abundance-at-age vector. Default is zero.

## Start list

Users can provide `R0_s` and `h_s` in the start list.
`make_parameters()` will make the appropriate transformation for the
starting values of `t_R0_s` and `t_h_s`, respectively.

## Movement setup for `make_map()`

If a single region model or `est_mov = "none"`: no movement parameters
are estimated.

If `est_mov = "dist_random"`: fix all values for `mov_x_marrs` and
`mov_v_ymas`. Fix `mov_g_ymars` for the first region for each year,
season, age, and stock. `mov_g_ymars` are random effects.

If `est_mov = "gravity_fixed"`: fix all values for `mov_x_marrs`. Fix
`mov_g_ymars` for the first region for each year, season, age, and
stock. Estimate all `mov_v_ymas`. Both `mov_g_ymars` and `mov_v_ymas`
are fixed effects.

By default `p$mov_x_marrs` is zero. Set to -1000 for areas for which
there is no abundance of a particular stock.

## See also

[MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
