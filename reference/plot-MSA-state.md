# Plotting functions for fitted MSA model

A set of functions to plot state variables (biomass, recruitment time
series, etc.)

## Usage

``` r
plot_S(
  fit,
  by = c("stock", "region"),
  r,
  s,
  prop = FALSE,
  facet_free = FALSE,
  figure = TRUE,
  ylab
)

plot_B(
  fit,
  by = c("stock", "region"),
  r,
  s,
  prop = FALSE,
  facet_free = FALSE,
  figure = TRUE
)

plot_R(fit, s, figure = TRUE)

plot_SRR(fit, s = 1, phi = TRUE)

plot_Rdev(fit, s = 1, log = TRUE, figure = TRUE)

plot_Fstock(fit, s, by = c("annual", "season"), figure = TRUE)

plot_self(fit, f = 1, type = c("length", "age"), figure = TRUE)

plot_seli(fit, i = 1, figure = TRUE)

plot_selstock(
  fit,
  s = 1,
  by = c("annual", "season"),
  plot2d = c("contour", "filled.contour"),
  ...
)

plot_N(fit, m = 1, r, s = 1, plot2d = c("contour", "filled.contour"), ...)

plot_V(fit, f = 1, by = c("stock", "region"), prop = FALSE, facet_free = FALSE)

plot_Ffleet(fit, f = 1)

plot_mov(fit, s = 1, y, a, palette = "Peach", figure = TRUE)

plot_recdist(fit, palette = "Peach", figure = TRUE)
```

## Arguments

- fit:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- by:

  Character to indicate whether to calculate selectivity from F per year
  or per season

- r:

  Integer for the corresponding region

- s:

  Integer for the corresponding stock

- prop:

  Logical, whether to plot proportions (TRUE) or absolute numbers

- facet_free:

  Logical, whether to allow the y-axis limits to vary by panel in
  facetted plots

- figure, :

  Logical, whether to generate the plot

- ylab:

  Optional character string for custom y-axis label

- phi:

  Logical, whether to plot unfished replacement line

- log:

  Logical, whether to plot the natural logarithm of the response
  variable

- f:

  Integer for the corresponding fleet

- type:

  For `plot_self`, indicates whether to plot the selectivity by age or
  length.

- i:

  Integer for the corresponding survey

- plot2d:

  Character, plotting function for either a
  [`contour()`](https://rdrr.io/r/graphics/contour.html) or
  [`filled.contour()`](https://rdrr.io/r/graphics/filled.contour.html)
  plot

- ...:

  Other argument to the base graphics function

- m:

  Integer for the corresponding season

- y:

  Integer, year for plotting the movement matrix (last model year is the
  default)

- a:

  Integer, corresponding age for plotting the movement matrix (age 1 is
  the default)

- palette:

  Character, palette name to pass to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).
  See
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  for options.

## Value

Invisible data frame of state variables that were plotted in base
graphics figures

## Details

- `plot_S` plots spawning output by stock or region (either whole
  numbers or proportions for the latter)

&nbsp;

- `plot_B` plots total biomass by stock or region (either whole numbers
  or proportions for the latter)

&nbsp;

- `plot_R` plots recruitment by stock

&nbsp;

- `plot_SRR` plots the stock-recruitment relationship and history
  (realized recruitment) by stock

&nbsp;

- `plot_Rdev` plots recruitment deviations by stock

&nbsp;

- `plot_Fstock` plots apical instantaneous fishing mortality (per year
  or per season) by stock

&nbsp;

- `plot_self` plots fishery selectivity

&nbsp;

- `plot_seli` plots index selectivity

&nbsp;

- `plot_selstock` plots the realized selectivity from total catch and
  total abundance at age

&nbsp;

- `plot_N` reports total abundance at age

&nbsp;

- `plot_V` plots vulnerable biomass, availability to the fishery

&nbsp;

- `plot_Ffleet` plots apical instantaneous fishing mortality (per
  season) by fleet

&nbsp;

- `plot_mov` plots movement matrices and the corresponding equilibrium
  distribution in multi-area models

&nbsp;

- `plot_recdist` plots the distribution of recruitment for each stock
