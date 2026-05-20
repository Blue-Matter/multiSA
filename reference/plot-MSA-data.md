# Plotting functions for data in MSA model

A set of functions to plot data variables and predicted values (catch,
age composition, etc.)

## Usage

``` r
plot_catch(
  fit,
  f = 1,
  by = c("region", "stock"),
  prop = FALSE,
  annual = FALSE,
  figure = TRUE
)

plot_index(fit, i = 1, zoom = FALSE, figure = TRUE)

plot_CAA(fit, f = 1, r = 1, do_mean = FALSE, figure = TRUE)

plot_CAL(fit, f = 1, r = 1, do_mean = FALSE, figure = TRUE)

plot_IAA(fit, i = 1, do_mean = FALSE, figure = TRUE)

plot_IAL(fit, i = 1, do_mean = FALSE, figure = TRUE)

plot_SC(fit, ff = 1, aa = 1, r = 1, prop = FALSE, figure = TRUE)

plot_tagmov(fit, s = 1, yy = 1, aa = 1, figure = TRUE)
```

## Arguments

- fit:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- f:

  Integer, indexes the fleet

- by:

  Character to indicate dimension for multivariate plots

- prop:

  Logical, whether to plot proportions (TRUE) or absolute numbers

- annual:

  Logical, whether to plot annual values (summed over seasons)

- figure:

  Logical, whether to generate the plot

- i:

  Integer, indexes the survey

- zoom:

  Logical, for `plot_index()`. If `TRUE`, plots a subset of years with
  observed data points. Otherwise, plots predicted values over all model
  years.

- r:

  Integer, indexes the region

- do_mean:

  Logical, whether to plot full compositions or time series of mean
  length or mean age

- ff:

  Integer, indexes the aggregate fleet (for stock composition data)

- aa:

  Integer, indexes the aggregate age class (for stock composition and
  tag data)

- s:

  Integer, indexes the stock

- yy:

  Integer, indexes the aggregate years (for the tag data)

## Value

Invisible data frame of observed and predicted values plotted in base
graphics figures

## Details

- `plot_catch` plots the fishery catch by stock or region (either whole
  numbers or proportions)

&nbsp;

- `plot_index` plots indices of abundance

&nbsp;

- `plot_CAA` plots the fishery catch at age

&nbsp;

- `plot_CAL` plots the catch at length

&nbsp;

- `plot_IAA` plots the index age composition

&nbsp;

- `plot_IAL` plots the index length composition

&nbsp;

- `plot_SC` plots the stock composition

&nbsp;

- `plot_tagmov` plots the tag movements
