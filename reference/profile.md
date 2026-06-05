# Profile parameters of MSA model

Evaluate change in objective function and likelihood components for up
to 2 parameters.

## Usage

``` r
# S4 method for class 'MSAassess'
profile(
  fitted,
  p1,
  v1,
  p2,
  v2,
  use_fitted = TRUE,
  return_models = TRUE,
  cores = 1,
  ...
)

get_likelihood_components(fitted)

# S3 method for class 'MSAprof'
plot(
  x,
  component = "objective",
  rel = TRUE,
  xlab,
  ylab,
  main,
  plot2d = c("contour", "filled.contour"),
  ...
)
```

## Arguments

- fitted:

  [MSAassess](https://blue-matter.github.io/multiSA/reference/MSAassess-class.md)
  object returned by
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

- p1:

  Character string that represents the first parameter to be profiled,
  including the parameter name and index of the vector/array. See
  "Parameters" section of
  [`make_parameters()`](https://blue-matter.github.io/multiSA/reference/make_parameters.md).
  Additionally, this function allows users to specify `R0_s` and `h_s`
  (in normal units).

- v1:

  Vector of values corresponding to `p1`

- p2:

  Character string that represents the optional second parameter to be
  profiled

- v2:

  Vector of values corresponding to `p2`

- use_fitted:

  Logical, whether to use estimated parameters or the starting values in
  `fitted` to start the profile

- return_models:

  Logical, whether to return fitted models in the profile

- cores:

  Integer for the number of cores to use for parallel processing

- ...:

  Other argument to the base graphics function, i.e., either plot() or
  contour()

- x:

  Output from `profile.MSAassess()`

- component:

  Character for the column in `x` to be plotted

- rel:

  Logical, whether the relative change in `component` is plotted (TRUE)
  or the raw values (FALSE)

- xlab:

  Optional character for the x-axis label

- ylab:

  Optional character for the y-axis label

- main:

  Optional character for the plot title

- plot2d:

  Character, plotting function for two-dimensional profiling (either a
  [`contour()`](https://rdrr.io/r/graphics/contour.html) or
  [`filled.contour()`](https://rdrr.io/r/graphics/filled.contour.html)
  plot)

## Value

Named list (length 2).

First, `profile` contains a data frame of the likelihood values that
correspond to fixed values of `p1` and `p2`. Other columns:

- Likelihood `loglike` refers to maximizing the probability of the
  observed data (higher values for better fit)

- Prior `logprior` refers to maximizing the probability of a parameter
  to their prior distribution (higher values are closer to the prior
  mode)

- Penalty `penalty` are values added to the objective function when
  parameters exceed model bounds (lower values are better)

- `fn` is the objective function returned by RTMB (lower values are
  better)

- `objective` is the objective function returned by the optimizer (lower
  values are better)

Second, `fits` contains a list of the `MSAassess` objects if
`return_models = TRUE`.

`get_likelihood_components()` returns a data.frame of the components to
the objective function (log-likelihoods, log-priors, etc.)

The accompanying plot function returns a line plot for a 1-dimensional
profile or a contour plot for a two dimensional profile. Will plot the
negative log likelihood or negative log prior (better fit with lower
values).

Relative values are obtained by subtracting from the fitted value. See
`attr(x$profile, "fitted")`
