# Likelihood for composition vectors

Returns the log-likelihood for composition data, e.g., length, age, or
stock composition, with various statistical distributions supported.

## Usage

``` r
like_comp(
  obs,
  pred,
  type = c("multinomial", "dirmult1", "dirmult2", "lognormal", "logitnormal"),
  N = sum(obs),
  theta,
  stdev
)
```

## Arguments

- obs:

  A vector of observed values. Internally converted to proportions.

- pred:

  A vector of predicted values. Same length as `obs`. Internally
  converted to proportions.

- type:

  Character for the desired distribution

- N:

  Numeric, the sample size corresponding to `obs` for multinomial or
  Dirichlet multinomial distributions.

- theta:

  Numeric, the linear (`type = "dirmult1"`) or saturating
  (`type = "dirmult2"`) Dirichlet-multinomial parameter, respectively.
  See Thorson et al. (2017)

- stdev:

  Numeric or vectorized for `obs`, the likelihood standard deviation for
  lognormal distribution.

## Value

Numeric representing the log-likelihood.

## Details

Observed and predicted vectors are internally converted to proportions.

For `type = "lognormal"` or `"logitnormal"`, zero observations are
removed from the likelihood calculation.

## References

Thorson et al. 2017. Model-based estimates of effective sample size in
stock assessment models using the Dirichlet-multinomial distribution.
Fish. Res. 192:84-93.
[doi:10.1016/j.fishres.2016.06.005](https://doi.org/10.1016/j.fishres.2016.06.005)

## Examples

``` r
M <- 0.1
age <- seq(1:10)
pred <- exp(-M * age)
obs <- pred * rlnorm(10, sd = 0.05)
like_comp(obs, pred, N = 10, type = "multinomial")
#> [1] -7.79363
like_comp(obs, pred, N = 100, type = "multinomial")
#> [1] -17.44753
like_comp(obs, pred, N = 10, type = "dirmult1", theta = 1)
#> [1] -11.33651
like_comp(obs, pred, N = 10, type = "dirmult1", theta = 20)
#> [1] -8.015671
```
