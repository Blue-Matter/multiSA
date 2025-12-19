# Priors for MSA model

Priors in MSA are set by providing character strings which are then
parsed into an expression and evaluated in the model environment (see
example). This provides flexibility to set a prior for any desired model
parameter or variable. See list of parameters in the documentation for
`[check_parameters()]` for options (note that priors for `log_rdev_ys`
and `log_initrdev_as` are not needed as they're hard-coded into the
model). Several functions below generate the character string for the
prior for important dynamics parameters, such as natural mortality and
steepness.

## Usage

``` r
prior_h(MSAdata, s = 1, m, stdev)

prior_M(MSAdata, s = 1, meanlog, sdlog)

prior_q(MSAdata, i = 1, meanlog, sdlog)
```

## Arguments

- MSAdata:

  Data object. Class
  [MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)

- s:

  Integer for stock

- m:

  Mean in un-transformed space

- stdev:

  Standard deviation in un-transformed space

- meanlog:

  Mean of the lognormal distribution on the log scale

- sdlog:

  Standard of the lognormal distribution on the log scale

- i:

  Integer for the corresponding index

## Value

Character.

## Details

- `prior_h` returns the log prior for stock-recruit steepness. Beta
  distribution for the Beverton-Holt function and normal distribution
  for Ricker function.

&nbsp;

- `prior_M` returns the log prior for natural mortality. Lognormal
  distribution.

&nbsp;

- `prior_q` returns the log prior for index catchability. Lognormal
  distribution.

## Examples

``` r
# Add M and steepness prior to model

dat <- new("MSAdata")
dat@Dmodel@ns <- 1
dat@Dstock@SRR_s <- "BH"

pr_M <- prior_M(dat, s = 1, log(0.25), 0.3)
pr_h <- prior_h(dat, s = 1, 0.8, 0.15)
dat@Dmodel@prior <- c(pr_M, pr_h)
```
