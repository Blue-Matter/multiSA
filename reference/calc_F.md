# Newton-Raphson search for fishing mortality

Performs a root finding routine to find the index of F that minimizes
the difference between observed catch and the value predicted by the
Baranov equation.

## Usage

``` r
calc_F(
  Cobs,
  N,
  sel,
  wt,
  M,
  q_fs,
  delta = 1,
  na = dim(N)[1],
  nr = dim(N)[2],
  ns = dim(N)[3],
  nf = length(Cobs),
  Fmax = 2,
  nitF = 5L,
  trans = c("log", "logit")
)
```

## Arguments

- Cobs:

  Observed catch. Matrix `[f, r]`

- N:

  Stock abundance at the beginning of the time step. Array `[a, r, s]`

- sel:

  Selectivity. Array `[a, f, s]`

- wt:

  Fishery weight at age. Array `[a, f, s]`

- M:

  Instantaneous natural mortality. Units of per year `[a, s]`

- q_fs:

  Relative catchability of stock `s` for fleet `f`. Defaults to 1 if
  missing. Matrix `[f, s]`

- delta:

  Numeric, the duration of time in years corresponding to the observed
  catch, e.g., 0.25 is a quarterly time step.

- na:

  Integer, number of age classes

- nr:

  Integer, number of regions

- ns:

  Integer, number of stocks

- nf:

  Integer, number of fleets

- Fmax:

  Numeric, the maximum Findex value

- nitF:

  Integer, number of iterations for the Newton-Raphson routine

- trans:

  Whether to perform the search in log or logit space

## Value

A list containing:

- `F_afrs` Fishing mortality array

- `F_ars` Fishing mortality array (summed across fleets)

- `Z_ars` Total mortality array

- `F_index` Index of fishing mortality. Matrix `[f, r]`

- `CB_frs` Catch (biomass) array

- `CN_afrs` Catch (abundance) array

- `VB_afrs` Vulnerable biomass at the beginning of the time step. Array

- `penalty` Penalty term returned by
  [`posfun()`](https://blue-matter.github.io/multiSA/reference/posfun.md)
  when `F_index` exceeds `Fmax`

- `fn` Difference between predicted and observed catch at the last
  iteration. Matrix `[f, r]`

- `gr` Gradient of `fn` with respect to `F_index` in either log or logit
  space at the last iteration. Vector by `[f, r]`

## Details

The predicted catch for fleet `f` in region `r` is \$\$
C^{\textrm{pred}}\_{f,r} = \sum_s \sum_a v\_{a,f,s} q\_{f,s} F\_{f,r}
\dfrac{1 - \exp(-Z\_{a,r,s})}{Z\_{a,r,s}} N\_{a,r,s} w\_{a,f,s} \$\$

The Newton-Raphson routine minimizes \\f(x\_{f,r}) =
C\_{f,r}^{\textrm{pred}} - C\_{f,r}^{\textrm{obs}}\\.

If `trans = "log"`, \\F\_{f,r} = \exp(x\_{f,r})\\.

If `trans = "logit"`, \\F\_{f,r} = F\_{\textrm{max}}/(1 +
\exp(x\_{f,r}))\\.

The gradient with respect to \\\vec{x}\\ is \$\$ f'(x\_{f,r}) = \sum_s
\sum_a v\_{a,f,s} q\_{f,s} N\_{a,r,s} w\_{a,f,s}
\left(\dfrac{\alpha\gamma}{\beta}\right)' \$\$

\$\$ \left(\dfrac{\alpha\gamma}{\beta}\right)' = \dfrac{(\alpha\gamma' +
\alpha'\gamma)\beta - \alpha\gamma\beta'}{\beta^2} \$\$

where

|  |
|----|
| \\\alpha\_{f,r} = F\_{f,r}\\ |
| \\\beta\_{a,r,s} = Z\_{a,r,s} = M\_{a,s} + \sum_f v\_{a,f,s} q\_{f,s} F\_{f,r}\\ |
| \\\gamma\_{a,r,s} = 1 - \exp(-Z\_{a,r,s})\\ |
| \\\beta'\_{a,f,r,s} = v\_{a,f} q\_{f,s} \alpha'\_{f,r}\\ |
| \\\gamma'\_{a,f,r,s} = \exp(-Z\_{a,r,s})\beta'\_{a,f,r,s}\\ |

If `trans = "log"`, \\\alpha'\_{f,r} = \alpha\_{f,r}\\.

If `trans = "logit"`, \\\alpha'\_{f,r} =
F\_{\textrm{max}}\exp(-x\_{f,r})/(1 + \exp(-x\_{f,r}))^2\\.

This function solves for \\\vec{x}\\ by iterating until \\f(\vec{x})\\
approaches zero, where the vector arrow indexes over fleet and region.
In iteration \\i+1\\: \$\$\vec{x}\_{i+1} = \vec{x}\_i -
\dfrac{f(\vec{x}\_i)}{f'(\vec{x}\_i)}\$\$.

## Author

Q. Huynh
