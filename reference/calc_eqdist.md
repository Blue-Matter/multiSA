# Equilibrium distribution from movement matrix

Applies the seasonal movement matrices several times in order to obtain
the equilibrium spatial distribution over the course of a year. Not used
in the model but useful for reporting.

## Usage

``` r
calc_eqdist(
  x,
  nm = dim(x)[1],
  nr = dim(x)[2],
  start = rep(1/nr, nr),
  m_start = 1L,
  nit = 20
)
```

## Arguments

- x:

  Movement array `[m, r, r]`. The second dimension corresponds to origin
  (sum to 1), and third dimension corresponds to destination

- nm:

  Number of seasons

- nr:

  Number of regions

- start:

  The initial distribution. Vector of length `nr`

- m_start:

  Integer, the season in which to apply the initial distribution and
  start the projection

- nit:

  Integer, the number of times the movement matrix will be applied

## Value

Matrix by season and region `[m, r]`
