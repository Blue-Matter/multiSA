# Calculate index at age

For indices of abundance, the function calculates the numbers vulnerable
to the survey.

## Usage

``` r
calc_index(
  N,
  Z,
  sel,
  na = dim(N)[1],
  nr = dim(N)[2],
  ns = dim(N)[3],
  ni = dim(sel)[2],
  samp = array(1, c(ni, nr, ns)),
  delta = rep(0, ni)
)
```

## Arguments

- N:

  Stock abundance at the beginning of the time step. Array `[a, r, s]`

- Z:

  Instantaneous total mortality. Array `[a, r, s]`

- sel:

  Index selectivity. Array `[a, i, s]`

- na:

  Integer, number of age classes

- nr:

  Integer, number of regions

- ns:

  Integer, number of stocks

- ni:

  Integer, number of indices

- samp:

  Boolean indicates which regions and stocks are sampled by the index.
  Array `[i, r, s]`

- delta:

  Fraction of time step when the index samples the population. Vector by
  `i`. Set to a negative number (-1) to sample the average population
  over the course of the time step, i.e. `N * (1 - exp(-Z))/Z`.

## Value

Index at age. Array `[a, i, s]`

## Details

The index is calculated as \$\$ I\_{a,i,s} = v\_{a,i,s} \sum_r
N\_{a,r,s} d\_{a,r,s} \times \mathbb{1}(r \in R_i) \mathbb{1}(s \in S_i)
\$\$

If the survey samples at a moment in time, then \$\$ d\_{a,r,s} =
\exp(-\delta_i Z\_{a,r,s}) \$\$

Otherwise, if the index samples the population over the duration of the
time step, then \$\$ d\_{a,r,s} = (1 - \exp(-Z\_{a,r,s}))/Z\_{a,r,s}
\$\$

where \\R_i\\ and \\S_i\\ denote the regions and stocks, respectively,
sampled by index \\i\\. For example, \\R_2 = 1\\ denotes that the second
index of abundance only samples region 1. These are informed by array
`samp` where `samp[i, r, s] = 1` indicates that stock `s` in region `r`
is sampled by index `i`.
