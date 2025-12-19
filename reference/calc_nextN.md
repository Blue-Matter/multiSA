# Project stock abundance to the next time step

This function applies survival of the current abundance, advances age
classes, re-distributes the stock using the movement matrix.

## Usage

``` r
calc_nextN(
  N,
  surv,
  na = dim(N)[1],
  nr = dim(N)[2],
  ns = dim(N)[3],
  advance_age = TRUE,
  mov = array(1/nr, c(na, nr, nr, ns)),
  plusgroup = TRUE
)
```

## Arguments

- N:

  Abundance at current time step. Array `[a, r, s]`

- surv:

  Survival during the current time step. Array `[a, r, s]`

- na:

  Integer, number of age classes

- nr:

  Integer, number of regions

- ns:

  Integer, number of stocks

- advance_age:

  Logical, whether the animals advance to their next age class

- mov:

  Movement array in the next time step. Array `[a, r, r, s]`. Rows
  denote region of origin and columns denote region of destination.

- plusgroup:

  Logical, whether the last age class is an accumulator plus group.

## Value

Abundance at the next time step. Array `[a, r, s]`
