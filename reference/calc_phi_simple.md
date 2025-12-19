# Simple spawners per recruit calculation

Calculate spawners per recruit by individual stock. Appropriate for a
population model with a single spatial area and an annual time steps,
i.e. single season.

## Usage

``` r
calc_phi_simple(Z, fec_a, mat_a, delta = 0)

calc_NPR(Z, na = length(Z), plusgroup = TRUE)
```

## Arguments

- Z:

  Total mortality at age

- fec_a:

  Fecundity at age. Vector

- mat_a:

  Maturity at age. Vector

- delta:

  Fraction of season that elapses when spawning occurs, e.g., midseason
  spawning when `delta = 0.5`.

- na:

  Integer, number of age classes

- plusgroup:

  Logical, whether the largest age class is a plusgroup accumulator age

## Value

`calc_phi_simple()` returns a numeric for the spawning biomass per
recruit.

`calc_NPR()` returns a numeric, numbers per recruit at age

## See also

[`calc_phi_project()`](https://blue-matter.github.io/multiSA/reference/calc_phi_project.md)
