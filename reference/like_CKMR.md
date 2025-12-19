# Likelihood for CKMR

Returns the log-likelihood for a set of pairwise comparisons. For a
parent-offspring pair, a comparison is defined by the capture year of
parent, capture age of parent, and birth year of offspring. For a
half-sibling pair, a comparison is defined by the cohort year of each
sibling. Binomial and Poisson distributions are supported (Conn et al.
2020).

## Usage

``` r
like_CKMR(n, m, p, type = c("binomial", "poisson"))
```

## Arguments

- n:

  The number of pairwise comparisons

- m:

  The number of kinship matches

- p:

  The probability of kinship match

- type:

  The statistical distribution for the likelihood calculation

## Value

Numeric representing the log-likelihood.

## References

Conn, P.B. et al. 2020. Robustness of close-kin mark-recapture
estimators to dispersal limitation and spatially varying sampling
probabilities. Ecol. Evol. 10: 5558-5569.
[doi:10.1002/ece3.6296](https://doi.org/10.1002/ece3.6296)

## See also

[`calc_POP()`](https://blue-matter.github.io/multiSA/reference/calc_POP.md)
[`calc_HSP()`](https://blue-matter.github.io/multiSA/reference/calc_POP.md)
