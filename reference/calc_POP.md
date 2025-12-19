# Predict the probability of CKMR kinship pairs

Calculate the probability of observing a parent-offspring pair
(`calc_POP`) and half-sibling pair (`calc_HSP`) for closed-kin mark
recapture (CKMR) for an age-structured model.

## Usage

``` r
calc_POP(t, a, y, N, fec)

calc_HSP(yi, yj, N, fec, Z)
```

## Arguments

- t:

  Vector, capture year of parent `i`

- a:

  Vector, age at capture of parent `i`

- y:

  Vector, birth year of offspring `j`

- N:

  Abundance of mature spawners. Matrix by `[y, a]`

- fec:

  Fecundity schedule of mature spawners. Matrix by `[y, a]`

- yi:

  Vector, birth year of sibling `i`. Must be older than sibling `j`.

- yj:

  Vector, birth year of sibling `j`.

- Z:

  Instantaneous total mortality rate. Matrix by `[y, a]`

## Value

Numeric, vector of probabilities

## Parent-offspring pairs

The parent-offspring probability is calculated from Bravington et al.
2016, eq 3.4:

\$\$p\_{\textrm{POP}} = 2 \times \dfrac{f(y_j,y_j - (t_i - a_i))}{\sum_a
f(y_j,a) N(y_j,a)}\$\$

where \\y_j - (t_i - a_i)\\ is the parental age in year \\y_j\\. Scalar
2 accounts for the fact that the parent could be either a mother or a
father. `calc_POP` is vectorized with respect to `t`, `a`, and `y`.

## Half-sibling pairs

The half-sibling probability is calculated from Bravington et al. 2016,
eq 3.10, and expanded by Hillary et al. 2018, Supplement S2.8.1 for
age-specific survival and fecundity of the parent:

\$\$p\_{\textrm{HSP}} = 4 \times \sum_a\left( \dfrac{N(y_i, a)f(y_i,
a)}{\sum\_{a'} N(y_i, a')f(y_i,a')}\times \exp(-\sum\_{t = 0}^{y_j -
y_i - 1} Z(y_i + t,a + t))\times \dfrac{f(y_j,a+y_j-y_i)}{\sum\_{a'}
N(y_j,a')f(y_j,a')} \right) \$\$

- The first ratio is the probability that a fish at age \\a\\ in year
  \\y_i\\ is the parent of \\i\\.

- The exponential term is that fish's survival from year \\y_i\\ to
  \\y_j\\.

- The second ratio is the probability that the parent of \\i\\, age
  \\a+y_j-y_i\\ in year \\y_j\\, is the parent of \\j\\.

The parent is not observed in the HSP, so we sum the probabilities over
all potential ages in year \\y_i\\. `calc_HSP` is vectorized with
respect to `yi` and `yj`.

## References

Bravington, M.V. et al. 2016. Close-Kin Mark-Recapture. Stat. Sci. 31:
259-274. [doi:10.1214/16-STS552](https://doi.org/10.1214/16-STS552)

Hillary, R.M. et al. 2018. Genetic relatedness reveals total population
size of white sharks in eastern Australia and New Zealand. Sci. Rep. 8:
2661.
[doi:10.1038/s41598-018-20593-w](https://doi.org/10.1038/s41598-018-20593-w)

## See also

[`like_CKMR()`](https://blue-matter.github.io/multiSA/reference/like_CKMR.md)

## Author

Q. Huynh with contribution from Y. Tsukahara (Fisheries Research
Institute, Japan)
