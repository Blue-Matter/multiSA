# Calculate von Bertalanffy length-at-age

Returns an array of length-at-age with seasonal dimension. Useful for
[MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
inputs.

## Usage

``` r
calc_growth(
  Linf_s,
  K_s,
  t0_s,
  ns = length(Linf_s),
  nm = 4,
  ny = 20,
  a = seq(1, 10) - 1
)
```

## Arguments

- Linf_s:

  Vector by stock `s` of asymptotic length

- K_s:

  Vector by stock `s` of the growth coefficient

- t0_s:

  Vector by `s` of the age at length zero.

- ns:

  Integer, number of stocks

- nm:

  Integer, number of seasons

- ny:

  Integer, number of years

- a:

  Integer vector of ages

## Value

Array `[y, m, a, s]`

## Examples

``` r
len_ymas <- calc_growth(c(30, 40), c(0.4, 0.2), c(-1, -1))

# Calculate stock weight at age
a_s <- rep(1e-6, 2)
b_s <- c(3, 3.1)

ns <- length(a_s)
swt_ymas <- sapply(1:ns, function(s) {
  a_s[s]*len_ymas[, , , s]^b_s[s]
}, simplify = "array")
```
