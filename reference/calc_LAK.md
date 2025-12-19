# Length-at-age key

Calculates the probability distribution of length-at-age using the
normal probability density function

## Usage

``` r
calc_LAK(len_a, sd_la, lbin, nl = length(lbin))
```

## Arguments

- len_a:

  Vector of length-at-age

- sd_la:

  Vector of standard deviation in length-at-age

- lbin:

  Vector of the lower boundary of the length bins

- nl:

  Integer, number of length bins (default is `length(lbin)`)

## Value

Matrix by age (rows) and length (columns)
