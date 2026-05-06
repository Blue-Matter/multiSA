# Changelog

## multiSA 0.2.0

- New selectivity options (constant over size and age range), mapping a
  subset of length or age from fleet to index

## multiSA 0.1.1

CRAN release: 2026-03-20

- More robust max F check when
  [`calc_F()`](https://blue-matter.github.io/multiSA/reference/calc_F.md)
  to prevent numerical overflow (check in log space rather than normal
  space)
- Various checks for NA’s in plotting functions
- Fix internal function `collapse_yearseason` that converts year and
  season dimensions of arrays into single time dimension
- Various fixes for cleaner console reporting
- Profile function exports model object and parameter values to
  individual cores in parallel mode
- Properly dispatch S4 generics for `MSAassess` (previously used S3
  methods)
- Remove [`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html) import
  for compatibility with RTMB 1.9

## multiSA 0.1.0 (alpha version)

CRAN release: 2026-02-03

- Initial CRAN release
