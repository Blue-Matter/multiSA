# Changelog

## Version 0.5.0

- Update check for catch and lambda factors in predicted stock
  composition
- Update selectivity names, using ‘total’ and ‘mature’ instead of ‘B’
  and ‘SB’, respectively. For good practice, units are not necessarily
  in biomass.
- Add option to fix index catchability in `Dsurvey@qest_i`
- Update
  [`plot_index()`](https://blue-matter.github.io/multiSA/reference/plot-MSA-data.md)
  to make facetted plots and
  [`plot_CAL()`](https://blue-matter.github.io/multiSA/reference/plot-MSA-data.md)
  to plot aggregate length composition

## Version 0.4.0

CRAN release: 2026-06-27

- Update
  [`calc_eqdist()`](https://blue-matter.github.io/multiSA/reference/calc_eqdist.md)
  where movement is indexed within the same season as stock distribution
  (previously applied to previous seasonal time step)

## Version 0.3.0

CRAN release: 2026-06-22

- Use the potentially faster default pipe `|>` instead of magittr’s pipe
  `%>%`
- Fix movement indexing in
  [`calc_population()`](https://blue-matter.github.io/multiSA/reference/calc_population.md)
  (model previously had a 1 season lag that was erroneous)
- Initial rec devs are length `na` if advance age after season 1
  (obviously in seasonal models), otherwise remains length `na-1`
- Remove some loops with [`apply()`](https://rdrr.io/r/base/apply.html)
  to speed up
  [`calc_F()`](https://blue-matter.github.io/multiSA/reference/calc_F.md)
- Length-age matrices, when modified by selectivity have a tiny number
  added to denominator to avoid division by zero
- Update how selectivity arrays are filled in, fixes issue when time
  blocks are used
- Use parallel package functions instead of snowfall for parallel
  computation with profiling and retrospectives
- Fix indexing with predictions of stock composition by coercing a
  vector of 1’s and 0’s to boolean
- Add
  [`do_jitter()`](https://blue-matter.github.io/multiSA/reference/do_jitter.md)
  function
- Clean up various internal functions:
  [`like_comp()`](https://blue-matter.github.io/multiSA/reference/like_comp.md),
  [`optimize_RTMB()`](https://blue-matter.github.io/multiSA/reference/optimize_RTMB.md),
  [`get_sdreport()`](https://blue-matter.github.io/multiSA/reference/get_sdreport.md)
- Export
  [`get_likelihood_components()`](https://blue-matter.github.io/multiSA/reference/profile.md)
- Reporet fitting time in
  [`fit_MSA()`](https://blue-matter.github.io/multiSA/reference/fit_MSA.md)

## Version 0.2.0

CRAN release: 2026-05-22

- New selectivity options: constant over size and age range, mapping a
  subset of length or age from fleet to index
- Add
  [`calc_init_population()`](https://blue-matter.github.io/multiSA/reference/calc_init_population.md)
  for spool-up in spatial or seasonal models
- Report most state variables and fits to data invisibly in plotting
  function
- Add multivariate logitnormal likelihood to comp data
- Update residual calculation for composition data
- Fix predictions of tag transitions. Model previously had a 1 season
  lag that was erroneous.

## Version 0.1.1

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

## Version 0.1.0

CRAN release: 2026-02-03

- Initial CRAN release
