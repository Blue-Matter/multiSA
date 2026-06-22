
## multiSA 0.3.0

- Use the potentially faster default pipe `|>` instead of magittr's pipe `%>%`
- Fix movement indexing in `calc_population()` (model previously had a 1 season lag that was erroneous)
- Initial rec devs are length `na` if advance age after season 1 (obviously in seasonal models), otherwise remains length `na-1`
- Remove some loops with `apply()` to speed up `calc_F()`
- Length-age matrices, when modified by selectivity have a tiny number added to denominator to avoid division by zero
- Update how selectivity arrays are filled in, fixes issue when time blocks are used
- Use parallel package functions instead of snowfall for parallel computation with profiling and retrospectives
- Fix indexing with predictions of stock composition by coercing a vector of 1's and 0's to boolean
- Add `do_jitter()` function
- Clean up various internal functions: `like_comp()`, `optimize_RTMB()`, `get_sdreport()`
- Export `get_likelihood_components()`
- Reporet fitting time in `fit_MSA()`

## multiSA 0.2.0

- New selectivity options: constant over size and age range, mapping a subset of length or age from fleet to index
- Add `calc_init_population()` for spool-up in spatial or seasonal models
- Report most state variables and fits to data invisibly in plotting function
- Add multivariate logitnormal likelihood to comp data
- Update residual calculation for composition data
- Fix predictions of tag transitions. Model previously had a 1 season lag that was erroneous. 

## multiSA 0.1.1

- More robust max F check when `calc_F()` to prevent numerical overflow (check in log space rather than normal space)
- Various checks for NA's in plotting functions
- Fix internal function `collapse_yearseason` that converts year and season dimensions of arrays into single time dimension 
- Various fixes for cleaner console reporting
- Profile function exports model object and parameter values to individual cores in parallel mode
- Properly dispatch S4 generics for `MSAassess` (previously used S3 methods)
- Remove `stats::uniroot` import for compatibility with RTMB 1.9

## multiSA 0.1.0 (alpha version)

- Initial CRAN release
