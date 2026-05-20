

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
