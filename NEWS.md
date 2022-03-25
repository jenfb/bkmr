# bkmr 0.2.2

## Bug fixes

* Corrected code that produced warning `length > 1 in coercion to logical`

* Update functions that use deprecated functions from `dplyr` package 

## Minor changes

* No longer export the following functions:

  * `CalcGroupPIPs`, `CalcWithinGroupPIPs`, and `CalcPIPs` as these should typically be calculated using the function `ExtractPIPs`

  * `ComputePostmeanHnew.approx` and `ComputePostmeanHnew.exact` as these should typically be calculated using the function `ComputePostmeanHnew`
  
  * `set_verbose_opts` as this is only called internally

* Expanded function documentation by adding example code

# bkmr 0.2.1

## Bug fixes

* allowable values for starting parameter for `r[m]` parameters updated as follows

  * no longer truncated to a single value (when `varsel = FALSE` and `rmethod = "varying"`)

  * can be equal to 0 (when `varsel = TRUE`)

* Error no longer generated if starting values for h.hat are not positive 

* When checking class of an object, use `inherits()` instead of `class()`

# bkmr 0.2.0

## Major changes

* Added ability to have binomial outcome `family` by implementing probit regression within `kmbayes()`

* Removed computation of the subject-specific effects `h[i]` within `kmbayes()`, as this is not always desired, and greatly slows down model fitting

  * This could still be done by setting the option `est.h = TRUE` in the `kmbayes` function
  
  * posterior samples of `h[i]` can now be obtained via the post-processing `SamplePred` function; alternatively, posterior summaries (mean, variance) can be obtained via the post-processing `ComputePostmeanHnew` function

* Added ability to use exact estimates of the posterior mean and variance by specifying the argument `method = 'exact'` within the post-processing functions (e.g., `OverallRiskSummaries()`, `PredictorResponseUnivar()`)

## Bug fixes

* Fixed `PredictorResponseBivarLevels()` when argument `both_pairs = TRUE` (#4)