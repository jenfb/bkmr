# bkmr 0.2.0

## Major changes

* Added ability to have binomial outcome `family` by implementing probit regression within `kmbayes()`

* Removed computation of the subject-specific effects `h[i]` within `kmbayes()`, as this is not always desired, and greatly slows down model fitting

  * This could still be done by setting the option `est.h = TRUE` in the `kmbayes` function
  
  * posterior samples of `h[i]` can now be obtained via the post-processing `SamplePred` function; alternatively, posterior summaries (mean, variance) can be obtained via the post-processing `ComputePostmeanHnew` function

* Added ability to use exact estimates of the posterior mean and variance by specifying the argument `method = 'exact'` within the post-processing functions (e.g., `OverallRiskSummaries()`, `PredictorResponseUnivar()`)

## Bug fixes

* Fixed `PredictorResponseBivarLevels()` when argument `both_pairs = TRUE` (#4)