## Test environments
* local OS X install, R-devel
* win-builder (devel)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

* New submission

  Package was archived on CRAN

  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2022-03-20 as check problems were not
      corrected in time.

  'length > 1 in coercion to logical' error in check of 'bkmrhat'.

These check problems have now been corrected.  

## Downstream dependencies
I have also run R CMD check on downstream dependencies with no issues.
