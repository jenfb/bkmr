## Resubmission
This is a resubmission. In this version I have:

* added the paper reference to the description field of the DESCRIPTION file
* added a \value field to the documentation for exported methods
* added small executable examples in Rd-files for exported functions

## Test environments
* local OS X install, R-devel
* win-builder (devel)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

* New submission

  Package was archived on CRAN

  Possibly misspelled words in DESCRIPTION:
    Bobb (9:41)
    al (9:49)
    et (9:46)
    
    * A reference for the methods in the package is included in the 
      description as requested by CRAN; these words are the author names 
      included as part of the reference.
    
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2022-03-20 as check problems were not
      corrected in time.

    'length > 1 in coercion to logical' error in check of 'bkmrhat'.

    * Check problems have now been corrected. 

  Found the following (possibly) invalid DOIs:
    DOI: 10.1093/biostatistics/kxu058
      From: DESCRIPTION
      Status: Forbidden
      Message: 403
      
    * This is the DOI provided on the article website <https://academic.oup.com/biostatistics/article/16/3/493/269719>.

## Downstream dependencies
I have also run R CMD check on downstream dependencies with no issues.
