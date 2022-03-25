## Resubmission

* If there are references describing the methods in your package, please
  add these in the description field of your DESCRIPTION file in the form
  authors (year) <doi:...>
  authors (year) <arXiv:...>
  authors (year, ISBN:...)
  or if those are not available: <https:...>
  with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
  auto-linking.
  (If you want to add a title as well please put it in quotes: "Title")
  
The paper reference has been added to the description field.

* Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar)
  
A \value field has been added to the documentation for exported methods. 

* Please add small executable examples in your Rd-files to illustrate the
  use of the exported function but also enable automatic testing.
  
Executable examples are now included in Rd-files for exported functions.

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
    
    * Author names were incuded as part of the methods reference, as requested by CRAN.
    
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
      
    * This was the DOI provided on the article website <https://academic.oup.com/biostatistics/article/16/3/493/269719>.

## Downstream dependencies
I have also run R CMD check on downstream dependencies with no issues.
