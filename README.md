The R package `bkmr` implements Bayesian kernel machine regression, a statistical approach for estimating the joint health effects of multiple concurrent exposures. Additional information on the statistical methodology and on the computational details are provided in [Bobb et al. 2015](https://academic.oup.com/biostatistics/article/16/3/493/269719). More recent extensions, details on the software, and worked-through examples are provided in [Bobb et al. 2018](https://ehjournal.biomedcentral.com/articles/10.1186/s12940-018-0413-y).

You can install the latest released version of `bkmr` from CRAN with:
```R
install.packages("bkmr")
```
Or the latest development version from github with:
```R
install.packages("devtools")
devtools::install_github("jenfb/bkmr")
```

For a general overview and guided examples, go to https://jenfb.github.io/bkmr/overview.html.

For examples from the software paper, please see

* https://jenfb.github.io/bkmr/SimData1
* https://jenfb.github.io/bkmr/ProbitEx