The R package `bkmr` implements Bayesian kernel machine regression, a statistical approach for estimating the joint health effects of multiple concurrent exposures. Additional information on the statistical methodology and on the computational details are provided in [Bobb et al. 2015](http://biostatistics.oxfordjournals.org/content/16/3/493).

You can install bkmr from github with:

```R
install.packages("devtools")
devtools::install_github("jenfb/bkmr")
```

The main function implementing BKMR is `kmbayes`. For more information, type

```R
?kmbayes
```

For a general overview and guided examples, go to [https://jenfb.github.io/bkmr/overview.html].
