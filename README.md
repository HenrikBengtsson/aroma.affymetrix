# aroma.affymetrix: Analysis of Large Affymetrix Microarray Data Sets

Please see [The Aroma Project](http://www.aroma-project.org/) for more
details on how to use the package and for joining the user forum where
you can ask questions.


## Installation

R package aroma.affymetrix is available on
[CRAN](http://cran.r-project.org/package=aroma.affymetrix).  The
easiest to install the package and all of its dependencies (of which
some are on Bioconductor), use
```r
source('http://callr.org/install#aroma.affymetrix')
```


## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/aroma.affymetrix/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/aroma.affymetrix) branch contains the code of the latest release, which is exactly what is currently on [CRAN](https://cran.r-project.org/package=aroma.affymetrix).

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [aroma.affymetrix repository](https://github.com/HenrikBengtsson/aroma.affymetrix).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/aroma.affymetrix">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-affymetrix">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | CRAN        | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://cran.r-project.org/web/checks/check_results_aroma.affymetrix.html"><img border="0" src="http://www.r-pkg.org/badges/version/aroma.affymetrix" alt="CRAN version"></a> | <a href="https://travis-ci.org/HenrikBengtsson/aroma.affymetrix"><img src="https://travis-ci.org/HenrikBengtsson/aroma.affymetrix.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-affymetrix"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/aroma.affymetrix?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/aroma.affymetrix"><img src="https://codecov.io/gh/HenrikBengtsson/aroma.affymetrix/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
