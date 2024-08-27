# pedtricks package

<!-- badges: start -->
[![R-CMD-check](https://github.com/juliengamartin/pedtricks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/juliengamartin/pedtricks/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/pedtricks)](https://cran.r-project.org/package=pedtricks)
<!-- badges: end -->

<!--[![Codecov test coverage](https://codecov.io/gh/r-lib/devtools/branch/main/graph/badge.svg)](https://app.codecov.io/gh/r-lib/devtools?branch=main)-->

This package is taking the core functions from pedantics, removed the dependencies on archived packages ("MasterBayes") and removed function based on C libraries that are outdated, providing an easier to maintain package

Install with

```r
remotes::install_github(
  "JulienGAMartin/pedtricks"
)
```

## Changes from pedantics

see https://juliengamartin.github.io/pedtricks/news/index.html for the list of all the changes

### Defunct functions
* `rpederr` and `fpederr` have been removed from the package

### major changes
* `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
* `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor

### functions renamed
* `fixPedigree()` is now `fix_ped()`
* `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`

## future for 1.0.0
* `ggpedigree` a ggplot2 style function to plot pedigrees to replace `drawPedigree`
* reimplement all plot methods from ped_stats using ggplot
* add a plotting function for relatedness matrices
* ...