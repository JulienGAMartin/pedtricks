# pedtricks package

<!-- badges: start -->
[![R-CMD-check](https://github.com/juliengamartin/pedtricks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/juliengamartin/pedtricks/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/pedtricks)](https://cran.r-project.org/package=pedtricks)
<!-- badges: end -->

<!--[![Codecov test coverage](https://codecov.io/gh/r-lib/devtools/branch/main/graph/badge.svg)](https://app.codecov.io/gh/r-lib/devtools?branch=main)-->

Functions for sensitivity and power analysis,
  for calculating statistics describing pedigrees from wild populations,
  and for visualizing pedigrees.
  
  This package is essentially a reboot of the now archived `pedantics` package. We took the `pedantics` package, removed the dependencies on archived packages (`MasterBayes`), renamed some functions and removed 2 (`fpederr()` and `rpederr()`). 
  
## Installation

Install with

```r
remotes::install_github(
  "JulienGAMartin/pedtricks"
)
```

## Changes from pedantics

see the [News page of pedtricks website](https://juliengamartin.github.io/pedtricks/news/index.html) for the list of all the changes

### Defunct functions
* `rpederr()` and `fpederr()` have been removed from the package.

The questions about rates of pedigree error typically encountered, and their effects on inference of elementary QG parameters are basically settled (Charmantier & Réale 2005, Morrissey *et al.* 2007, Firth *et al.* 2015, Bourret & Garant 2017), and outstanding questions will not plausibly be supported by these functions or any straightforward modifications thereof.

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

## References

Bourret, A., and D. Garant. 2017. An assessment of the reliability of quantitative genetics estimates in study systems with high rate of extra pair paternity and low recruitment.  Heredity 118: 229-38.

Charmantier, A., and D. Réale. 2005. How do misassigned paternities affect estimation of heritability in the wild? Molecular Ecology 14: 2839-50

Firth, J.A. J.D. Hadfield, A.W. Santure, J. Slate and B.C. Sheldon. 2015. The influence of non-random extra pair paternity on heritability estimates derived from wild pedigrees. Evolution 69: 1336-44.

Morrissey, M.B., A.J. Wilson, J.M. Pemberton, and M.M. Ferguson. 2007. A framework for power and sensitivity analysis for quantitative genetic studies of natural populations, and a case study in Soay sheep (Ovis aries). Journal of Evolutionary Biology 20: 2309-21.


