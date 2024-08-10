# pedtricks package

It seems that pedantics is not maintained anymore. This package is taking the core functions from pedantics and removed the dependencies on archived packages ("MasterBayes") and removed function based on C libraries that are outdated.

New functions will include a ggplot function to plot pedigree `ggpedigree()` and more

Install with

```r
remotes::install_github(
  "JulienGAMartin/pedtricks"
)
```

## Changes from pedantics

* `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
* `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor
* `fixPedigree()` is now `fix_ped()`
* `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`

For backward compatibility with code using `pedantics` older names are still usable but not recommended.