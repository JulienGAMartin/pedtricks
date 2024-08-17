# pedtricks package

This package is taking the core functions from pedantics, removed the dependencies on archived packages ("MasterBayes") and removed function based on C libraries that are outdated, providing an easier to maintain package

New functions will include a ggplot function to plot pedigree `ggpedigree()` and more

Install with

```r
remotes::install_github(
  "JulienGAMartin/pedtricks"
)
```

## Changes from pedantics

### Defunct functions
* `rpederr` and `fpederr` have been removed from the package

### major changes
* `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
* `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor

### functions renamed
* `fixPedigree()` is now `fix_ped()`
* `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`
