# pedtricks 0.2.0

## New

  - `ggpedigree()`
    - compatible with dplyr/tidyverse use
    - return a ggplot object that can be modified with ggplot functions

  - `draw_pedA()`
    - returns a trellis plot object for visualizing relatedness among pedigreed individuals
     
# pedtricks 0.1.0

## NEW
  - Initial transfer and update of CRAN archived [`pedantics`](https://CRAN.R-project.org/package=pedantics) package under the new name `pedtricks`.

## DEPRECATED from `pedantics`

  - `rpederr()`and `fpederr()` are  not implemented
  - `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
  - `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor
      
## Small changes

  - `fixPedigree()` is now `fix_ped()`
  - `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`

