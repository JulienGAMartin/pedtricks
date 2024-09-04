# pedtricks 0.4

## 0.4.2

- more fix for CRAN submission

## 0.4.1

-`ped_stats`: more optimisation to run even faster.  Now running in 0.7s on the gryphons pedigree instead of 160s so 210x faster. Loops are slow


## 0.4.0

- `ped_stats`: optimized to run much faster (I think it can even be improved).  Now running in 1.6s on the gryphons pedigree instead of 160s so 100x faster. Now using sparse matrix format for A and removing some sluggish loops in estimations.

# pedtricks 0.3

## 0.3.2

- `ggpedigree()`
  - adding functionality to highlight phenotyped individuals with `pheno` argument

## 0.3.1

- `ggpedigree()`
  - draw pedigree using ggplot2
  - compatible with dplyr/tidyverse use
  - return a ggplot object that can be modified with ggplot functions
  - optimized graph layout with igrah
- `draw_pedA()`
  - returns a trellis plot object for visualizing relatedness among pedigreed individuals

## 0.3.0

- `ggpedigree()`
  - draw pedigree using ggplot2
  - compatible with dplyr/tidyverse use
  - return a ggplot object that can be modified with ggplot functions


# pedtricks 0.1

- Initial transfer and update of CRAN archived [`pedantics`](https://CRAN.R-project.org/package=pedantics) package under the new name `pedtricks`.

### DEPRECATED from `pedantics`

- `rpederr()`and `fpederr()` are  not implemented
- `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
- `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor
      
### Small changes

- `fixPedigree()` is now `fix_ped()`
- `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`

