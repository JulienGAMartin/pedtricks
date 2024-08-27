# TO DO list

- [ ] document and fix the plot function for ped_stats
- [ ] better documentation overall
- [ ] clean all linter warnings
- [ ] convert plots from plot.ped_stats to ggplot format
- [X] add ORCID for all of us
- [X] add functions to plot relatedness matrix


## ggpedigree.R
- [ ] replace draw_ped by ggpedigree - make sure functionality matches
- [X] optimize ggpedigree x-coordinates (similar to PedView programme)
- [X] add sex information (took this out temporarily)
- [X] add data checks e.g. are vectors the same length etc.

## nadiv integration

since we are already dependent on nadiv why not remove the redundant functions

- internal
  - `prune` -> `prunePed`
  - `orderPed` -> `prepPed`
- exported
  - `convert_ped` -> `numPed`
  - `fix_ped` -> `prepPed`