# pedtricks: Visualize, Summarize and Simulate Data from Pedigrees

Sensitivity and power analysis, for calculating statistics describing
pedigrees from wild populations, and for visualizing pedigrees. This is
a reboot of the methods developed by Morrissey and Wilson (2010) [doi:
10.1111/j.1755-0998.2009.02817.x](https://doi.org/%2010.1111/j.1755-0998.2009.02817.x)
.

## Details

- [`pedigreeStats()`](https://juliengamartin.github.io/pedtricks/reference/ped_stats.md)
  and
  [`pedStatSummary()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md)
  have been replaced by
  [`ped_stats()`](https://juliengamartin.github.io/pedtricks/reference/ped_stats.md)
  with a [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods to
  simplify the workflow and allow to get the plots without running the
  statistics each time

- `makePedigreeFactor()` and
  [`makePedigreeNumeric()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md)
  have been combined in
  [`convert_ped()`](https://juliengamartin.github.io/pedtricks/reference/convert_ped.md)
  which convert a pedigree to numeric or factor

- [`fixPedigree()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md)
  is now
  [`fix_ped()`](https://juliengamartin.github.io/pedtricks/reference/fix_ped.md)

- [`genomesim()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md),
  [`microsim()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md)
  and
  [`phensim()`](https://juliengamartin.github.io/pedtricks/reference/pedantics-deprecated.md)
  have been renamed as
  [`genome_sim()`](https://juliengamartin.github.io/pedtricks/reference/genome_sim.md),
  [`micro_sim()`](https://juliengamartin.github.io/pedtricks/reference/micro_sim.md)
  and
  [`phen_sim()`](https://juliengamartin.github.io/pedtricks/reference/phen_sim.md)

For backward compatibility with code using `pedantics` older names are
still usable but not recommended.

## See also

Useful links:

- <https://juliengamartin.github.io/pedtricks/>

- Report bugs at <https://github.com/JulienGAMartin/pedtricks/issues>

## Author

**Maintainer**: Julien Martin <julien.martin@uottawa.ca>
([ORCID](https://orcid.org/0000-0001-7726-6809))

Authors:

- Matthew Wolak <mew0099@auburn.edu>
  ([ORCID](https://orcid.org/0000-0002-7962-0071))

- Susan Johnston <Susan.Johnston@ed.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-5623-8902))

- Joel Pick <joel.l.pick@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6295-3742))

- Michael Morrissey <michael.morrissey@st-andrews.ac.uk>
  ([ORCID](https://orcid.org/0000-0001-6209-0177))
