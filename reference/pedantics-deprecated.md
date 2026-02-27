# Deprecated functions from package pedantics.

The functions listed below are deprecated and are available for backward
compatibility only. They might be defunct in a future version When
possible, alternative functions with similar functionality are also
mentioned.

## Usage

``` r
makemakePedigreeFactor(id, sire, dam, key = NULL)

makePedigreeNumeric(id, sire, dam, missingVal = NA)

drawPedigree(
  Ped,
  cohorts = NULL,
  sex = NULL,
  dat = NULL,
  dots = "n",
  plotfull = "y",
  writeCohortLabels = "n",
  links = "all",
  sexInd = c(0, 1),
  dotSize = 0.001,
  dataDots = "n",
  dataDots.cex = 2,
  cohortLabs.cex = 1,
  retain = "informative",
  focal = NULL,
  sexColours = c("red", "blue"),
  ...
)

fixPedigree(ped, dat = NULL)

genomesim(
  pedigree,
  founders = NULL,
  positions = NULL,
  initHe = NULL,
  mutationType = NULL,
  mutationRate = NULL,
  phenotyped = NULL,
  founderHaplotypes = NULL,
  genotyped = NULL,
  returnG = "n",
  initFreqs = NULL
)

microsim(
  pedigree,
  genFreqs = NULL,
  genotypesSample = NULL,
  knownGenotypes = NULL,
  records = NULL,
  eRate1 = 0,
  eRate2 = 0,
  eRate3 = 0
)

pedigreeStat()

pedStatSummary()

phensim(
  pedigree,
  traits = 1,
  randomA = NULL,
  randomE = NULL,
  parentalA = NULL,
  parentalE = NULL,
  sampled = NULL,
  records = NULL,
  returnAllEffects = FALSE
)
```

## Value

dependant on the deprecated function see pedantics help for details

## `makePedigreeFactor`

the function is now a special case of the `convert_ped` function when
type is factor

## `makePedigreeNumeric`

the function is now a special case of the `convert_ped` function when
type is numeric

## `drawPedigree`

the function has just been renamed with no other changes for the moment,
but will soon be replace by "ggpedigree" function

## `fixPedigree`

the function has just been renamed with no other changes for the moment

## `genomesim`

the function has just been renamed with no other changes for the moment

## `microsim`

the function has just been renamed with no other changes for the moment

## `pedigreeStats`

the function has been simplified but only the functionality are still
available via `ped_stats` and its summary and plot methods

## `pedStatSummary`

the function has been simplified but only the functionality are still
available via `ped_stats` and its summary and plot methods

## `phensim`

the function has just been renamed with no other changes for the moment
