# Simulates microsatellite data across a pedigree.

Uses a pedigree with parents identified for all non-founding individuals
and simulates microsatellite genotypes

## Usage

``` r
micro_sim(
  pedigree,
  genFreqs = NULL,
  genotypesSample = NULL,
  knownGenotypes = NULL,
  records = NULL,
  eRate1 = 0,
  eRate2 = 0,
  eRate3 = 0
)
```

## Arguments

- pedigree:

  A pedigree

- genFreqs:

  (optional) A list of allele frequencies, can be produced with
  `extractA`

- genotypesSample:

  (required if `genFreqs` is not supplied) a sample of genotypes from
  which to estimate population allele frequencies

- knownGenotypes:

  (not yet implemented) a data frame of genotypes for (potentially a
  subset) of founder individuals

- records:

  Record availability, see details.

- eRate1:

  The rate of genotypic substitution errors, i.e., when a true genotype
  at a given locus is replaced by a pair of alleles selected at random
  based on the population allele frequencies

- eRate2:

  The rate of allelic substitution errors, i.e. when an allele is
  erroneously replaced at a given locus by an allele chosen at random
  based on the population allele frequencies

- eRate3:

  The rate of large allele dropouts, simulated by setting the value of
  the larger allele at a locus to the value of the smaller allele

## Value

- trueGenotypes:

  A data frame of true genotypes

- observedGenotypes:

  A data frame of plausible observed genotypes, given specified patterns
  of missingness and errors.

## Details

Error rates and data availability rates can be specified as either (1)
single values to be applied to all individuals and all loci, (2) as a
vector the same length as the number of loci, representing
locus-specific rates to be applied uniformly to all individuals, or (3)
as data frames with rows for each individual and columns for each locus.
In the third option, observed patterns of data availability can be
simulated by supplying 0s and 1s for missing and available individual
genotypes, respectively.

## See also

[`phen_sim`](https://juliengamartin.github.io/pedtricks/reference/phen_sim.md),
[`genome_sim`](https://juliengamartin.github.io/pedtricks/reference/genome_sim.md)

## Examples

``` r
pedigree <- as.data.frame(matrix(c(
  "m1",   NA,     NA,
  "m2",   NA,     NA,
  "m3",   NA,     NA,
  "d4",   NA,     NA,
  "d5",   NA,     NA,
  "o6",   "m1",   "d4",
  "o7",   "m1",   "d4",
  "o8",   "m1",   "d4",
  "o9",   "m1",   "d4",
  "o10",  "m2",   "d5",
  "o11",  "m2",   "d5",
  "o12",  "m2",   "d5",
  "o13",  "m2",   "d5",
  "o14",  "m3",   "d5",
  "o15",  "m3",   "d5",
  "o16",  "m3",   "d5",
  "o17",  "m3",   "d5"
), 17, 3, byrow = TRUE))
names(pedigree) <- c("id", "dam", "sire")
for (x in 1:3) pedigree[, x] <- as.factor(pedigree[, x])

## some sample genotypes, very simple, two markers with He = 0.5
sampleGenotypes <- as.data.frame(matrix(c(
  1, 2, 1, 2, 2, 1, 2, 1
), 2, 4, byrow = TRUE))
## locus names
names(sampleGenotypes) <- c("loc1a", "loc1b", "loc2a", "loc2b")

## simulate some genotypes
micro_sim(pedigree = pedigree, genotypesSample = sampleGenotypes)
#> $trueGenotypes
#>     loc1b loc1b_b loc2b loc2b_b
#> m1      2       2     1       2
#> m2      1       1     1       1
#> m3      1       1     2       1
#> d4      2       2     2       2
#> d5      1       2     2       2
#> o6      2       2     2       2
#> o7      2       2     2       2
#> o8      2       2     2       2
#> o9      2       2     1       2
#> o10     1       2     1       2
#> o11     1       1     1       2
#> o12     1       2     1       2
#> o13     1       2     1       2
#> o14     1       1     1       2
#> o15     1       2     2       2
#> o16     1       1     2       2
#> o17     1       2     2       2
#> 
#> $observedGenotypes
#>     loc1b loc1b_b loc2b loc2b_b
#> m1      2       2     1       2
#> m2      1       1     1       1
#> m3      1       1     2       1
#> d4      2       2     2       2
#> d5      1       2     2       2
#> o6      2       2     2       2
#> o7      2       2     2       2
#> o8      2       2     2       2
#> o9      2       2     1       2
#> o10     1       2     1       2
#> o11     1       1     1       2
#> o12     1       2     1       2
#> o13     1       2     1       2
#> o14     1       1     1       2
#> o15     1       2     2       2
#> o16     1       1     2       2
#> o17     1       2     2       2
#> 
```
