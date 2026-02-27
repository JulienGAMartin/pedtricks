# A function to simulate QTL and/or SNP data.

Simulates a chromosome of arbitrary length with arbitrary numbers,
types, and spacings of genetic loci over arbitrary pedigrees.

## Usage

``` r
genome_sim(
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
```

## Arguments

- pedigree:

  A pedigree

- founders:

  A vector of indicator variables denoting founder status (1=founder,
  0=non-founder)

- positions:

  Genome locations in cM for markers

- initHe:

  Initial levels of expected heterozygosity

- mutationType:

  A vector of locus types - see details

- mutationRate:

  A vector of mutation rates

- phenotyped:

  A vector of IDs of those individuals for which to return phenotypic
  data

- founderHaplotypes:

  A matrix or dataframe containing founder haplotypes

- genotyped:

  A vector of IDs of those individuals for which to return genotypic
  data

- returnG:

  If 'y' then genotypic data for all loci (including `cIAM` loci) will
  be returned.

- initFreqs:

  A list of allele frequencies for all loci. If `initFreqs` is
  specified, it will override information from `initHe`. `extractA` from
  package `MasterBayes` can be used to obtain obtain `initFreqs` form a
  sample of genotypes. For `cIAM` loci, allele names in `initFreqs`
  should be allelic substitution effects.

## Value

- Phenotypes:

  A vector of phenotypes. Calculated as the sum of all allelic effects.
  Scaling is currently left to be done post-hoc.

- MarkerData:

  A vector of marker genotypes, i.e. alleles at all loci except those
  designated \`cIAM'

## Details

Valid mutation types are `Micro', `Dom', `dIAM' and `cIAM', for
microsatellite, dominant (AFLP), discrete infinite alleles mutation
model loci (SNPs), and continuous infinite alleles mutation model loci
(polymorphisms effecting phenotypic variation). cIAM loci have
mutational allelic substitution effects taken drawn from a normal
distribution with mean 0 and variance 1.

## See also

[`phen_sim`](https://juliengamartin.github.io/pedtricks/reference/phen_sim.md),
[`micro_sim`](https://juliengamartin.github.io/pedtricks/reference/micro_sim.md)

## Examples

``` r
testData <- as.data.frame(matrix(
  c(
    1,      NA,     NA,     1,      1,      1,      2,      2,
    2,      NA,     NA,     1,      1,      1,      2,      2,
    3,      NA,     NA,     1,      1,      1,      2,      2,
    4,      NA,     NA,     1,      0,      1,      2,      2,
    5,      NA,     NA,     1,      0,      1,      2,      2,
    6,      1,      4,      0,      -1,     2,      3,      3,
    7,      1,      4,      0,      -1,     2,      3,      3,
    8,      1,      4,      0,      -1,     2,      3,      3,
    9,      1,      4,      0,      -1,     2,      3,      3,
    10,     2,      5,      0,      -1,     2,      3,      3,
    11,     2,      5,      0,      -1,     2,      3,      3,
    12,     2,      5,      0,      -1,     2,      3,      3,
    13,     2,      5,      0,      -1,     2,      3,      3,
    14,     3,      5,      0,      -1,     2,      3,      3,
    15,     3,      5,      0,      -1,     2,      3,      3,
    16,     3,      5,      0,      -1,     2,      3,      3,
    17,     3,      5,      0,      -1,     2,      3,      3
  ),
  17, 8,
  byrow = TRUE
))

names(testData) <- c(
  "id", "dam", "sire", "founder", "sex",
  "cohort", "first", "last"
)
pedigree <- as.data.frame(cbind(
  testData$id, testData$dam,
  testData$sire
))
for (x in 1:3) pedigree[, x] <- as.factor(pedigree[, x])
names(pedigree) <- c("id", "dam", "sire")
pedigree
#>    id  dam sire
#> 1   1 <NA> <NA>
#> 2   2 <NA> <NA>
#> 3   3 <NA> <NA>
#> 4   4 <NA> <NA>
#> 5   5 <NA> <NA>
#> 6   6    1    4
#> 7   7    1    4
#> 8   8    1    4
#> 9   9    1    4
#> 10 10    2    5
#> 11 11    2    5
#> 12 12    2    5
#> 13 13    2    5
#> 14 14    3    5
#> 15 15    3    5
#> 16 16    3    5
#> 17 17    3    5

## make up some microsatellite and gene allele frquencies:
sampleGenotypes <- as.data.frame(matrix(c(
  1, 2, -1.32, 0.21, 2, 1, 0.21, 0.21
), 2, 4, byrow = TRUE))
testFreqs <- extractA(sampleGenotypes)

## note that alleles at the gene locus are given as their
## allelic substitution effects:
testFreqs
#> $V2
#>   2   1 
#> 0.5 0.5 
#> 
#> $V4
#>  0.21 -1.32 
#>  0.75  0.25 
#> 

## simulate data for these indivdiuals based on a single QTL
## with two equally alleles with balanced frequencies in the
## founders, linked (2 cM) to a highly  polymorphic microsatellite:
genome_sim(
  pedigree = pedigree, founders = testData$founder, positions = c(0, 2),
  mutationType = c("Micro", "cIAM"), mutationRate = c(0, 0),
  initFreqs = testFreqs, returnG = "y"
)
#> Processing pedigree... 
#> 0%                     50%                     100% 
#> |                       |                       | 
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> -
#> 
#>  ...done. 
#> 
#>  Calculating phenotypes...
#> done. 
#> 
#>  Tabulating marker genotypes...
#> done. 
#> $Phenotypes
#>    V1    V2
#> 1   1  0.42
#> 2   2 -1.11
#> 3   3  0.42
#> 4   4  0.42
#> 5   5  0.42
#> 6   6  0.42
#> 7   7  0.42
#> 8   8  0.42
#> 9   9  0.42
#> 10 10 -1.11
#> 11 11 -1.11
#> 12 12 -1.11
#> 13 13  0.42
#> 14 14  0.42
#> 15 15  0.42
#> 16 16  0.42
#> 17 17  0.42
#> 
#> $MarkerData
#> [1] doubleids markers  
#> <0 rows> (or 0-length row.names)
#> 
#> $genomes
#>    doubleids V2    V3
#> 1          1  2  0.21
#> 2          1  2  0.21
#> 3          2  2  0.21
#> 4          2  2 -1.32
#> 5          3  2  0.21
#> 6          3  1  0.21
#> 7          4  1  0.21
#> 8          4  2  0.21
#> 9          5  1  0.21
#> 10         5  2  0.21
#> 11         6  2  0.21
#> 12         6  1  0.21
#> 13         7  2  0.21
#> 14         7  1  0.21
#> 15         8  2  0.21
#> 16         8  1  0.21
#> 17         9  2  0.21
#> 18         9  2  0.21
#> 19        10  2 -1.32
#> 20        10  2  0.21
#> 21        11  2 -1.32
#> 22        11  1  0.21
#> 23        12  2 -1.32
#> 24        12  2  0.21
#> 25        13  2  0.21
#> 26        13  2  0.21
#> 27        14  1  0.21
#> 28        14  2  0.21
#> 29        15  2  0.21
#> 30        15  2  0.21
#> 31        16  2  0.21
#> 32        16  1  0.21
#> 33        17  2  0.21
#> 34        17  1  0.21
#> 
## since we specified returnG='y', we can check that
## the phenotypes add up to the
## allelic substitution effects for the second locus.
```
