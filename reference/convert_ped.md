# Converts a pedigree with individuals specified as factors to a numeric pedigree

Some internal `pedtricks` modules require that pedigrees be specified
only by numerical values, or including numerical values for missing
data. This function provides the conversion to numeric but also back to
factors if needed

## Usage

``` r
convert_ped(type = "numeric", id, sire, dam, missingVal = NA, key = NULL)
```

## Arguments

- type:

  define how to convert the pedigree so "numeric" or "factor"

- id:

  Individual identifiers - pass using
  [`as.character()`](https://rdrr.io/r/base/character.html)

- sire:

  Sire codes - pass `using as.character()`

- dam:

  Dam codes - pass `using as.character()`

- missingVal:

  the indicator that should be substituted for missing values

- key:

  A dataframe, as produced by `convert_ped`, specifying factor codes for
  numeric values in id, sire, and dam

## Value

- numericPedigree:

  The factor pedigree in numeric form

- idKey:

  A key to facilitate conversion back to the original identifiers

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

## make the test pedigree numeric with NAs denoted by -1
convert_ped(
  type = "numeric",
  id = as.character(pedigree[, 1]),
  dam = as.character(pedigree[, 2]),
  sire = as.character(pedigree[, 3]),
  missingVal = -1
)
#> $numericPedigree
#>    id sire dam
#> 1   3   -1  -1
#> 2   4   -1  -1
#> 3   5   -1  -1
#> 4   1   -1  -1
#> 5   2   -1  -1
#> 6  14    1   3
#> 7  15    1   3
#> 8  16    1   3
#> 9  17    1   3
#> 10  6    2   4
#> 11  7    2   4
#> 12  8    2   4
#> 13  9    2   4
#> 14 10    2   5
#> 15 11    2   5
#> 16 12    2   5
#> 17 13    2   5
#> 
#> $idKey
#>    pn   pf
#> 1   3   m1
#> 2   4   m2
#> 3   5   m3
#> 4   1   d4
#> 5   2   d5
#> 6  14   o6
#> 7  15   o7
#> 8  16   o8
#> 9  17   o9
#> 10  6  o10
#> 11  7  o11
#> 12  8  o12
#> 13  9  o13
#> 14 10  o14
#> 15 11  o15
#> 16 12  o16
#> 17 13  o17
#> 18 -1 <NA>
#> 
```
