# Manipulating pedigrees to prepare them for requirements of subsequent analyses Prepares a pedigree to conform with requirements of many softwares used in quantitative genetic analysis, as well as for many of the functions in pedtricks.

Manipulating pedigrees to prepare them for requirements of subsequent
analyses Prepares a pedigree to conform with requirements of many
softwares used in quantitative genetic analysis, as well as for many of
the functions in pedtricks.

## Usage

``` r
fix_ped(ped, dat = NULL)
```

## Arguments

- ped:

  An ordered pedigree with 3 columns: id, dam, sire

- dat:

  An optional data frame, the same length as the pedigree

## Value

Returns a pedigree in which all individuals that exist in the dam and
sire columns are represented by their own record lines, occurring before
the records of their first offspring. If data are supplied, then fix_ped
will return a dataframe, the first three columns are the 'fixed'
pedigree, and the following columns of which contain appropriately
reordered data.

## Examples

``` r
##  a valid pedigree, i.e., no loops, no bisexuality, etc.,
## but where not all parents have a record line, and where
## parents do not necessarily occur before their offspring:
pedigree <- as.data.frame(matrix(c(
  10, 1, 2,
  11, 1, 2,
  12, 1, 3,
  13, 1, 3,
  14, 4, 5,
  15, 6, 7,
  4, NA, NA,
  5, NA, NA,
  6, NA, NA,
  7, NA, NA
), 10, 3, byrow = TRUE))
names(pedigree) <- c("id", "dam", "sire")
pedigree
#>    id dam sire
#> 1  10   1    2
#> 2  11   1    2
#> 3  12   1    3
#> 4  13   1    3
#> 5  14   4    5
#> 6  15   6    7
#> 7   4  NA   NA
#> 8   5  NA   NA
#> 9   6  NA   NA
#> 10  7  NA   NA
fixed_pedigree <- fix_ped(ped = pedigree)
fixed_pedigree
#>    id  dam sire
#> 7   4 <NA> <NA>
#> 8   5 <NA> <NA>
#> 9   6 <NA> <NA>
#> 10  7 <NA> <NA>
#> 31  1 <NA> <NA>
#> 61  2 <NA> <NA>
#> 63  3 <NA> <NA>
#> 1  10    1    2
#> 2  11    1    2
#> 3  12    1    3
#> 4  13    1    3
#> 5  14    4    5
#> 6  15    6    7
```
