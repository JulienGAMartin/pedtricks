# Allele Frequencies

Extracts allele frequencies from genotype data

## Usage

``` r
extractA(G, marker.type = "MSW")
```

## Arguments

- G:

  data frame or list of `genotype` objects

- marker.type:

  `"MSW"` or `"MSC"` for co-dominant markers with Wang's (2004) model of
  genotyping error or CERVUS's model of genotyping error
  (Marshall, 1998) or `"AFLP"` for dominant markers.

## Value

list of allele frequencies at each loci

## Author

Jarrod Hadfield <j.hadfield@ed.ac.uk>

## Examples

``` r
# \donttest{
data(WarblerG)

A <- extractA(WarblerG)
A[[1]]
#>        272        279        284 
#> 0.50822368 0.45230263 0.03947368 
# }
```
