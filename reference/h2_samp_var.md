# Calculates expected sampling variance of a given \$h^2\$ from the A matrix of a pedigree, based on the method of Visscher and Goddard 2015 Genetics

Calculates expected sampling variance of a given \$h^2\$ from the A
matrix of a pedigree, based on the method of Visscher and Goddard 2015
Genetics

## Usage

``` r
h2_samp_var(Ped, h2 = 0.3, plot = FALSE)
```

## Arguments

- Ped:

  A pedigree

- h2:

  vector of heritabilities for which the sampling variance is
  calculated.

- plot:

  logical. If TRUE, a plot of the sampling variance as a function of h2
  is produced (only if length(h2)\>1).

## Value

Returns a vector of sampling variances for each given h2.

## Examples

``` r
# \donttest{
data(gryphons)
pedigree <- gryphons[, 1:3]

h2_samp_var(pedigree,h2=0.3)
#>         0.3 
#> 0.001321648 
# }
```
