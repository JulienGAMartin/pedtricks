# A function to simulated phenotypic data

Simulates phenotypic data across arbitrary pedigrees. ` phen_sim`
simulate direct, maternal and paternal genetic and environmental effects
for an arbitrary number of traits with arbitrary patterns of missing
data.

## Usage

``` r
phen_sim(
  pedigree,
  traits = 1,
  randomA = NULL,
  randomE = NULL,
  parentalA = NULL,
  parentalE = NULL,
  sampled = NULL,
  records = NULL,
  returnAllEffects = FALSE,
  verbose = TRUE
)
```

## Arguments

- pedigree:

  A pedigree

- traits:

  The number of traits for which data should be simulated.

- randomA:

  An additive genetic covariance matrix, with dimensions a multiple of
  traits - see details

- randomE:

  An additive environmental covariance matrix, with dimensions a
  multiple of traits - see details

- parentalA:

  A vector indicating which effects in `randomA` (if any) to treat as
  parental effects

- parentalE:

  A vector indicating which effects in `randomE` (if any) to treat as
  parental effects

- sampled:

  A vector indicating which individuals are sampled

- records:

  A single value, array of matrix specifying data record availability -
  see details

- returnAllEffects:

  If `TRUE` then all individual breeding values and environmental
  effects are returned

- verbose:

  If `TRUE` provide a progress bar and messages, Default: TRUE

## Value

- phenotypes:

  A dataframe containing phenotypes for all individuals specified to
  have records.

- allEffects:

  (optional) A dataframe with all direct and indirect genetic and
  environmental effects.

## Details

`randomA` and `randomE` are square matrices with dimension equal to the
sum of the number direct and indirect effects. This must be a multiple
of the number of traits, i.e. if an indirect effect is to be simulated
for only one of multiple traits, those traits with no indirect effect
should be included with (co)variances of zero.

`parentalA` and `parentalE` are optional vectors of characters
indicating which trait positions in `randomA` and `randomE` are to be
treated as indirect effects, and which effects to treat as maternal or
paternal. Valid values are 'd', 'm', and 'p', for direct, maternal
indirect and paternal indirect effects, respectively.

`records` can be specified either (1) as a single value to be applied to
all individuals and traits, (2) as a vector the same length as the
number of traits, representing trait-specific rates to be applied
uniformly to all individuals, or (3) as data frames with rows for each
individual and columns for each trait. In the third option, observed
patterns of data availability can be simulated by supplying 0s and 1s
for missing and available individual genotypes, respectively.

## See also

[`micro_sim`](https://juliengamartin.github.io/pedtricks/reference/micro_sim.md),
[`genome_sim`](https://juliengamartin.github.io/pedtricks/reference/genome_sim.md)

## Examples

``` r
## make up a pedigree
id <- c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9")
dam <- c(NA, NA, NA, "a1", "a1", "a1", "a4", "a4", "a4")
sire <- c(NA, NA, NA, "a2", "a2", "a2", "a5", "a6", "a6")
pedigree <- as.data.frame(cbind(id, sire, dam))

traits <- 2
## no correlations
randomA <- diag(4)
randomE <- diag(4)
parentalA <- c("d", "d", "m", "m")
parentalE <- c("d", "d", "m", "m")

## generate phenoypic data based on this architecture
phen_sim(
  pedigree = pedigree, traits = 2, randomA = randomA, randomE = randomE,
  parentalA = parentalA, parentalE = parentalE
)
#> Simulating breeding values and environmental effects...
#> done. 
#> Calculating phenotypes...
#> done. 
#> $phenotypes
#>    id     trait_1    trait_2
#> a1 a1 -0.86946859  0.9596854
#> a2 a2 -0.48857480 -0.1537290
#> a3 a3 -3.61599848 -0.5482039
#> a4 a4 -0.03733646  2.8267077
#> a5 a5  2.44432122  0.2253911
#> a6 a6 -0.19259122  4.1744163
#> a7 a7  2.93255113 -0.1206532
#> a8 a8  2.31069940  1.5540842
#> a9 a9  1.59106545 -0.3333576
#> 

## let's do it again but see how the phenotypes were composed
phen_sim(
  pedigree = pedigree, traits = 2, randomA = randomA, randomE = randomE,
  parentalA = parentalA, parentalE = parentalE, returnAllEffects = TRUE
)
#> Simulating breeding values and environmental effects...
#> done. 
#> Calculating phenotypes...
#> done. 
#> $phenotypes
#>    id    trait_1    trait_2
#> a1 a1 -3.0074621 -2.3891237
#> a2 a2  1.9386390 -3.6175357
#> a3 a3  0.2506263 -0.2853697
#> a4 a4  3.9857360  0.6870305
#> a5 a5  1.1865358 -0.2253990
#> a6 a6  1.7904845 -2.1249351
#> a7 a7  1.0614037 -2.8693030
#> a8 a8 -0.9673147 -3.3481675
#> a9 a9 -2.5694630 -1.0926158
#> 
#> $allEffects
#>    id sire  dam       a_tr1      a_tr2   bv_m_tr1   bv_m_tr2       e_tr1
#> a1 a1 <NA> <NA> -1.57457485 -0.1963194  0.9150601  1.0786365 -0.25537509
#> a2 a2 <NA> <NA>  1.45280285 -1.6532811 -2.7618675  0.2643692  1.13866752
#> a3 a3 <NA> <NA>  1.02457159 -0.6591855 -0.3863044 -0.2424989 -2.13709831
#> a4 a4   a2   a1  0.09509698 -1.9720891 -0.4109696 -1.0096539  1.76313061
#> a5 a5   a2   a1 -0.10713729 -1.2569317 -1.1217587  1.0779737 -0.83383531
#> a6 a6   a2   a1 -0.20632149 -0.7587147 -0.8451777  0.8653931 -0.13070247
#> a7 a7   a5   a4  1.10761619 -1.2130379 -0.8085164 -0.0716146 -0.06286163
#> a8 a8   a6   a4  0.13259750 -0.9917643 -1.0377629  0.5941052 -1.11656132
#> a9 a9   a6   a4 -0.55508516 -2.2672473 -0.1095806 -0.3692617 -2.03102701
#>          e_tr2     e_m_tr1     e_m_tr2 P_bv_m_tr1   P_bv_m_tr2   P_e_m_tr1
#> a1 -0.36923822  1.21244833 -0.44536301  0.5280695 -0.968152923 -1.70558168
#> a2  0.19052348 -1.36639361 -1.29810260 -0.4802665 -0.918715157 -0.17256490
#> a3 -0.09620553  0.69450905 -0.47459123  1.3305972  0.008731295  0.03255579
#> a4  2.02584606  0.42761881 -0.56647935  0.9150601  1.078636538  1.21244833
#> a5  0.39825915 -1.47453793 -0.39544120  0.9150601  1.078636538  1.21244833
#> a6 -1.99949395 -0.99593086  0.42011400  0.9150601  1.078636538  1.21244833
#> a7 -0.08013185 -0.03228341 -0.71898093 -0.4109696 -1.009653917  0.42761881
#> a8 -0.78026990 -1.77695853 -0.42783487 -0.4109696 -1.009653917  0.42761881
#> a9  2.75076475  1.51366972  0.03398894 -0.4109696 -1.009653917  0.42761881
#>     P_e_m_tr2   Phen_tr1   Phen_tr2
#> a1 -0.8554131 -3.0074621 -2.3891237
#> a2 -1.2360629  1.9386390 -3.6175357
#> a3  0.4612901  0.2506263 -0.2853697
#> a4 -0.4453630  3.9857360  0.6870305
#> a5 -0.4453630  1.1865358 -0.2253990
#> a6 -0.4453630  1.7904845 -2.1249351
#> a7 -0.5664793  1.0614037 -2.8693030
#> a8 -0.5664793 -0.9673147 -3.3481675
#> a9 -0.5664793 -2.5694630 -1.0926158
#> 
```
