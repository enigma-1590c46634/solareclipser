
<!-- README.md is generated from README.Rmd. Please edit that file -->

# solareclipser

<!-- badges: start -->
<!-- badges: end -->

## Description

`solareclipser` is an R package to interface SOLAR and to run its main
models: polygenic, association and linkage.

| Model                                                                                                 | SOLAR cmd                                                                                           | `solareclipser` function | Input data format                             |
|-------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|--------------------------|-----------------------------------------------|
| [Polygenic](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#polygenic-model-in-solar)     | [polygenic](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#polygenic)   | solarPolygenic           | phen: R data.frame (ID/FA/MO pedigree-fields) |
| [Linkage](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#linkage-model-in-solar)         | [multipoint](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#multipoint) | solarMultipoint          | MIBD: SOLAR                                   |
| [Association](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#association-model-in-solar) | [mga](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#mga)               | solarAssoc               | SNP: R matrix, SOLAR, PLINK (markers/imputed) |

The `solareclipser` package also enables parallel computation of
[association](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#parallel-computation)
/
[linkage](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#parallel-computation-1)
models to make genome-wide scans more efficient.

## Prerequisites 

The package requires a version of [SOLAR-Eclipse](https://www.nitrc.org/projects/se_linux).

## Installation

You can install the development version of `solareclipser` like so:

``` r
library(devtools)
install_github("enigma-1590c46634/solareclipser")
```

## Minimal example

The following code shows how the main functions (`solarPolygenic`,
`solarMultipoint` and `solarAssoc`) of the package work.

``` r
# load library
library(solareclipser)

# load data set
data(dat30)

# univariate polygenic model
mod1 <- solarPolygenic(trait1 ~ 1, dat30)

# bivariate polygenic model
mod2 <- solarPolygenic(trait1 + trait2 ~ 1, dat30,
  polygenic.options = '-testrhoe -testrhog')

# specify directory with IBD matrices and run linkage model
mibddir <- system.file('extdata', 'solarOutput',
  'solarMibdsCsv', package = 'solareclipser')
link <- solarMultipoint(trait1 ~ 1, dat30,
  mibddir = mibddir, chr = 5)

# run association model in parallel
assoc <- solarAssoc(trait1 ~ 1, dat30, cores = 2,
  snpcovdata = genocovdat30, snpmap = mapdat30)
#> [1] "snp.geno-list"
```

## Tutorial

See [tutorial](inst/doc/tutorial.md) for more information.

## SOLAR references

- The new [SOLAR web page](https://solar-eclipse-genetics.org/)
  (SOLAR-Eclipse)
- [Appendix 1. SOLAR Command
  Descriptions](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html)
