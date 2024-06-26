[![cran version](http://www.r-pkg.org/badges/version/solareclipser)](https://cran.r-project.org/web/packages/solareclipser)
[![downloads](http://cranlogs.r-pkg.org/badges/solareclipser)](http://cranlogs.r-pkg.org/badges/solareclipser)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/solareclipser)](http://cranlogs.r-pkg.org/badges/grand-total/solareclipser)
[![research software impact](http://depsy.org/api/package/cran/solareclipser/badge.svg)](http://depsy.org/package/r/solareclipser)

This is a fork of: <https://github.com/ugcd/solareclipser>

## Table of content

* [About solareclipser](#about-solareclipser)
  * [solareclipser references](#solareclipser-references)
  * [Install](#install)
  * [Quick start](#quick-start)
  * [Citation](#citation)
  * [Rationale](#rationale)
* [FAQ](faq)
* [SOLAR references](#solar-references)

## About solareclipser

![](docs/figures/solareclipser-models.png)

SOLAR is known as the old-school player in the quantitative trait loci (QTLs) mapping field ([>2600](https://scholar.google.es/citations?view_op=view_citation&hl=en&user=AjEIQ3MAAAAJ&citation_for_view=AjEIQ3MAAAAJ:u5HHmVD_uO8C) citations),
which variance components or linear mixed models are implemented especially for the analysis of related individuals in pedigrees.

`solareclipser` is an R package to interface SOLAR and to run its main models: polygenic, association and linkage. 


| Model |	SOLAR cmd |	`solareclipser` function |	Input data format |
|-------|---------------|---------------------|-----------------------|
| [Polygenic](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#polygenic-model-in-solar) | [polygenic](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#polygenic) | solarPolygenic | phen: R data.frame (ID/FA/MO pedigree-fields) |
| [Linkage](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#linkage-model-in-solar) | [multipoint](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#multipoint) | solarMultipoint | MIBD: SOLAR |
| [Association](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#association-model-in-solar) | [mga](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#mga) | solarAssoc |	SNP: R matrix, SOLAR, PLINK (markers/imputed) |


The `solareclipser` package also enables parallel computation of [association](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#parallel-computation) / [linkage](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#parallel-computation-1)
models to make genome-wide scans more efficient.


### solareclipser references

* The stable release on [CRAN](https://cran.r-project.org/package=solareclipser)
* The article describing `solareclipser` in [Bioinformatics](http://bioinformatics.oxfordjournals.org/content/32/12/1901)
    * The preprint version is available in [biorxiv](http://biorxiv.org/content/early/2015/12/25/035378) (the content of Section 2 is different)
* Vignettes 
  * R code [vignettes/](vignettes/)
  * hmlt output:
     1. [tutorial.html](http://ugcd.github.io/solareclipser/vignettes/tutorial.html)
     2. [minimal.html](http://ugcd.github.io/solareclipser/vignettes/minimal.html)
     3. [modelsGAIT1.html](http://ugcd.github.io/solareclipser/vignettes/modelsGAIT1.html)
* Documentation [http://ugcd.github.io/solareclipser/doc/](http://ugcd.github.io/solareclipser/doc/) (out of date)
* Project web [http://ugcd.github.io/solareclipser/](http://ugcd.github.io/solareclipser/) (out of date)

### Install

To install the official release from [CRAN](https://cran.r-project.org/package=solareclipser):

```
install.packages("solareclipser")
```

To install the latest development version from source on GitHub (master branch): 

```
library(devtools)
install_github("ugcd/solareclipser")
```

_Note 1_: Starting from version 3.*, `solareclipser` is not supported for Windows. 
[DESCRIPTION](https://github.com/ugcd/solareclipser/blob/master/DESCRIPTION) file has a special line `OS_type: unix`.
This is a clear limitation of the `solareclipser` package that comes from the dependency on SOLAR.
See also the SOLAR [FAQ](http://solar-eclipse-genetics.org/faq.html) 
and the question 4 "Can SOLAR-Eclipse only be downloaded in the linux version?".

_Note 2_: The `solareclipser` user needs to install and _register_ SOLAR, that might be the most annoying part of the installation process.
Please see the Installation [section](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#installation) of the tutorial vignette.

### Quick start

Please see the vignette [minimal.html](http://ugcd.github.io/solareclipser/vignettes/minimal.html).


### Citation

To cite the `solareclipser` package in publications use:

```
  Ziyatdinov et al., solareclipser: an R interface to SOLAR for variance
  component analysis in pedigrees, Bioinformatics (2016)
```

A BibTeX entry for LaTeX users is:

```
  @article{ziyatdinov2016solareclipser,
    title = {solareclipser: an R interface to SOLAR for variance component analysis in pedigrees},
    author = {Andrey Ziyatdinov and Helena Brunel and Angel Martinez-Perez and Alfonso Buil and Alexandre Perera and Jose Manuel Soria},
    year = {2016},
    publisher = {Oxford Univ Press},
    journal = {Bioinformatics},
    pages = {btw080},
    url = {http://bioinformatics.oxfordjournals.org/content/32/12/1901},
  }
```

### Rationale

The rationale behind the `solareclipser` software:

* do not automate things in R, which `SOLAR` has already automated
    * call `SOLAR` by `system` R base function passing options/settings to SOLAR as parameters
* make it more R self-content and independent on other programs
    * phenotypes format as `data.frame`
    * make use of R plots like plotting pedigrees
    * make use of parallelization insfrastructure available in R
    * do not rely on tcl  scripts anymore
* get rid of `salamboR` artifacts (ancestor of `solareclipser`)
    * GAIT-specific functions
    * interface with other programs than `SOLAR`
    * lost version-control traces
    * dependence on old-school code from previous mantainers
    * dependence on (many) tcl scripts
* get rid of `SOLAR` artifacts
    * store models/phenos in folders/files
* make use of github infrastructure: collaborative coding, issues, gh-pages, etc
 
On the SOLAR side:

* Designed for the family-based studies (HHID, PROBND, FAMID descriptors of individuals)
    * support of extended pedigrees
* Stable routines for optimization of VC models
* Advanced polygenic models
    * support of multivariate models
    * liability threshold model (probit)
    * LRT applied to both covariates and variance components
* Elaborated linkage models
    * Multi-pass
    * Multivariate
    * Adjustment of LOD scores
* Association models
    * Speed-up based on residuals
* Advanced VC models (custom scripts)
    * Sex-specificity 
    * Longitudinal

On the R side:

* Interactive environment for data manipulation
* Graphics
    * Plot residuals, QQ-plot, Manhattan plot
* Parallel computing

## FAQ

Q1: I use GEMMA for my qtl mapping analysis. Are there any special reasons to switch to SOLAR / `solareclipser`?

A1: That depends on the type of analysis. Some of the features SOLAR has and GEMMA does not:

* GEMMA allows for a single random effect (the additive genetic), while SOLAR is flexible in the number of such effects.
   * By default, a SOLAR model has (1) the additive genetic (it is compulsory) (ID/FA/MO fields); and (2) the house-hold grouping (HHID field) effect. The `solareclipser` package can be easily used to specify such models with the two random effects.
   * Note that if HHID variable is not measured in a sample under study, SOLAR and GEMMA perform the same modeling.
   * If one is interested in a more general model with >2 random effects, then `solareclipser` can **not** help here and the only way is to deal with SOLAR tcl scripts, etc.
* SOLAR has the liability threshold (probit) model for binary traits, while GEMMA seems not.
  Consequently, a bi- or multi-variate trait model, where some traits are binary, is easy to define in SOLAR
  (SOLAR computes conditional likelihoods).
* Generally, SOLAR might be more stable and robust to ambiguous models in the case of multi-trait analysis,
  because of many years of improvement in this direction.
  
Q2: I was wondering whether there was any support for user-defined omega functions, or any slightly more complicated variance components models (e.g. to allow different variance components for males and females, or to report results as variance components rather than proportions of variance)?  If not, are there any plans to incorporate these in future? 

A2: The `solareclipser` package supports the three basic models: `polygenic`, `multipoint` and `mga` commands of SOLAR. 

In future versions of `solareclipser`, we are **not** going to include any complex models. In this sence, `solareclipser` is a partial port of SOLAR into R.

Also note that GxE tcl scripts are not publicly available. In our group, we did some specific analyses like GxE with custom tcl scripts (requested from the SOLAR authors).

Q3: When I tried to install the package, this was the message I got:

```
> install.packages(“solareclipser”)

Warning message:
package ‘solareclipser’ is not available (for R version 3.2.2)
``` 

Could you please let me know what I have done wrong and correct way to install the package?

A3: I suspect that you might have the Windows system, for which `solareclipser` is not supported. That is because of SOLAR.

Please see the installation notes on [http://ugcd.github.io/solareclipser/vignettes/tutorial.html#installation](http://ugcd.github.io/solareclipser/vignettes/tutorial.html#installation). 

## SOLAR references

* The new [SOLAR web page](http://solar-eclipse-genetics.org/) (SOLAR-Eclipse)
    * The old [SOLAR web page at txbiomedgenetics.org](http://solar.txbiomedgenetics.org/) (depreciated)
* [Appendix 1. SOLAR Command Descriptions](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html)
