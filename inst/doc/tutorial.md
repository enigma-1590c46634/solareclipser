Tutorial on `solareclipser` R package {#tutorial-on-solareclipser-r-package .title .toc-ignore}
=====================================

#### Andrey Ziyatdinov, Helena Brunel {#andrey-ziyatdinov-helena-brunel .author}

#### 2024-05-10 {#section .date}

::: {#cb1 .sourceCode}
``` {.sourceCode .r}
library(solareclipser)
```
:::

::: {#introduction .section .level1}
Introduction
============

This vignette is a tutorial on the R package `solareclipser`. The
document contains a brief description of the main statistical models
(polygenic, association and linkage) implemented in `SOLAR` and
accessible via `solareclipser`, installation instructions for both
`SOLAR` and `solareclipser`, reproducible examples on synthetic data
sets available within the `solareclipser` package. The examples are
intended to cover all the necessary steps in the quantitative genetic
analysis with `SOLAR`: data loading, data exploration, running
statistical models and interpreting the obtained results by printing,
summarizing and plotting the models in R.

In human genetics research, the study of a complex disorder is targeted
to discover its underlying genetic architecture. Whereas simple
Mendelian diseases are results of a mutation in a single gene, for
example [Huntington's
disease](http://en.wikipedia.org/wiki/Huntington%27s_disease), complex
disorders or phenotypes are the product of the action of multiple genes
and environmental factors as well as their interactions. For instance,
obesity is a complex disease which genetic basis might be better
dissected via measurable quantitative variation in phenotypes, such as
body-mass index, closely related to the disease risk.

One of the study designs to identify the quantitative trait loci (QTLs)
that affect complex phenotypes is based on examination of related
individuals grouped by families (pedigrees). Pedigree-based designs
allow for both QTL linkage and association mappings. In addition, QTL
analysis with pedigrees possess more power to detect rare variants,
permits explicit control for population stratification in association
studies, and suit for examination of parent-of-origin and
sex-specificity effects.

The analysis of the pedigree data is traditionally based on variance
component (VC) models, also known as linear mixed models. `SOLAR` is one
of the well-known and widely used tools that implements the variance
component models and other routines needed for the QTL analysis in
extended pedigrees [\[\@Almasy1998\]]{.citation}
(<http://solar.txbiomedgenetics.org/>). The name `SOLAR` stands for
Sequential Oligogenic Linkage Analysis Routines. The analyses available
in `SOLAR` include customized quantitative genetic analysis, advanced
linkage analysis, SNP association analysis (QTN, QTLD, and MGA). Other
operations are implemented for calculation of marker-specific or
multipoint identity-by-descent (IBD) matrices in pedigrees of arbitrary
size and complexity. The linkage analysis allows for multiple
quantitative traits and/or discrete traits which might involve multiple
loci (oligogenic analysis), dominance effects, household effects, and
interactions.

There are other tools similar to `SOLAR`, such as `Mendel`, `MMM`,
`GEMMA` and `FAST-LMM`, but the R package `solareclipser` presented here
has been especially developed to work with `SOLAR`.

`solareclipser` provides an interface from the program `SOLAR` to the
environment for statistical computing `R`. On one side, main `SOLAR`
procedures can be called in `R` by one or a few commands from the
`solareclipser` package, that makes the analysis more user-friendly with
the same power of computations implemented in `SOLAR`. On the other
side, the `solareclipser` user benefits from the infrastructure
available with the R environment and its packages, such as interactive
data manipulation, data visualization and parallel computation.
:::

::: {#methodology .section .level1}
Methodology
===========

::: {#linear-mixed-model .section .level2}
Linear mixed model
------------------

Consider that a particular trait, such as body-mass index, is observed
in [\\(n\\)]{.math .inline} individuals grouped by families or clusters
of related individuals. A [\\(n \\times 1\\)]{.math .inline} vector
[\\(y\\)]{.math .inline} contains the phenotypic values measured for the
trait. One assumes that these observations are described adequately by a
linear model with a [\\(p \\times 1\\)]{.math .inline} vector of fixed
effects ([\\(\\beta\\)]{.math .inline}) and a [\\(q \\times 1\\)]{.math
.inline} vector of random effects ([\\(u\\)]{.math .inline})
[\[\@Lynch1998, Chapter 26, pp. 746\]]{.citation}. In matrix form,

[\\\[ y = X \\beta + Z u + e \\\]]{.math .display}

where [\\(X\\)]{.math .inline} and [\\(Z\\)]{.math .inline} are
respectively [\\(n \\times p\\)]{.math .inline} and [\\(n \\times
q\\)]{.math .inline} incidence matrices, and [\\(e\\)]{.math .inline} is
a [\\(n \\times 1\\)]{.math .inline} vector of residual deviations
assumed to be distributed independently of the random effects. Denote
the [\\(n \\times n\\)]{.math .inline} covariance matrix for the vector
[\\(e\\)]{.math .inline} by [\\(R\\)]{.math .inline} and the the [\\(q
\\times q\\)]{.math .inline} covariance matrix for the vector
[\\(u\\)]{.math .inline} by [\\(G\\)]{.math .inline}.

The first element of the vector [\\(\\beta\\)]{.math .inline} is
typically the population mean, and other elements may be age, gender,
diet and so on. The elements of the incidence matrices are usually equal
to 0 or 1, that indicate the effect contribution to the individuals's
phenotype. The model is referred to as a *mixed model*, because this
model jointly accounts for fixed and random effects.

Giving the two terms of the random effects [\\(u\\)]{.math .inline} and
[\\(e\\)]{.math .inline}, the first term typically accounts for the
contribution from random genetic effects such as additive genetic
values. The second term accounts for the residual variance, and one
usually assumes that residual errors have constant variance and are
uncorrelated, so that [\\(R\\)]{.math .inline} is a diagonal matrix,
with [\\(R = \\sigma\_e\^2 I\\)]{.math .inline}.

Two more equations for the means and variances of the component vectors
of the mixed model can be derived on the basis of its definition. Since
[\\(E(u) = E(e) = 0\\)]{.math .inline} by definition,

[\\\[ E(y) = X \\beta \\\]]{.math .display}

Excluding the difference among individuals due to fixed effects and
calling the assumption that [\\(u\\)]{.math .inline} and [\\(e\\)]{.math
.inline} are uncorrelated, the covariance matrix for the vector of
observations [\\(y\\)]{.math .inline} is

[\\\[ Var(y) = V = Z G Z\^{T} + R \\\]]{.math .display}

Now the mixed model can be rewritten in a compact form

[\\\[ y = X \\beta + Z u + e\\\]]{.math .display}

[\\\[ \\mbox{ where } \\mbox{ } \\mbox{ } u \\sim (0, G) \\mbox{ and } e
\\sim (0, R) \\\\ \\mbox{ implying } \\mbox{ } \\mbox{ } y \\sim (X
\\beta, V) = (X \\beta, Z G Z\^T + R) \\\]]{.math .display}

For the mixed model, [\\(y\\)]{.math .inline}, [\\(X\\)]{.math .inline}
and [\\(Z\\)]{.math .inline} are observed variables, while
[\\(\\beta\\)]{.math .inline}, [\\(u\\)]{.math .inline}, [\\(R\\)]{.math
.inline} and [\\(G\\)]{.math .inline} are generally unknown variables.
Thus, the analysis consists of two complementary estimation procedures:
(1) estimation of the covariance matrices [\\(G\\)]{.math .inline} and
[\\(R\\)]{.math .inline}, and (2) estimation of the vectors of fixed and
random effects, [\\(\\beta\\)]{.math .inline} and [\\(u\\)]{.math
.inline}. The covariance matrices are usually functions of a few unknown
parameters or variance components. For instance, the covariance matrix
of the residuals errors typically has the form [\\(R = \\sigma\_e\^2
I\\)]{.math .inline}, where [\\(\\sigma\_e\^2\\)]{.math .inline} is the
parameter to be estimated.

The description of the two estimation procedures is out of scope of this
document, and the reader is referred to [\[\@Lynch1998, Chapter
26-27\]]{.citation}, where the necessary theoretical material and
accompanying examples are given. On the other hand, the implementation
of the estimation procedures in `SOLAR` is discussed elsewhere
[\[\@Almasy1998\]]{.citation}.

The next three sections will present `SOLAR` models, which follow the
general definition of linear mixed model given above, but
parametrization of the models and estimation procedures are implemented
in a `SOLAR`-specific way.
:::

::: {#polygenic-model-in-solar .section .level2}
Polygenic model in `SOLAR`
--------------------------

The polygenic model is a basis for any quantitative genetic analysis in
`SOLAR`. The model includes fixed effects ([\\(\\beta\\)]{.math
.inline}) and three random effects (the genetic additive effect
[\\(g\\)]{.math .inline}, the house-hold effect [\\(c\\)]{.math
.inline}, and the residual effect [\\(e\\)]{.math .inline}).

[\\\[y = X\\beta + g + c + e\\\]]{.math .display}

[\\\[\\mbox{ where } \\mbox{ } \\mbox{ } g \\sim (0, \\sigma\_g\^2 A)
\\mbox{ and } \\mbox{ } \\mbox{ } c \\sim (0, \\sigma\_c\^2 H) \\mbox{
and } \\mbox{ } \\mbox{ } e \\sim (0, \\sigma\_e\^2 I)\\\]]{.math
.display}

[\\\[ \\mbox{ implying } \\mbox{ } \\mbox{ } y \\sim (X\\beta, V) =
(X\\beta, \\sigma\_g\^2 A + \\sigma\_c\^2 H + \\sigma\_e\^2
I)\\\]]{.math .display}

The first [\\(g\\)]{.math .inline} term accounts for the additive
genetic effect with the covariance matrix [\\(A = 2 K\\)]{.math
.inline}, where [\\(K\\)]{.math .inline} is the kinship matrix. `SOLAR`
estimates the *empirical* kinship matrix based on the identifier (ID)
fields in the pedigree data, such as individual ID `id`, gender `sex`,
family ID `fam`, father ID `fa`, mother ID `mo`, and indicator of
identical individuals, e.g. monozygotic twins, `mztwin` [\[\@SOLARDOC,
Appendix 1,
\[file-pedigree\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#file-pedigree)
command\]]{.citation}.

The second [\\(c\\)]{.math .inline} term represents an annotated
environmental effect. It is named as a house-hold effect, but it is
generally related to the sharing of any non-genetic factor among
individuals. The house-hold effect is included in the analysis if the
pedigree data contain a house-hold ID field (the expected field name is
`hhid`). `SOLAR` automatically creates `Pedigree-Household Groups` based
on the `hhid` identifiers. That means that pedigrees and household
groups are merged into larger groups, which contain every individual in
the same pedigree as well as every individual that has the same `hhid`
as any individual in the pedigree [\[\@SOLARDOC, Section 9.3 \[Household
Group
Analysis\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/09.chapter.html\#household)\]]{.citation}.
As a consequence of this grouping operation by `SOLAR`, the covariance
matrix [\\(H\\)]{.math .inline} of the [\\(c\\)]{.math .inline} random
effect contains the information related to the pedigree-household groups
rather than relationships given by the `hhid` field.

The third [\\(e\\)]{.math .inline} term is the residual error with the
diagonal covariance matrix. It is important to note that the `residual`
command in `SOLAR` [\[\@SOLARDOC, Appendix 1,
\[residual\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#residual)
command\]]{.citation} returns the covariate-free value of a trait under
study, expressed as [\\(y - X \\beta\\)]{.math .inline}, and the
returned value is not the residual vector [\\(e\\)]{.math .inline} given
in the equation of the polygenic model.

The polygenic term [\\(g\\)]{.math .inline} represents the additive
genetic effect by default. The incorporation of the dominance effect
into `SOLAR` models is also possible, but it is attributed to advanced
modeling topics [\[\@SOLARDOC, Section 9.4 \[Dominance
Analysis\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/09.chapter.html\#dominance)\]]{.citation}.
The `SOLAR`'s developers suggest that the dominance variance component
is unlikely useful in quantitative trait mapping, particularly in human
genetic analysis, as only bi-lineal relatives contribute to this
component.

`SOLAR` is able to work with discrete traits, which possible values are
coded as two consecutive integers, for example, `0` and `1`
[\[\@SOLARDOC, Appendix 1,
\[discrete-notes\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#discrete-notes)\]]{.citation}.
`SOLAR` uses a *liability threshold model* to handle discrete traits
[\[\@Duggirala1997\]]{.citation}, and more details about implementation
and practical notes are available in [\[\@SOLARDOC, Section 9.1
\[Discrete
Traits\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/09.chapter.html)\]]{.citation}.
`SOLAR` authors recommend to work with quantitative traits whenever it
is possible. The user is able to disable the transformation of discrete
traits, if it is necessary.

The likelihood-based statistics, in particular Likelihood Ratio Tests
(LRTs), is used for the inference of the linear mixed models implemented
in `SOLAR`. The following table shows the alternative models used to
test the significance of the random effects [\\(g\\)]{.math .inline} and
[\\(c\\)]{.math .inline}.

  ------------------------------------------------------------------------------------
  Alternative model        Free parameters               Restricted parameters
  ------------------------ ----------------------------- -----------------------------
  No polygenic effect      [\\(\\beta\\)]{.math          [\\(\\sigma\_g\^2\\)]{.math
                           .inline},                     .inline} = 0
                           [\\(\\sigma\_c\^2\\)]{.math   
                           .inline},                     
                           [\\(\\sigma\_e\^2\\)]{.math   
                           .inline}                      

  No house-hold effect     [\\(\\beta\\)]{.math          [\\(\\sigma\_c\^2\\)]{.math
                           .inline},                     .inline} = 0
                           [\\(\\sigma\_g\^2\\)]{.math   
                           .inline},                     
                           [\\(\\sigma\_e\^2\\)]{.math   
                           .inline}                      
  ------------------------------------------------------------------------------------

Under the assumption of multivariate normality, the LRT statistic is
asymptotically distributed as a [\\(\\frac{1}{2}:\\frac{1}{2}\\)]{.math
.inline} mixture of of a [\\(\\chi\_1\^2\\)]{.math .inline} variable and
a point mass at the origin [\[\@Blangero2000\]]{.citation}. This
statement is applied to any LRT performed for parameters like
[\\(\\sigma\_g\^2\\)]{.math .inline}, which is greater or equal to `0`
by definition.

`SOLAR` takes the decision of rejecting the house-hold effect by
considering two factors, the test results and the value of
[\\(\\sigma\_c\^2\\)]{.math .inline} parameter. If value of
[\\(\\sigma\_c\^2\\)]{.math .inline} is non-zero, the house-hold effect
is retained in the polygenic model, nevertheless the effect is not
significant. The user has an option to force the rejection of the
house-hold effect if it is necessary.

The tests applied to fixed effects (covariates) are also performed by
means of LRTs. In the table given below the alternative model for the
test of significance of `age` covariate is shown.

  ------------------------------------------------------------------------------------
  Alternative model        Free parameters               Restricted parameters
  ------------------------ ----------------------------- -----------------------------
  No age effect            [\\(\\sigma\_g\^2\\)]{.math   [\\(\\beta\_{age}\\)]{.math
                           .inline},                     .inline} = 0
                           [\\(\\sigma\_c\^2\\)]{.math   
                           .inline},                     
                           [\\(\\sigma\_e\^2\\)]{.math   
                           .inline}                      

  ------------------------------------------------------------------------------------

::: {#multivariate-polygenic-model-in-solar .section .level3}
### Multivariate polygenic model in `SOLAR`

In principle, one can estimate the multivariate polygenic model by
simply performing the univariate analysis in an appropriate way.
Consider two traits [\\(y\_1\\)]{.math .inline} and [\\(y\_2\\)]{.math
.inline} measured in each of [\\(n\\)]{.math .inline} individuals,
assuming that each trait follows the same univariate polygenic analysis.
Then the [\\(2 \\times n\\)]{.math .inline} dimensional column vector of
observations is denoted by the stack of univariate vectors, and the
polygenic model has the matrix form:

[\\\[ \\begin{pmatrix} y\_1 \\\\ y\_2 \\end{pmatrix}\_{2n \\times 1} =
\\begin{pmatrix} X\_1 & 0 \\\\ 0 & X\_2 \\end{pmatrix}\_{2n \\times nk}
\\begin{pmatrix} \\beta\_1 \\\\ \\beta\_2 \\end{pmatrix}\_{2n \\times 1}
+ \\begin{pmatrix} g\_1 \\\\ g\_2 \\end{pmatrix}\_{2n \\times 1} +
\\begin{pmatrix} e\_1 \\\\ e\_2 \\end{pmatrix}\_{2n \\times 1}
\\\]]{.math .display}

[\\\[ G\_{2n \\times 2n} = C\_{2 \\times 2} \\otimes A\_{n \\times n} =
\\begin{pmatrix} c\_{11} A & c\_{12} A \\\\ c\_{21} A & c\_{21} A
\\end{pmatrix} = \\begin{pmatrix} \\sigma\_{1}\^2 A & \\rho \\sigma\_{1}
\\sigma\_{2} A \\\\ \\rho \\sigma\_{1} \\sigma\_{2} A & \\sigma\_{2}\^2
A \\end{pmatrix}\_{g} \\\]]{.math .display}

[\\\[ R\_{2n \\times 2n} = E\_{2 \\times 2} \\otimes I\_{n \\times n} =
\\begin{pmatrix} e\_{11} I & e\_{12} I \\\\ e\_{21} I & e\_{21} I
\\end{pmatrix} = \\begin{pmatrix} \\sigma\_{1}\^2 I & \\rho \\sigma\_{1}
\\sigma\_{2} I \\\\ \\rho \\sigma\_{1} \\sigma\_{2} I & \\sigma\_{2}\^2
I \\end{pmatrix}\_{e} \\\]]{.math .display}

The model can be generalized to [\\(k\\)]{.math .inline} traits, and
other random effects such as house-hold can be added in the similar way.
The *Kronecker product* notation [\\(\\otimes\\)]{.math .inline} is
convinient for denoting the covariance matrices. More details on the
multivariate analysis are given in [\[\@Lynch1998, Chapter 26, Multiple
Traits, pp. 774\]]{.citation}.

The multivariate analysis implemented in `SOLAR` allows for any number
of traits, but the tests for significance of the correlation
coefficients are performed only for the bivariate models. The following
table describes the main LRTs available in `SOLAR` for the bivariate
analysis.

  -----------------------------------------------------------------------------
  Alternative model       Free parameters               Restricted parameters
  ----------------------- ----------------------------- -----------------------
  No genetic correlation  [\\((\\sigma\_1\^2,           [\\(\\rho\_g =
                          \\sigma\_2\^2)\_g\\)]{.math   0\\)]{.math .inline}
                          .inline}, [\\((\\sigma\_1\^2, 
                          \\sigma\_2\^2,                
                          \\rho)\_e\\)]{.math .inline}  

  Complete pleiotropy     [\\((\\sigma\_1\^2,           [\\(\\rho\_g =
                          \\sigma\_2\^2)\_g\\)]{.math   1\\)]{.math .inline}
                          .inline}, [\\((\\sigma\_1\^2, 
                          \\sigma\_2\^2,                
                          \\rho)\_e\\)]{.math .inline}  

  No environmental        [\\((\\sigma\_1\^2,           [\\(\\rho\_e =
  correlation             \\sigma\_2\^2,                0\\)]{.math .inline}
                          \\rho)\_g\\)]{.math .inline}, 
                          [\\((\\sigma\_1\^2,           
                          \\sigma\_2\^2)\_e\\)]{.math   
                          .inline}                      
  -----------------------------------------------------------------------------

It is important to note that the transformation of discrete traits in
the multivariate analysis is automatically performed by `SOLAR`.
:::
:::

::: {#linkage-model-in-solar .section .level2}
Linkage model in `SOLAR`
------------------------

`SOLAR` introduces sequential oligogenic linkage analysis, as stated in
its name abbreviation. Consider the linkage model with a single QTL,
which is described by another random effect [\\(q\\)]{.math .inline}.

[\\\[y = X\\beta + g + c + q + e\\\]]{.math .display}

[\\\[ \\mbox{ where } \\mbox{ } \\mbox{ } g \\sim (0, \\sigma\_g\^2 A)
\\mbox{ and } \\mbox{ } \\mbox{ } c \\sim (0, \\sigma\_c\^2 H) \\mbox{
and } \\mbox{ } \\mbox{ } q \\sim (0, \\sigma\_q\^2 M) \\mbox{ and }
\\mbox{ } \\mbox{ } e \\sim (0, \\sigma\_e\^2 I) \\\]]{.math .display}

[\\\[ \\mbox{ implying } \\mbox{ } \\mbox{ } y \\sim (X\\beta, V) =
(X\\beta, \\sigma\_g\^2 A + \\sigma\_c\^2 H + \\sigma\_q\^2 M +
\\sigma\_e\^2)\\\]]{.math .display}

The QTL under study [\\(q\\)]{.math .inline} corresponds to a specific
position in the genome and an identity-by-descent (IBD) matrix. The
genome-wide linkage analysis assumes that the highly informative
markers, such as micro-satellites, are available in the study, in order
to estimate the two-point or multipoint IBD matrices by `SOLAR` or other
specific program like `MERLIN`. Once the IBD matrices are calculated,
the genome-wide scan is conducted by sequentially running the linkage
model for each genome position.

The alternative model for testing the significance of a single QTL
effect is the same as the polygenic model.

  -----------------------------------------------------------------------------------
  Alternative model       Free parameters               Restricted parameters
  ----------------------- ----------------------------- -----------------------------
  No linkage              [\\(\\beta, \\sigma\_g\^2,    [\\(\\sigma\_q\^2\\)]{.math
                          \\sigma\_c\^2\\)]{.math       .inline} = 0
                          .inline},                     
                          [\\(\\sigma\_e\^2\\)]{.math   
                          .inline}                      

  -----------------------------------------------------------------------------------

Traditionally, the logarithm (base 10) of the odds (LOD) values are used
to score the QTLs. If the LRT statistics for the linkage model is
denoted as [\\(\\Lambda\\)]{.math .inline}, then the correspondent LOD
score is equal to [\\(\\Lambda / (2 ln(10))\\)]{.math .inline}.

`SOLAR` also allows for multi-pass linkage scan, that means that any
further scan can be conditioned to QTLs found in the previous scans. For
example, the second-pass scan conditioned to the previously found QTL
[\\(q\_1\\)]{.math .inline} has the form:

[\\\[y = X\\beta + g + c + q\_1 + q + e\\\]]{.math .display}

The second-pass LRT test now looks like in the table below.

  ----------------------------------------------------------------------------------------
  Alternative model       Free parameters                    Restricted parameters
  ----------------------- ---------------------------------- -----------------------------
  No linkage              [\\(\\beta, \\sigma\_g\^2,         [\\(\\sigma\_q\^2\\)]{.math
                          \\sigma\_c\^2\\)]{.math .inline},  .inline} = 0
                          [\\(\\sigma\_{q\_1}\^2\\)]{.math   
                          .inline},                          
                          [\\(\\sigma\_e\^2\\)]{.math        
                          .inline}                           

  ----------------------------------------------------------------------------------------

The linkage model is extended to deal with many traits or discrete
traits in the same way as for the polygenic model.

Additional information related to the linkage modeling in `SOLAR` can be
found in related publications [\[\@Almasy1998\]]{.citation},
[\[\@Almasy1997\]]{.citation}, [\[\@Williams1999\]]{.citation}, the
`SOLAR` tutorial [\[\@SOLARDOC, Chapter 3
\[Tutorial\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/03.chapter.html)\]]{.citation},
and the `SOLAR` manual, for example, [\[\@SOLARDOC, Apendix 1,
\[multipoint\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#multipoint)
command\]]{.citation}.
:::

::: {#association-model-in-solar .section .level2}
Association model in `SOLAR`
----------------------------

A basic association model in `SOLAR` is extended from the polygenic
model by simply adding a new covariate (fixed effect) `snp`.

[\\\[ y = X\\beta + \\beta\_{snp} \* snp + g + c + e\\\]]{.math
.display}

[\\\[\\mbox{ where } \\mbox{ } \\mbox{ } g \\sim (0, \\sigma\_g\^2 A)
\\mbox{ and } \\mbox{ } \\mbox{ } c \\sim (0, \\sigma\_c\^2 H) \\mbox{
and } \\mbox{ } \\mbox{ } e \\sim (0, \\sigma\_e\^2 I)\\\]]{.math
.display}

[\\\[ \\mbox{ implying } \\mbox{ } \\mbox{ } y \\sim (X\\beta, V) =
(X\\beta, \\sigma\_g\^2 A + \\sigma\_c\^2 H + \\sigma\_e\^2
I)\\\]]{.math .display}

In the genome-wide association scan `SOLAR` passes through SNPs variants
under the study one by one according to the equation given above. The
p-value of association between trait [\\(y\\)]{.math .inline} and SNP
[\\(snp\\)]{.math .inline} is again based on LRT test with the
alternative model outlined in the table below. This alternative model is
the same as the polygenic model.

  -----------------------------------------------------------------------------
  Alternative model       Free parameters               Restricted parameters
  ----------------------- ----------------------------- -----------------------
  No association          [\\(\\beta, \\sigma\_g\^2,    [\\(\\beta\_{snp} =
                          \\sigma\_c\^2\\)]{.math       0\\)]{.math .inline}
                          .inline},                     
                          [\\(\\sigma\_e\^2\\)]{.math   
                          .inline}                      

  -----------------------------------------------------------------------------

More details on implementation of the association analysis in `SOLAR`
can be found in the `SOLAR` manual, for example, [\[\@SOLARDOC, Apendix
1,
\[mga\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#mga)
and
\[snp\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix\_1\_text.html\#snp)
commands\]]{.citation}.
:::
:::

::: {#software .section .level1}
Software
========

The idea of the proposed software (the R package `solareclipser`) is to
perform the calculus in `SOLAR` at the back-end and to offer the user
interface in R at the front-end. The table given below shows a basic
description of the software via implementation of the main three models:
polygenic, association and linkage.

  -------------------------------------------------------------------------------------
  Model          `SOLAR`        `solareclipser`     `solareclipser` S3   Parallel
                 command        function            classes of returned  computation
                                                    object               
  -------------- -------------- ------------------- -------------------- --------------
  Polygenic      `polygenic`    `solarPolygenic`    `solarPolygenic`     None

  Association    `mga`          `solarAssoc`        `solarAssoc`,        Automatic or
                                                    `solarPolygenic`     custom (by SNP
                                                                         file)

  Linkage        `multipoint`   `solarMultipoint`   `solarMultipoint`,   Custom (by
                                                    `solarPolygenic`     chromosome)
  -------------------------------------------------------------------------------------

One can see how the software benefits from both `SOLAR` and R sides. On
the side of `SOLAR`, the `SOLAR` commands `polygenic`, `mga` and
`multipoint` do all the hard job on estimating the variance component
models with control via a set of parameters and settings. On the side of
`R`, the `solareclipser` functions `solarPolygenic`, `solarAssoc` and
`solarMultipoint` prepare the input data for the analysis, pass the
arguments to the commands and set up other necessary settings, and
finally read the output `SOLAR` files and store the results into R
objects of S3 classes in R.

The main advantage of running the `SOLAR` models in `solareclipser`
instead of `SOLAR` is that the analysis is conduced by only one function
call, while the same analysis in `SOLAR` could require typing more
commands with appropriate options and settings (to be learned from the
`SOLAR` manual). The `solareclipser` functions just automate the
analysis flow and keep configuration options via the R function
arguments, if such configuration is required. Another advantage of using
the R package `solareclipser` is that the results are stored in objects
of especially designed S3 classes with associated methods like print,
summary, plot and others.

The last column of the table given above contains the information
related to parallel computation available for the models. The
implementation of parallel calculations is straightforward in
`solareclipser` package, since the association and linkage analyses are
*embarrassingly parallel* problems (easily divisible into different
identical independent task) and the R environment offers a number of
packages with parallel interfaces.

More information about the `solareclipser` functions is available in the
software manual of the package, and the practical guide of using the
functions is presented in the next section dedicated to examples.
Meantime, some points should to be stressed to describe some other
features and drawbacks of the `solareclipser` package.

-   The `solareclipser` functions, by default, do not create any
    directories visible to the user, though the analysis flow in `SOLAR`
    is heavily based on storing the data in directories and files.
    However, the input phenotype and genotype data might be stored in
    CSV data format and passed to the functions. The user can specify a
    special directory (`dir` argument for the `solarPolygenic`,
    `solarAssoc` and `solarMultipoint` functions), where the `SOLAR`
    output will be saved in addition to the output R object.
-   The `tcl` scripts used by `SOLAR` for automation are not visible to
    the user. The routine of executing either `SOLAR` separate commands
    or a block of commands in a `tcl` script is internally organized in
    the functions of `solareclipser`.
-   The work on parsing the output of `SOLAR` analysis is implemented in
    the functions in `solareclipser`, and the user does not need to
    manually look for the analysis results.
-   The `SOLAR` commands used to compute the analysis and some `SOLAR`
    output files are saved in the returned objects. This feature might
    be necessary for those users who are familiar with `SOLAR`. It is
    also helps to go back to `SOLAR` and repeat or extend the analysis
    in the `SOLAR` environment, if it is necessary.
-   The package `solareclipser` offers only three analyses via the most
    powerful `SOLAR` commands `polygenic`, `mga` and `multipoint`. If a
    more specific genetic analysis is required, the user needs to go
    back to `SOLAR`.

::: {#installation .section .level2}
Installation
------------

The installation process has three steps: install `SOLAR`, register
`SOLAR` and install `solareclipser`.

::: {#solar-installation .section .level3}
### `SOLAR` installation

The `SOLAR` distribution is available on the official web page
<http://solar.txbiomedgenetics.org/> for the most popular operating
systems (Linux, Mac OS, Windows).

In the case of
[Linux](http://solar.txbiomedgenetics.org/solarlinux.html), the
installation files are compressed in `solar_linux.tar.gz` archive. Once
the archive is downloaded to a local machine of the user and the
installation files are unpacked, one can find the `README` file and
`install_solar` script. In particular, the script copies the necessary
library, binary and documentation files into an installation directory
defined by the user. A special bash script called `solar` is also
created and copied to the system path defined by the user,
e.g. `/usr/local/bin/solar`, so that the program can be executed by
simply typing `solar` in the terminal.

An example of the installation command might be the following.

    ./install_solar /opt/appl/solar/7.6.6 /usr/local/bin solar

As a result of execution of the command, the library, binary and
documentation files of `SOLAR` version `7.6.6` will be copied to
`/opt/appl/solar/7.6.6` directory, and `solar` execution script will be
copied to `/usr/local/bin` system directory.

After the installation process the new user is required to obtain a
registration key, as it is explained by the message that appears when
`SOLAR` program starts. The registration key has to be stored in
`.solar_reg` file in the home directory of the user.

Once both installation and registration processes are finished,
something similar to the following message will be displayed when
`SOLAR` is executed.

    SOLAR Eclipse version 7.6.6 (Experimental), last updated on December 11, 2014
    Copyright (c) 1995-2014 Texas Biomedical Research Institute
    Enter help for help, exit to exit, doc to browse documentation.

    solar>
:::

::: {#installation-problems-with-solar .section .level3}
### Installation problems with SOLAR

**SOLAR version for Windows**. Currently, SOLAR does not have an
appropriated version for Windows. As stated in the official web page of
SOLAR <http://solar.txbiomedgenetics.org/solarwindows.html>:

> Because SOLAR requires a Unix OS or equivalent, we have prepared a
> VMware Virtual Machine containing SOLAR plus the Xubuntu edition of
> linux which runs on Windows XP or later using the free VMware Player.
> Beware this requires a very large multi 1.4 Gb download, and will
> require more than 10 gigabytes free disc space during installation.

This is a clear limitation of `solareclipser` package that comes from
the dependency on SOLAR. Starting from version 3.\*, `solareclipser` is
not supported for Windows.
[DESCRIPTION](https://github.com/ugcd/solarius/blob/master/DESCRIPTION)
file has a line `OS_type: unix`.
:::

::: {#solareclipser-installation .section .level3}
### `solareclipser` installation

`solareclipser` package is available in the CRAN repository for R
packages <http://cran.r-project.org/web/packages/solareclipser/>, so it
can be install by calling the usual R function `install.packages`.

    install.packages(solareclipser)

The code repository of `solareclipser` package is hosted in the public
GitHub repository of the UGCD group <https://github.com/ugcd/solarius>.
The R package `devtools` allows to install the development version of
the package hosted at GitHub. The installation commands are the
following.

    library(devtools)
    install_github("ugcd/solareclipser")
:::

::: {#installation-problems-with-solareclipser .section .level3}
### Installation problems with solareclipser

`solareclipser` has several R packages as dependencies (that can be seen
in
[DESCRIPTION](https://github.com/ugcd/solareclipser/blob/master/DESCRIPTION)
file). Installation of these packages should be straightforward using
the standard R function `install.packages`.

In the case of problems in the installation process, web pages of the
packages are the main source of information.

Imported packages:

-   <https://github.com/hadley/plyr> for `plyr`.
-   <https://github.com/hadley/ggplot2> for `ggplot2`.
-   <https://github.com/Rdatatable/data.table> for `data.table`.
    -   See also
        [CRAN\_Release.cmd](https://github.com/Rdatatable/data.table/blob/master/CRAN_Release.cmd)
        document.

**Installation of data.table R package**. If the user experiences
problems installing `data.table` package, we recommend to upgrade the R
environment to the latest stable version.
:::
:::

::: {#solareclipser-r-package-structure .section .level2}
`solareclipser` R package structure
-----------------------------------

`solareclipser` is organized as a standard R package, and its structure
can be seen in the code repository web site
<https://github.com/ugcd/solarius>.

-   [R](https://github.com/ugcd/solarius/tree/master/R) directory
    contains `*.R` files with functions, classes and methods.
-   [data](https://github.com/ugcd/solarius/tree/master/data) directory
    has data sets distributed with the package and stored in `*.RData`
    files.
-   `man` directory includes documentation organized in `Rd` files
    (**R** **d**ocumentation files).
-   [inst](https://github.com/ugcd/solarius/tree/master/inst) directory
    includes directories, where some extra data in any format can be
    stored for the R package. Sub-directories in `inst` directory will
    be copied to the installation directory of the package and, thus,
    will be available within the package via `system.file` function.
    -   [inst/examples](https://github.com/ugcd/solarius/tree/master/inst/examples)
        has `*.R` scripts with code examples.
    -   [inst/extdata](https://github.com/ugcd/solarius/tree/master/inst/extdata)
        mostly contains `SOLAR` data files for data sets of the package.
    -   [inst/tests](https://github.com/ugcd/solarius/tree/master/inst/tests)
        has `*.R` scripts for testing units during development of the
        package (R package `testthat`).
-   [vignettes](https://github.com/ugcd/solarius/tree/master/vignettes)
    directory includes vignettes, which are known as a long-form
    documentation for an R package. Typically, vignettes has a form of a
    technical report or an academic paper, that describes the problem
    that your package is designed to solve. Files with `Rmd` extension
    are the source files of the vignettes.
:::

::: {#documentation .section .level2}
Documentation
-------------

The source files of documentation are stored in `man` directory as `Rd`
files. R package `staticdocs` converts this documentation into html
pages. The documentation for `solarius` is published on
<http://ugcd.github.io/solarius/doc/>.

In the R environment the help for a function of the package, for example
`solarPolygenic` function, can be obtained by the standard
`?solarPolygenic` R comand. The content of the documentation for
`solarPolygenic` function is the same as on
<http://ugcd.github.io/solarius/doc/solarPolygenic.html>.
:::

::: {#vignettes .section .level2}
Vignettes
---------

The vignetts are created by `rmarkdown` R package and published as html
pages on the web.

-   <http://ugcd.github.io/solarius/vignettes/tutorial.html>: The
    vignette presents a tutotial of the package.
-   <http://ugcd.github.io/solarius/vignettes/minimal.html>: The
    vignette gives a minial example of the package. It includes three
    analyses, polygenic, association and linkage, applied to `dat30`
    data set.
-   <http://ugcd.github.io/solarius/vignettes/modelsGAIT1.html>: The
    vignette reports results on the real family-based data set GAIT1.
    The traits under study are the FXI levels in blood, the BMI and the
    Thrombosis affection. Although the original calculus were conducted
    with `SOLAR`, the code presented in the vignette is based on
    `solareclipser`. The produced results are the same as those reported
    in the published articles.
:::

::: {#load-solareclipser .section .level2}
Load `solareclipser`
--------------------

To start using `solareclipser` package, the user needs to load the
package with `library` function.

::: {#cb6 .sourceCode}
``` {.sourceCode .r}
library(solareclipser)
```
:::

The version of `solareclipser` loaded here is:

::: {#cb7 .sourceCode}
``` {.sourceCode .r}
packageVersion("solareclipser")
#> [1] '0.3.3'
```
:::

Other R packages will be needed for the tutorial.

::: {#cb8 .sourceCode}
``` {.sourceCode .r}
library(plyr)
```
:::
:::

::: {#minimal-example .section .level2}
Minimal example
---------------

The following code shows how the main functions (`solarPolygenic`,
`solarMultipoint` and `solarAssoc`) of the package work.

::: {#cb9 .sourceCode}
``` {.sourceCode .r}
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
:::

The user can think of this block of code as a quick start in using
`solareclipser`. The material given further will provide more details
and all necesseary information on the usage of `solareclipser`.
:::

::: {#data-sets .section .level2}
Data sets
---------

The genome-wide QTL analysis typically involves large-scale data sets.
Data sets distributed within `solareclipser` package contain relatively
small simulated data, which are used for reproducible demonstrations of
how functions of the package work. There are three data sets available
in `solareclipser`:

-   An adapted version of the simulated data from the R package `multic`
    [\[\@Andrade2006\]]{.citation} (polygenic and linkage models);
    -   A smaller subset of these data is given in `dat30` data set
        (polygenic, linkage and association models);
-   An adapted version of the simulated data `dat50` from the R package
    [FFBSKAT](http://mga.bionet.nsc.ru/soft/FFBSKAT/) (polygenic model
    and association model, both with custom kinship matrix);
-   GAW10 data distributed within `SOLAR` (polygenic and linkage
    models).

The simulated data set (in `SOLAR` format) from `multic` R package is
placed in
[inst/extdata/solarOutput](https://github.com/ugcd/solarius/tree/master/inst/extdata/solarOutput)
directory of `solareclipser` package. 1200 individuals in 200 families
are included in the study. The simulated phenotypes include two traits
`trait1` and `trait2` and two covariates `age` and `sex`. The two traits
possess a high genetic correlation. For the linkage analysis 6 IBD
matrices for chromosome 5 6 IBD matrices for chromosome 2 are available.
The six matrices of chromosome 5 are the original data from `multic`
package, while the six matrices of chromosome 2 are the exact copies of
those of chromosome 5. Data for the two chromosomes are needed to test
parallel calculation (separated by chromosome) of the linkage model.

A smaller subset of the simulated data from `multic` package is provided
in `dat30` data set (`dat30.RData` file in
[data](https://github.com/ugcd/solarius/tree/master/data) directory).
This smaller data set contains only 174 individuals in 29 families.

-   The phenotypes (`dat30` data frame) are the same as in the complete
    data described above.
-   A hundred of synthetic SNPs (`genocovdat30` matrix) were randomly
    generated to run the association analysis on these data (see the
    script
    [inst/scripts/generate-genocovdat30.R](https://github.com/ugcd/solarius/blob/master/inst/scripts/generate-genocovdat30.R)).
-   Annotation information (mapdat30 data frame) also was generated,
    mainly in order to plot the association results with Manhattan plot.

Documentation of `dat30` data set is available on the page
[dat30](http://ugcd.github.io/solarius/doc/dat30.html).

The simulated data from the R package `FFBSKAT` are stored in `dat50`
data set (`dat50.RData` file in
[data](https://github.com/ugcd/solarius/tree/master/data) directory).
This data set is adapted from `example.data` data set in the R package
`FFBSKAT`. The following variables are included.

-   `phenodata`: A data.frame with the phenotypes `trait`, `sex` and
    `age` given for 66 individuals (related and unrelated).
-   `genodata`: A matrix of 50 genetic variants (columns) given for 66
    individuals (rows). The genotypes are coded in the format such as
    `1/1`, `1/2` and `2/2`.
-   `genocovdata`: A matrix of covariates derived from the genotype data
    (the additive model). The covariates are coded as integers 0, 1
    and 2.
-   `snpdata`: A data.frame of annotation for genetic variants. The
    variables name of the variant `name`, chromosome `chrom`, gene name
    `gene` and position in bp `position`.
-   `kin`: A kinship matrix for 66 individuals.

Documentation of `dat50` data set is available on the page
[dat50](http://ugcd.github.io/solarius/doc/dat50.html).

The simulated data from the workshop GAW10 are described elsewhere
[\[\@SOLARDOC, Chapter 3
\[Tutorial\](http://helix.nih.gov/Documentation/solar-6.6.2-doc/03.chapter.html)\]]{.citation}.
:::
:::

::: {#examples .section .level1}
Examples
========

::: {#loading-data .section .level2}
Loading data
------------

To load both data sets `dat30` and `dat50`, the user needs to use the
standard R function `data` A quick view to the data tables might be
performed by means of `str` function.

::: {#cb10 .sourceCode}
``` {.sourceCode .r}
data(dat30)
str(dat30)
#> 'data.frame':    174 obs. of  10 variables:
#>  $ famid : int  1 1 1 1 1 1 2 2 2 2 ...
#>  $ id    : int  11 12 13 14 15 16 21 22 23 24 ...
#>  $ fa    : int  0 0 11 11 11 11 0 0 21 21 ...
#>  $ mo    : int  0 0 12 12 12 12 0 0 22 22 ...
#>  $ sex   : int  1 2 1 2 1 1 1 2 2 1 ...
#>  $ affect: int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ class : logi  NA NA NA NA NA NA ...
#>  $ trait1: num  11.96 7.1 10.32 9.76 9.46 ...
#>  $ trait2: num  13.58 5.37 6.4 8.98 9.21 ...
#>  $ age   : int  50 25 35 49 51 45 37 29 39 41 ...
str(genocovdat30)
#>  num [1:174, 1:100] 1.978 0.795 0.231 0.139 0.487 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:174] "11" "12" "13" "14" ...
#>   ..$ : chr [1:100] "snp_1" "snp_2" "snp_3" "snp_4" ...
str(mapdat30)
#> 'data.frame':    100 obs. of  4 variables:
#>  $ SNP : Factor w/ 100 levels "snp_1","snp_10",..: 1 12 23 34 45 47 48 49 50 2 ...
#>  $ chr : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ pos : num  2105324 2105467 2106094 2108138 2109262 ...
#>  $ gene: Factor w/ 12 levels "gene1","gene2",..: 1 1 1 1 1 1 1 1 1 2 ...
```
:::

::: {#cb11 .sourceCode}
``` {.sourceCode .r}
data(dat50)
str(phenodata)
#> 'data.frame':    66 obs. of  4 variables:
#>  $ id   : num  1 2 3 4 5 6 7 8 9 10 ...
#>  $ sex  : int  0 1 0 0 0 1 0 1 1 0 ...
#>  $ age  : int  80 77 56 44 75 79 75 82 77 76 ...
#>  $ trait: num  -1.763 -1.423 -0.805 0.268 -1.334 ...
str(genodata)
#>  chr [1:66, 1:50] "1/1" "1/1" "1/1" "1/2" "1/2" "1/1" "1/1" "1/1" "1/1" ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:66] "1" "2" "3" "4" ...
#>   ..$ : chr [1:50] "s1" "s2" "s3" "s4" ...
str(genocovdata)
#>  int [1:66, 1:50] 0 0 0 1 1 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:66] "1" "2" "3" "4" ...
#>   ..$ : chr [1:50] "s1" "s2" "s3" "s4" ...
```
:::

`loadMulticPhen` function in `solareclipser` package loads the complete
version of the data from the R package `multic` (1200 individuals).

::: {#cb12 .sourceCode}
``` {.sourceCode .r}
mdat <- loadMulticPhen()
str(mdat)
#> 'data.frame':    1200 obs. of  10 variables:
#>  $ ID    : chr  "1001" "1002" "1003" "1004" ...
#>  $ affect: int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ class : logi  NA NA NA NA NA NA ...
#>  $ trait1: num  9.61 11.5 11.05 12.37 7.97 ...
#>  $ trait2: num  9.61 11.66 13.62 13.88 7.29 ...
#>  $ age   : int  45 37 45 48 31 41 50 38 43 54 ...
#>  $ FAMID : chr  "100" "100" "100" "100" ...
#>  $ FA    : chr  "0" "0" "1001" "1001" ...
#>  $ MO    : chr  "0" "0" "1002" "1002" ...
#>  $ SEX   : chr  "1" "2" "1" "2" ...
```
:::

`loadExamplesPhen` function returns the GAW10 data sets available in
`SOLAR`.

::: {#cb13 .sourceCode}
``` {.sourceCode .r}
#gawdat <- loadExamplesPhen()
#str(gawdat)
```
:::

Both `loadMulticPhen` and `loadExamplesPhen` functions are based on a
more general function `readPhen`, which reads the phenotype data from
text files of a table format. The data from the R package `multic` can
be also loaded by means of `readPhen` function.

::: {#cb14 .sourceCode}
``` {.sourceCode .r}
dat.dir <- system.file("extdata", "solarOutput", package = "solareclipser")

mdat <- readPhen(phen.file = file.path(dat.dir, "simulated.phen"), sep.phen = ",", 
  ped.file = file.path(dat.dir, "simulated.ped"), sep.ped = ",")
```
:::
:::

::: {#data-exploration .section .level2}
Data exploration
----------------

::: {#plot-kinship-matrix .section .level3}
### Plot kinship matrix

The two times kinship coefficients describe the relatedness among
individuals in the family-based sample. Two functions `solarKinship2`
and `plotKinship2` allows to compute and plot the kinship matirx of a
given data set. `solarKinship2` function runs `SOLAR` for estimation of
the kinship coefficients based on the ID fields in the data. The
returned value of the function is two times kinship coefficients in a
matrix form. `plotKinship2` plots the two times kinship matrix by using
the R package `Matrix`.

For example, the kinship matrix for `dat50` data sets has the following
image.

::: {#cb15 .sourceCode}
``` {.sourceCode .r}
plotKinship2(2*kin)
```
:::

![](data:image/png;base64,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)

One can see that both unrelated and related individuals appear in
`dat50` data set.

The relationships among individuals in `dat30` data set is examined in
the following two lines of code.

::: {#cb16 .sourceCode}
``` {.sourceCode .r}
kin2.dat30 <- solarKinship2(dat30)
plotKinship2(kin2.dat30)
```
:::

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAASAAAAEgCAMAAAAjXV6yAAACu1BMVEUAAAABAQECAgIDAwMEBAQFBQUGBgYHBwcICAgJCQkKCgoLCwsMDAwNDQ0ODg4PDw8QEBARERETExMUFBQVFRUWFhYYGBgZGRkaGhobGxscHBwdHR0eHh4fHx8hISEiIiIkJCQlJSUmJiYoKCgpKSkqKiorKyssLCwuLi4vLy8wMDAxMTEyMjIzMzM0NDQ1NTU2NjY3Nzc4ODg5OTk6Ojo7Ozs9PT0+Pj4/Pz9AQEBBQUFCQkJDQ0NERERFRUVGRkZHR0dISEhKSkpLS0tMTExNTU1OTk5PT09QUFBRUVFSUlJTU1NUVFRVVVVWVlZXV1dZWVlaWlpbW1tcXFxeXl5fX19gYGBhYWFiYmJjY2NkZGRlZWVmZmZnZ2doaGhpaWlqampra2tsbGxubm5vb29wcHBxcXFzc3N0dHR1dXV2dnZ3d3d4eHh5eXl6enp7e3t8fHx9fX1+fn5/f3+AgICCgoKDg4OEhISGhoaHh4eIiIiJiYmKioqLi4uMjIyOjo6Pj4+QkJCRkZGSkpKTk5OUlJSVlZWWlpaXl5eYmJiZmZmbm5udnZ2enp6fn5+goKChoaGioqKjo6OkpKSlpaWmpqanp6eoqKipqamqqqqrq6usrKytra2vr6+wsLCxsbGysrKzs7O0tLS1tbW2tra3t7e4uLi5ubm6urq7u7u9vb2+vr6/v7/AwMDBwcHCwsLDw8PExMTFxcXGxsbHx8fIyMjJycnKysrLy8vMzMzNzc3Pz8/Q0NDR0dHS0tLT09PU1NTV1dXW1tbX19fY2NjZ2dna2trb29vc3Nzd3d3f39/g4ODh4eHi4uLj4+Pk5OTl5eXm5ubn5+fo6Ojp6ens7Ozt7e3u7u7v7+/w8PDx8fHy8vLz8/P09PT19fX29vb39/f4+Pj5+fn6+vr7+/v8/Pz9/f3+/v7///8Q9npJAAAACXBIWXMAAA7DAAAOwwHHb6hkAAANFElEQVR4nO2d+2MUxR3AvwkJXB4QMBEhiAUp4SGRmAQCoaAtyEOk2mLFFkul9d2ilAqpqAilpGgtoggooUBCq6AGqW2pqRGoRoqNDySIaAJEIMn8Gd3du+Que3v7nd2d3Z1lv58f9o7br7MzH2e/czvZuQVGmAJ+V0B2SBACCUIgQQgkCIEEIZAgBBKEQIIQSBACCUIgQQgkCIEEIZAgBBKEQIIQSBACCUIgQQgkCIEEIZAgBBKEQIIQSBACCUIgQQgkCIEEIZAgBBKEQIIQSBACCUIgQQgkCEEGQQc6+WPr+UMvvm29KsnIIGjI59yh7Vn8xR4dY6MuSUgh6AR36HkLgo6QIHNIEAIJQiBBCCQIQW5B7xTz03ccd+iEPvzFFmUbftxgrSFuCarJbVi7dm0DF2/whYmJrdxprSGuCaqohdpalwp3wi3yCFIMuVS4E2QSVCvDCKBHKkHUg1JTMz9aunSa5BLEgASlICpITdQkyJBuQdKNZfIJkmwsk0tQ7BgyJSIZBUk1lMkoSKpMLaUgmS7L5BQk0VAmo6DYgeQ4z6QVxCRJRNIKkiURSSyIelAiRoKkyNQyC5LiokNaQbGj+X6eyS1IgkwtuyDfE5H0gvxORPILoh6kYSwodkg/NQVAkL8TaAEQ5O9MfiAE+TnYB0OQj2NZAAR1H9if8ywwgvyaHwqMIL++EgVIkD+ZOkCC/LksC5IgX4aywAiKHd7zRBQsQT4kIk8FLQCFFva/6XlTj+t28QryPBF5Kmj89ubm5k5WvOJ01TjdLm5Bl3UPyv1M3TYMusQ6C/7eexe/II8ztZeCTqbdmFNUw7ZNVt5P3Rz7cO+zGj+dZqEKXvajyiXRCu7ljHci6HDZvraayLvrZyrv5zwV+3D13RozbuAvx9P5oZIZ0Qqu5ox33MHnrNpSobxUPtf7Y85hXsXb6Q8vT7FDajedv+ZQfifrGvxW733WBHk4lnkp6EC/vefrBn7EitZ3Vl/b1XufRUHepWpPR7Gd43JL3mTsoyn5k/+r22VBUAyPElGwvkkn4FWmDqwgr2bygyvIo6EswIK8ydSBFaThwUx+sAV5MD8UbEEeXN0HXZDrmTrogly/6Ai8ILeHsmALiuJqpr4MBLmbiIQKerfLbK8pjgS5mYiECoIrF75w0l49nAkKSg/a/VBZXxj30P526/VwKMi9TC06B52vr5oEFn7VqBsnSVrFtZl8wYJO7bx/YkZmqfV6OBXk2vyQUEF3XQvZ037z+jkb9XAqyLWZfLFJGr6z+4y9eggQ5M5YJlTQW6tnF6SNuXuzfsKZAxGCXEnVwr8oHts0BWxU1HEOiiI+EYkV1HW0+rbBULjIej3ECHIhUwsVNO8K6H/z+g/s1EOMIBdm8oUKKl9x8JLNeggSJH4oE52DLnxs42s0EydIeKYWK2jX+HRIG7vLRj0EJWkVsTP5Yq/FYOYL+7fMtvN9RKAgsfNDQgWV3KW9/MSPS404Yq/uhQrKrdNeXu1vvR5CBYnM1EIFjYrehvXUaOv1ENuDBF50CBW0MrK+patlQ2SV9XoIPsXEDWVCBXUuzYB0yLinw3o9RCZpFWFDmeDvQZ/ueWbPJ3bqIVqQsEQkUFDbay9/rFyOnfhP/Uzr9RAvSFAiEifo/WEA8LP3rla2acl7T31mXq4LgsQkInGCbr7ixcPbCoaWbd936Gzy3sXrlE1skYb9tRr8CPsyJE5QwRPK5gloNtq3Y36aKii2SMP+Wg2LiJjJFycIXlY2O433120Yu65nkYaDtRpWETA/JFBQjbLZlWr/XEVQbJFG77UaG3+pMbfcWj24EPGbDeXzohXcyBnvRFBskUbvtRpbV2ssrOCsgBVEzORPXhit4FbOeCeCYos0nK7V4EfETL7AUywzEolkQkQlea8qKLZIw+laDas4S0TiBC1PIHmvKqh7kYbDtRqWkUSQOZqg2CINMWs1uHGWqS+DG6gwnF2WhUKQk6EsBII0bE9/hEWQ7Zn8sAiynYjCI4h6kCm2M3WIBNm76AiLIA07Q1moBNnJ1CETZD0RhU2Q5UQUOkHUgzAszuSHT5DFCbTwCbI4kx9GQZbGsnAKsjCWhU9QFO5EFFZB3ENZWAVxZ+rQCuK9LAuvIM6hLKyCNHimP0ItiGcmP9SCeBJRyAVRDzKFJ1OHXRB60RFqQRrI/BAJQqY/SBDymw0kCJkfIkHI/JBXgrxeq2GZVJnaK0HaWg3Hz9Vwj5QTaN4Iiq3VcPxcDfdIOZR5Iyi6VsP5czXcI+VQ5uVtwIbP1Vi2QKOs2GbBgkiZqSeURyu4jLMgJ4IMn6tx8BWNB2bYLFgkRvND0x+IVvAgdxn20G4kV3H+XA0XMUhEXp5igp6r4R5G0x9eChL0XA33MMrUnq7VEPpcDRcwuuigS40EjIYyEpRE70xNgvToEhEJ0qNLRCRIj26sJ0F6dJmaBBkSn8knQcYACTIlPpNPggyJz+STIEPiM/kkyAzFUgUJSg2QIHPURE2CTFATtTSCpjVw8wZ/qLPYtQq5kgj6VzE/fcdxh07ow19sUbbhx+9Ya4gXD6bEGHKCO/S8hYfEHBljoy5JkCAEEoRAghBIEAIJQpBC0Ofcoe0WBB29bAQd6OSPrecPvfi29aokI4MgqSFBCCQIgQQhkCAEEoQQGEHYnet2YzF8FmRyK7oO7CkDSbH8RZvisyCTW9F7gT9lICmWt2gEnwWZ3IreC/wpA/pY7qIR/BVkeCu6MSmfMpAi1kLRpvgryPBWdGNSPmUgRayFok2RYBRLcSu6DuwpA/pYC0Wb4q8gs1vRdWBPGUiK5S/aFH8Fmd2KrgN7yoA+1kLRpvh8ipnciq4De8pAUix/0aZIkIPkhgQhkCAEEoRAghBIEAIJQiBBCCQIgQQhkCAEEoRAghBIEAIJQiBBCFIJ2n3T0AHFf+xI/Ch/g1+ViSGToAdg/oaN30/7ReJnJCjOQXhWfXkaDid8aCrokrsV0pBIUHmp9tJWuZWxzqrR/W9Ql1CqgiJblDd3zmKscG1Z9rf+8PmsgcNfUv5RPR3yFpx2u1byCOro93T8H0sjj9ctgt16QZkrXp+VVrj2LyWRNlY46If1T2Ysdbta8gg6Bn/uef9pxnplO/M6vaAFjDXBvYzVwRFWOL6LsZsnuV0teQR9kCBoD7Qo283pF3SCnlA6GmxVLTWyQvVn2h4sdbta8gi61HdN9M2GavZMuvqnvv3wsU7QGlVQTUyQ+gf3MAli10/UXr7Ju0/pQaeUd1vS2uOC5iYJUn2GStA+0O45WAYH2CcZ6ug+e2w0B2U/xti5K0kQWwKzf189F5Yob+/JWvPqYtgVFVR61eb6m9JJEGPbpxQMKt2kpp+Old/OLdnDooLer8iC639OgkzoOuXboYMhyEdIEAIJQiBBCCQIgQQhkCAEEoRAghBIEAIJQiBBCCQIgQQhkCAER4JKASBt6NQd6vuB6vwfN5aiG7dfYuxBiLJOXCwXjgWpPMbcFHRmJJyNN1r7U3RlJPKh41g+HAq6o+Pi0e9C5ARjHR1WFmXxR9c8NhzURnd2KLyaPvxr5bMXlcY3OYzlxaGgO5Xt2QGwXOsTZwGeuTan7N/7i3NL/snYhUdGZ1+nLqcF2HtbfuEjnYydWDQsa9Sq9lgPqpmUN/zW48YBX9bUfKkdI1/tCmejx/u6MO1N9WXwDd2Nfh7gDfYkZDRyxNpDgCB2K9zSLQgyAfIylc2QDvY9iIwBWK22vyC9D8BzrGsc5I1Ohx9HBT2pROVATpNhQANAg3aMI40v9TT6flisvUQO9TR6GpS15sPDjCfWFiIE3QtjuwX9oHUHwO2t2wCO7YesFrYJIqeV9ld+cnIU3M6aAI6z1wZf06FGt2TDKtZaCnOYUUCPIM1VtNHNmZnNysvRjJVf9DT6WASmwzXnGE+sLQQLqmdtaqc/A9C4AkZUVf06DfYq7d+l/h+dzlojcNWSXa1Mi66DgcqAsxsGMaOABHoa/bDatRibOvKbhEavUrrtXzlj7SDmFJvXLaiRaZuvlM2i2FDyvNL+vzG2XGk/+9MI5YOcp7ToalB/hO0IQKtRQALdje66WlHP2DZ4jSU0+qu+UMQbawcBgs7lwaMGgh6FO7qP0dN+xg5XKXmpWY2u1XpQLeQZBiTQ3ei3tf+A3RkTvz2691fK2xrOWDs4FtT14SyIfGYgqBYKv2DvTixpibd/06g5HawlHf6hRp/MgsdZWxnMYkYBPaNYvNG/VUOVk2fEiBHXAFxdp+18LwPmweAzjCfWFiK/KPYWxCogv7wf3JfQ/qYcGFA8AIouaKPYaoDhAyC7yTDAIEkvUG+ditJz2nSWwIz2odEBC4u1h3NBQ6Z0X2roBJ27b2T2ePW+zPgZdGjO0H7DfvRp7HvQK+V5w+YfNw4wEHRd/NKhp9G/A3iHPQtp9Ryx9qCLVQQShECCEEgQAglCIEEIJAiBBCH8Hx9GPyI0748hAAAAAElFTkSuQmCC)

The structure of patterns in the kinship matrix can be explored more
closely if the first 30 individuals are plotted.

::: {#cb17 .sourceCode}
``` {.sourceCode .r}
plotKinship2(kin2.dat30[1:30, 1:30])
```
:::

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAASAAAAEgCAMAAAAjXV6yAAAC+lBMVEUAAAABAQECAgIDAwMEBAQFBQUGBgYHBwcICAgJCQkKCgoLCwsMDAwNDQ0ODg4PDw8QEBARERESEhITExMUFBQVFRUWFhYXFxcYGBgZGRkaGhobGxscHBwdHR0eHh4fHx8gICAhISEiIiIjIyMkJCQlJSUmJiYnJycoKCgpKSkqKiorKyssLCwtLS0uLi4vLy8wMDAxMTEyMjIzMzM0NDQ1NTU2NjY3Nzc4ODg5OTk6Ojo7Ozs8PDw9PT0+Pj4/Pz9AQEBBQUFCQkJDQ0NERERFRUVGRkZHR0dISEhJSUlKSkpLS0tMTExNTU1OTk5PT09QUFBRUVFSUlJTU1NUVFRVVVVWVlZXV1dYWFhZWVlaWlpbW1tcXFxdXV1eXl5fX19gYGBhYWFiYmJjY2NkZGRlZWVmZmZnZ2doaGhpaWlqampra2tsbGxubm5vb29wcHBxcXFycnJzc3N0dHR1dXV2dnZ3d3d4eHh5eXl6enp7e3t8fHx9fX1+fn5/f3+AgICCgoKDg4OEhISFhYWGhoaHh4eIiIiJiYmKioqLi4uMjIyNjY2Ojo6Pj4+QkJCRkZGSkpKTk5OUlJSVlZWWlpaXl5eYmJiZmZmampqbm5ucnJydnZ2enp6fn5+goKChoaGioqKjo6OkpKSlpaWmpqanp6eoqKipqamqqqqrq6usrKytra2urq6vr6+wsLCxsbGysrKzs7O0tLS1tbW2tra3t7e4uLi5ubm6urq7u7u8vLy9vb2+vr6/v7/AwMDBwcHCwsLDw8PExMTFxcXGxsbHx8fIyMjJycnKysrLy8vMzMzNzc3Ozs7Pz8/Q0NDR0dHS0tLT09PU1NTV1dXW1tbX19fY2NjZ2dna2trb29vc3Nzd3d3e3t7f39/g4ODh4eHi4uLj4+Pk5OTl5eXm5ubn5+fo6Ojp6enq6urr6+vs7Ozt7e3u7u7v7+/w8PDx8fHy8vLz8/P09PT19fX29vb39/f4+Pj5+fn6+vr7+/v8/Pz9/f3+/v7///9090knAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAUFUlEQVR4nO2dfXwURZ6Ha97f8koyIa8CISQhrxMG8kISCBoBBY1IcrucdwuymNNFV7PnIYq4DueCeHDLgqBxYd1FQYEVSNxFoizyoqJxxReU5cgu+HLEWwUugJwQUp/P/ao6r0NXV2Z6ZjJh6vnjS01Xd0/PQ1dPZ7qrGmGBImigNyDYEYI4CEEchCAOQhAHIYiDEMRBCOIgBHEQgjgIQRyEIA5CEAchiIMQxEEI4iAEcRCCOAhBHIQgDkIQByGIgxDEQQjiIARxEII4CEEchCAOQhAHIYiDEMRBCOIgBHEQgjgIQRyEIA5CEIdACtrbwag4+VfmIqyKz1oZFZcOshZpbmPVKBFIQTHfMiqWLGJUXLCx1lX7LKPi0yzWIlNeY9UoEVBB3zAqXCxB59mCnmFUHGEL2sWqUUII4iAEcRCCOAhBHIQgDn4RlGklZDudziRSMJogTEarxmo2ktcJTjeSEt2ndFKgY1Q47cMZFdkW1iKR6c5azz+LXwSZEeHp5ubmm26uq6uLq4AYXgiRkwMxo6zZjYMH3Kd0sYdVse9tjxfZ8+6mPM8/ix8FNUHpjp81NDSMnAcxphqi8kaIxdP88Z794UMhSBkhiEOQCRpTWVkZaQ4PD9eaIHRGCL0BwsL8nvE3QSMo/dmmpqYx1FKNy+UaehNE6gSIfCfEP0/yx3v2h6ARlPsRRCUVJJqYHEIQByGIQ/zkmpqaCCooPDY2VmeDMFggjGaICC820zcEjSDzjAULFph+CBFGLd0ApeuKIbLyIKrL/fGe/SF4BC2ElhTughhKBYkm5oYQxEEI4mBIzMjI0CZD6KkgO5QsQyDCIyGSnf54z/4QPIImzZ492zAZwjIBIopayoXSiDSIG4r98Z79IWgE9Wpi9DyIChJNrAchiIMQ1Isa8uHdLpQbRjqdTu0oCMMICOnXj3goRcdApBWpfU9vGRBBjs0nTpy40neaoaS6ulpfDmEqhAjLg4ihloZDqTS0vsXCvrx6mnsTo+dBVFDoNbFTmsm27K2dL/ZES2iCVNBefecG7un/MmoFfVSy+9xW8wfSi47TEqYgFfRhducGsm5UksEn32JVrr6vjTkVFRX6fAjjaAhbGoR0tpgApbxSX7ynNwxEEzv0KsTMp/pO1BdMnTpVNw7CmAdhzYSITIWIo5byVb6n1wyEoDdNr37XGNXSd6LsQZqeB1FBIXWQxtvywwr/5DZNCOIgBHEwz3K5XJY5EFHkulgiuS42ilwXk77KipqAT2C+jj+RUj0JciWt6bndEK9f8scmUYJGkL24srIyvRQiqwgi1wkxJh+iMBuigFoqhPm+MDocjlwNhAPlQmhGQ9gO+GOTKEEjiF72YdJIBY2F0kk7NLuXDBANms0Q+nUQ2fv8sUkUIYiDEMRhUAnKam5uboyoq6u7VwdRh+6F0M2FSGJ2P1DNIBH08Vin0zmaWtKnpaWlaiHSUCqJ6yDMW5UWVsUgEUR5V7rcIZoYCyGIgxDEpuNvLS0tr1BBUU888cRjeognNI9B6Ooghoe8oCZTfHx8hAFiiHSkhlI8iiMRC2F82R+bRBkkghrHQUO6fxTEcipINDE3QkfQ1zIXNChCkMRd/wnxt8rIiuNuFUIQYctMDRHkXPztEvcfmZUFraVaNCR06Go018pBunF1Hghqjr6Mr9jf6lvV7z2IngdJXq69PQjj20DQpnIoVDwvTTixVCLhXaXFBk7QoYTODTzR/2XUClpFfmCuWi5NaFkgMfSQ0mIDJ+idOGn7Hmnhz9uFWkEbJ0Bh0q/7VogmJkEEHYq9gjvi3T6TsqCnkQboFdKxufvltXKQlgTh7FVX1qS7Xevm7EEFmzZtuicV4vFYiPUGiE2SphVQyry29iDcMjG2/L/cKjxrYtJ5EBV0jTUxJkIQByGIA+dMWqPT6bQkdIiEVKLQl9fOQZqJsqCGvPr6+rnDIRYOgfiVHqJeswpC0hTye5B8E6PnQUIQQQgKKUEfeHC/YzfKgv4Y4XA4hlkh0undHdqeuzuooHGVALktdD0p5OVBOLMgitLJ60e92J4u/CEIxf3Tb095ulJlQZf2NDU17X6B3BC00S3I/UG51BK5vHrPLS6Xa/gkiNxCiLJsiHllnm5ML/whaMeDJUaU/2DTRU9W2o/rYmzKegTdA80uZxbEhJshbh8PsSzYBAHf7V1SiiyerDS0BH297YFxeoNHneBCSdDcdGS9/rE3Lni0UlWCyqmg8TU1NUPDYmNj9VbS6Z70tzebICJHqFi1Xw7S6IYdZzxdqSpB778MFFFLkxcsWJBUBpFeADE2A+IOh4pV+0PQvqW32jU5tc+7/6KhiCpBlJlU0GBoYoRjGyYiufm8vHDYHwaRoI5P1vwgHiXfKVNFLxzK9TgMJUEzYlD4Las+k6npvHAo1+NQlaD6WkDq9RKdkpJiiIAwh0HYLBBxaSpW7Q9B4xfvvyxfI104lO1xqO5r/vb58+dbpkJESn95QGlEDkTeSIjbc1Ws2k/HoO9PMk6jyW/SfXoc4jNSf7VsVYKWQUuKIb09rqOCfNnEcqTt+18PluEK2u7QIk3edrkqIqhPj8OuLpn6Nz3YAnf8KMgfXTJ3oGm/bdp4K2qQqaNXNQjuPQ7VNTE/7kF+aGKFc+k/8+T+1KAXDuV6HKoSlJXhdDr1pOe9iQoi47xG2iHioiEy1Iwx6A9BYY30n13hMnVEkGyPQ1WCMiuqq6uNxRC2Agipr2sqlFKTICoyVazaH4Iyl9J/lsv9x9EmJtfj0FdNjJ4HUUHB28QeN69q7WhdbXZx5utDKAm6Ml+PtEj/k3ZPVhpKgjD+Yue6nZ97tlJVgnILSH/7bAhLOoQ02mAilJLtEAXBdaJ47rWXTsKfY18d2etRR1xVgtKKpk6danBAWLIgItIgYqmloVAqCqo/NT5NIX8xfkjORzSerNSnTYyeB0k/ogVdE7sl5ncfbbInlWzefei8JysNGUH2ZRDLkAe3PEqEjCD0EsQ2zy++qhLkKK+qqjIVQoTlQwzJhEiUjtRQKg+qn1zptYXtARY0Iq+srEyfCWEaCWFLgYhKgIinloLqR/sBESTbxOh5EBUUXE1MCOIIMpjNZgMyEzxZacgIWtQLT1aqTtA8l8sV8UOIOLebF6igXHKXA/2N/BApvdAdL5I7H/YojowySO4P4vAoucklnQyKkuW8+vaXcdQS+c3jnI7eNpNHbpvJgtCmQ0T8QWnV14YgZfZTQaOgdNZGbrxCv4fQrocwr4AYJ/fLZzdCkBA0mAQdLg8f9TtvumSqQhI0oqWl5bCZ9LdHLgjtAgjjfIjMnUoLB1TQ2agV37xgO+RFl0xVHEtLTU2VLnfoaH/7oSTsJGIgTGuUFg6ooAbyw0zVUi+6ZKrnL1RQkDexti8xPp+8p0+XzCO1EjFvKS2pHm8FvRXTuYFH+v9eag7Sh8fMa+/TJbP1GYnk91WstR94K6g5uXMDWQ+KlsF7QRfuS/y1V10y1XH29OnT71FB1ucAtBZCuwzCtAgiP3ia2JWKmeQeAC+6ZKriz7qwsDALCRu1hKAUhmwkLBD61UoLB1TQzsQvW1tbz3vRJVMV+7OgIT2eBLFeEhS0B+lH6eYt9aJLpioGjyAlhCAOQhAHPwragpQJnoO0EoHZg6TzIIrYg7oRgjgIQRyEIA5bNN1d7d3620ud7oPn5w4l/LkHpdXX1z8wFGKZhfS3R2shtNTSI1DKC54fzJQIZBOj50FUkGhiBCGIgxDEYZvWZDIZNBAm1BMUI5S0a5UWDgVBb6auXLnyHjvEY2aIlWg5hHYxhIFaku1U0kUoCJJvYvQ8iAoK+SY2iAR1Xjj0T5dMJoNHUNeFQ593yVTmoM3hcIw0QWS7391BBTkqKytvJAOBbic3iThzyV0h5NaQUnJrSHGc52+o9sKh77tkKnPpdTIySs/zkqQgtwY9twtC6m9PLrI8dKPL5Ro9HqIoD2LySIhZHt0GJqH2wmGfLpmXWiQy/Hxllc20HkGzodkVV0FMq4CY44B4yNS5gR48gErthcM+XTL3kQvngElNl0xVcARppO1L82AMJ7UXDik+7ZKpCo6gQDaxrguHPu+SqQpJ0MTa2tqRkSkpKaZwCKsVItwMEWvwfI1qLxz6vkumKt4hV95LqaWS+fPnX+eAGJ0OUZQCMdXo+RpVXzj0eZdM9cylmzbQTUwJIYiDEMRh4AS9QkZkT6WC7BkZGdZoiIgIiBgrRGIgD9JKDOC32A2zZ8+2lUJEU0tZUBqZCuFIgJgkBE1bDC0p6T6ITCpINDE3hCAOQhCH4jSn02kcBmGhguKgFDMEIiECIjWQJ4pKDJwgZ2l1dbXFCRGRA2GnllKglBkHUSIO0r2aGD0PooJEE+tBCOIgBHEoL66qqrIWQERmQQylgpKhNNoOMU4cpPOyysrKjGkQluEQkUkQ0pE6CkqZes/XeI0Jcm9i9DxIutwhmhhBCOIgBHGYVDlnzpywUogYJ0RyNkQ6FZQIpUkmz9d4jQnKSnU4HPoUCFMChM0OERUNIY2eKw7Ssk2MngdRQQPdxN6tlIjY78u1egJHkFbavpve6/8afSqorUki9bAv1+oJHEHGzg1s6/8a/dPENjfLsUd2KnDwAKuGuci+t2Unl82oq6uLngyRWA6RNgbCMRpiQlA0sS7mOWUxOOSnO5MSGRUFOkaF0z5cdnKC1aK1Gs1Wq9VogjAZIcwGCIvBqjNbvRgI1i+CGMR8w6hwsQbfOW9jrav2GUbFEeZwwVN2sWqUEII4CEEchCAOQhAHIYhDQAV9y6hYwhJ0gS3oWUbFp2xBr7FqlAikoL2s5yWe/CtzEVbFZ6zxSS4dZC3S7MH5cw+BFDQoEYI4CEEchCAOQhAHIYjDwAtiPgiQWRFQAidIrmMZgT4IUG64StYTApkjXCr0YPOewAmS61jW/SDAq4erZD4hkDnCpVIPNu8JnCC5jmVdDwKUGa6S+YRA5giXSj3YvCdggvo+668X5BlTfYar7F0hs5TcCJe9K5hv5B0BE9SnY1lviIc+w1X2rpBf6qoRLntXMN/IOwL7LebesYxAPMgOV8l6QiBzhEuFHmzeEzBBsh3LCPRBgHLDVTKeEMgc4VKpB5v3BEyQbMcyAt1R5IarZDwhkDnCpVIPNu8JXBOT61hGoILkhqtkPCGQOcKlUg827xn4M+kgRwjiIARxEII4CEEchCAOQhAHIYiDEMRBCOIgBHEQgjgIQRyEIA5CEAchiEPQCNoxJSnCWd/ee1Ks4jMgAkSwCPoZmrl67T9o7us9TQjqYT+iN2WuQL07UikKuuzfDeomSASNL6b/nJv0IsZXlmSFF5ExWYkg80YozJmOcfLKEmvq0/89PWrYC/BiTSWKrGHdE+pTgkNQu2lFz4v55l803ol2uAsyLH5juiZ55R8KzedwcvQde5/Uzw/EpgWHoGPole7yF/pVkNMK3AXVYHwU/RTjRvQxTnZ0YHxLaSA2LTgEfdZL0E5668rz2u/dBC2DHQ29SCwdxskLYeq/Fgdi04JD0GVj55XQ1WvwOi25GtiETroJeooI2topiFyTDyVBeOw4+s//Rd4Pe9DXUNqoudgj6LarBBGfISVoN6J3KixEb+LP6YMcb82TjkHWn2N8IU4IwnejW3+15jZ0NxR/Ynlq113k8QZEUHHC83unaIUgjDdPtEcXbyCHn/bHR4cVkmf0EEGfTrCgsfcKQQp0fD2gbx/8ggYYIYiDEMRBCOIgBHEQgjgIQRyEIA5CEAchiIMQxEEI4iAEcRCCOAhBHHwjiAxsrUmq2ELKUcrPOnWj/3P//cfxtoINpLRqTNSUN7ycxWN8J4jwc+w3QZfzkdmO0HOkV48mEZmavJrFc3wl6Eftlz6ZisxfYdzezhoFR45+z70b2U+116BS3BaGXscLUIVXs3iOrwTNgTwfgRbRfeI8QuvSbSXvNznDCt/B+PuHs6wFpPstQq/+IDb54SsYf3VniiXTdbFzD9paGjms+rj8DKe3bj1N3uL3Of+C8VJUhF9CIzE+jrRnydTfILQHP4n0h9mzqMSXgnA1ur1LEDIgFGmASGzHNyNzDunpBp/frtWRRxB25KPILC36sSToSZjLhmxHZWdoRqi5620+H43W4X9HZRhfQOgYnXQ9KmmLRf+GFWZRh08F/RTldQn6x7YtCM1q2wRb2YQsrXgDMn8Ln3/S56cy0Sx8FKHj+LX4Ee1k7lYrcuG2YlSF5WboJehhhG7uwA+hKVDWorfptGNmVIlGXMAKs6jDX4L24nNk7z+D0OHFKG3Jkkc16FVMH+b9AKrEbWaUcPd2MmQWzN2Ioi5jvANFY7kZetE4x4hq4QBMPr0OdQ7F5YK99Y/Ks6jCx01sRpegw5jGWYg7O7/ifgOf/wDGi+Dz4/VpMMG2nM69BuXAsh8j1CY3Qx8akPH8L9F4jC+SPYxy1oiyObOowpeCLkSiR2QEPYJ+1PVm3Z8f44+WwHHpBJm7ge5BDShSdoYuNk+HA83/IPTnl1ECxkcQOidNfwhEblWeRRW+E9Txl+nI/KWMoAaU/A3+YFxha8/n35BZ1Y5byUEC5j5lQb/A50rQdCw3Q9e32BakP/Ddg0h3/owRbbtUi66X3vlDPZqB4s9ghVnU4ZcTxb6C8AQUO96E7u/1+Y/aUIQzAmV/T7/FliI0LAJZj8rO0HWQvpRHvxkXYjjBQSakl0aNvFKIbryYhO7C7FlU4kNBiRO7/tRwE3Th/lFWB7lJs6cFHapKMqXM/qLzPOjl8ZEpM4/Lz9D9Lfb3u4bbHGvJsDf/4YyY8Lr0xr9E6D38DNLsZc+iEvHHKgchiIMQxEEI4iAEcRCCOAhBHIQgDv8PCI+6+e/pp68AAAAASUVORK5CYII=)

It seems that the simulated data include only sib-pair relationships.
Such an observation can be visually confirmed by `histKinship2`
function.

::: {#cb18 .sourceCode}
``` {.sourceCode .r}
histKinship2(kin2.dat30)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```
:::

![](data:image/png;base64,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)

One see that all the related individuals have the twice kinship
coefficient equal to `0.5`. The unique values of `kin2.dat30` matrix can
also be extracted by the following code.

::: {#cb19 .sourceCode}
``` {.sourceCode .r}
unique(as.vector(kin2.dat30))
#> [1] 1.0 0.0 0.5
```
:::
:::

::: {#plot-pedigrees .section .level3}
### Plot pedigrees

The family trees are plotted by `plotPed`, which reuses some routines
from the R package `kinship2`. The family with index `2` from `dat30`
data set is plotted in the following code.

::: {#cb20 .sourceCode}
``` {.sourceCode .r}
plotPed(dat30, 2)
```
:::

![](data:image/png;base64,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)
:::

::: {#transform-traits .section .level3}
### Transform traits

To show an example of trait transformation, a trait with an anomalous
distribution `exp_trait1` is created in `dat30` data set.

::: {#cb21 .sourceCode}
``` {.sourceCode .r}
dat30 <- mutate(dat30,
  exp_trait1 = exp(trait1))
```
:::

The list of available transforms in `solareclipser` is returned by
`availableTransforms` function.

::: {#cb22 .sourceCode}
``` {.sourceCode .r}
availableTransforms()
#> [1] "none"    "inormal" "log"     "out"
```
:::

The following code plots the distribution of the original values of
`exp_trait1` and its transformed values.

::: {#cb23 .sourceCode}
``` {.sourceCode .r}
library(ggplot2)
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:solareclipser':
#> 
#>     annotate
library(gridExtra)

theme_set(theme_bw())

with(dat30, {
  grid.arrange(
    qplot(exp_trait1),
    qplot(transformTrait(exp_trait1, "out")),
    qplot(transformTrait(exp_trait1, "log")),
    qplot(transformTrait(exp_trait1, "inormal"))
  )
})
#> Warning: `qplot()` was deprecated in ggplot2 3.4.0.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> Warning: Removed 18 rows containing non-finite values (`stat_bin()`).
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```
:::

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAASAAAAEgCAMAAAAjXV6yAAAC2VBMVEUAAAABAQECAgIDAwMEBAQFBQUGBgYHBwcICAgJCQkKCgoLCwsMDAwNDQ0ODg4PDw8QEBARERESEhITExMUFBQVFRUWFhYXFxcYGBgZGRkaGhobGxscHBwdHR0eHh4fHx8gICAiIiIjIyMkJCQmJiYoKCgpKSkqKiorKystLS0uLi4wMDAxMTEyMjIzMzM0NDQ1NTU2NjY3Nzc4ODg5OTk8PDw+Pj4/Pz9AQEBBQUFCQkJERERFRUVGRkZHR0dISEhJSUlLS0tMTExNTU1OTk5PT09QUFBRUVFSUlJTU1NUVFRVVVVWVlZXV1dYWFhZWVlaWlpbW1tcXFxdXV1eXl5fX19gYGBhYWFiYmJjY2NkZGRlZWVmZmZnZ2doaGhpaWlqampra2tsbGxtbW1ubm5vb29wcHBxcXFycnJzc3N0dHR1dXV2dnZ3d3d4eHh5eXl6enp7e3t8fHx9fX1+fn5/f3+AgICBgYGCgoKDg4OFhYWGhoaHh4eIiIiJiYmKioqLi4uMjIyNjY2Ojo6Pj4+QkJCRkZGSkpKTk5OUlJSVlZWWlpaXl5eYmJiZmZmampqbm5ucnJydnZ2enp6fn5+goKChoaGioqKjo6OkpKSlpaWmpqanp6eoqKipqamqqqqrq6usrKytra2urq6vr6+wsLCysrKzs7O0tLS1tbW2tra3t7e4uLi5ubm6urq7u7u8vLy9vb2+vr6/v7/AwMDBwcHCwsLDw8PExMTFxcXGxsbHx8fIyMjJycnKysrLy8vMzMzNzc3Ozs7Pz8/Q0NDR0dHS0tLT09PU1NTV1dXW1tbX19fY2NjZ2dnb29vc3Nzd3d3e3t7f39/g4ODh4eHi4uLj4+Pk5OTl5eXm5ubn5+fo6Ojp6enq6urr6+vs7Ozt7e3u7u7v7+/w8PDx8fHy8vLz8/P09PT19fX29vb39/f4+Pj5+fn6+vr7+/v8/Pz9/f3+/v7///9khLolAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAYiklEQVR4nO2d+2MU1dnHR0mhUGmrpUCqVvFesVyUAkotLZdALLAQUQQJIBLlKkUaK8EWo8hFRSUkUnx9RVBRVJCQEJKwQRGseRHEGiHVwt4v2d0km33+gnfOmduZy+6Zze4yG3K+P+TMPDNnzjefndmdy3POcMCUUJzVBrJdDBBFDBBFDBBFDBBFDBBFnQNUvWLFiuVLlq8gtFQ1s4ScW6ZaMU611zttZMVSou3ly4h2SRfkOnHjhLVtKQGqOAYQdbSSIR85E/CQcy3txEyHI2JULfp4p42Aj2i7LUS061KmY35lOuLo0BsA9X8ULRInGCBlkgESdZEA2Ww2MmQhINJIlgAqqwsGeUBBQm5yxuNSzfmJmYDDZ1TNV0Rv1UDZCqj8aHs7D6idkJec8bnJuUCEmGl1tBhVizzWKSPZCogdYiZ8MUAUXwwQxRcDRPHFAFF8MUAUXwwQxVdWACqrDQQCNltAls+jTHucynTArUx6HX7DuN/hVdZ5VGyhiwMq/zQajdpsUVmRoDIddCnT7T5lOuRoU2a8ymSrI6xMLxJb6OKA2CFmwoiVgOra+D+NC2e+1iEXsi8GCKB5agu/6qwvwsurpELxxQDB3+8bzwM6tgzgcLFUKL4YIF75PKCPNgI0LZAKvtndu3c/cyQc5gGFCXnIGa+LnPO1EDMhR8CoWrBr3g9CgHa+BnC+QCoAfhgyZMjsgw4HD8iRPv330YRG4ikbAH24id915kmF4osdYiAAOvYEQH2xVCi+GCAQAEUf+Kbjr5VSofiyAlD4A9Cfb1gNCL5aVFgWkwvZlxWANs4F/fkGO5OWdaRoLujPNxggSe5FjXNBOd9oO3ny5Kaj7ej5ivKUJBwknqa4lOk2nzLd4mg1fBrT6gjJ05EueLEaKz7+3VzQn2840nq+IZ12LBRb7UKAPngVECD5fIPtQRqtnzZtyoRpXt35BvsOUoT2IN35BgOkCAHSnW8wQBQjDBDFCANEMcIAUYxkIyCWQEURehxlIx9HqZ4vRaN+NzkXaCVm2hwho2qtLIFKVPc4xBggE74YIIqvrABULlysKgml3fli1Ujo5zSIAYm/hn6v8gvpdRr/yvocAcM4+QPru0SyO8RDTPHCDjG9EQaIYoQBohhhgChGGCCKEcsBLRk/fvyW7E6gshZQgb+tLZrdCVSWAgoX4CKrE6gsBfRtwaPT17izIYHqyMLpqy8YJi9YCujUWlf0uWeyIIHq/NTT4c1/M0xesPxX7PTULEigqioBODfdMHnBUkCnTwN8U5AFCVThAED1auVYl4xYDah21vmODZuzI4Gq5oGvlGPdWVBQUHTI4/FgQB5BbpdHlsuhTHucZNxtGHcTcWcSyQu75zywviUbEqj8a4qaiOQFX0lJySqhM4tN7s+SdGcW16XTmaXtMfwBGSUvWP4lrZUlgA6tFNY0SF5ggJDK+Eue8TMMkxcYoERGGCCKEQaIYoQBohhhgChGGCCKEQaIYoQBohjJRkCoP7+N7M+vuqbhL3xUcx7i0kd1wUNcInlTSqDKOkBiAhWZQmVpAlXWARIPMfIgs+YQ21rj9XqxE68gj9sry+1Qpr0uMu4xjHscSmVXWq7mrQek2oPQX7YH6Y0wQBQjDBDFCAnIZmOAdEYEJ5IYIJ0RBohihAGiGDECJM3JKzJAhEAFSDAZ8ytuGSBrABkkUJGMLh4g4+wOY0Cyw8wDMkqgsgRQnOyOBEIrioDwHx0gcZ0ocUgkD8gogYpEJP2neD6TgOJkd1gOSE6qCFZUVKypD4U0LjwhQfxkKOR1hQj5gsRMi8NPLpOqBUzfD9IZ0TrRCTfrIQKOFsWrR/IcCgWJGoGkr+b1I0KkVeZHoMqwEclP0kNT6BKorHorgs5IlrwVQZdUYRUgfXZHdgBSJ1XsaGg4sr+ugVAVOVNdSc7V1BMz9v21RtXqzf+KaY00VBFt19UQ7R4gmiXWqd1vN/RN/kf1yZ8HkUkVX/Lfjy/Pfr6CUBk5s6pItaicmHl19jpy2VZpojJh44mMVCxdqmyvfJsyXTyfiBP+np39qoEBzX+0P2lAWp0bUhd32dPT4y5yDdkfd1kntWiRcfzFscbxT4a4DONnjf4jBkgRAyTqIgFqbw7HXeY+H3dRR3Mo7rJOyuk0jvu+N46HmjsM44b/UecBdRMxQBQlAYi8zaANvjV54kP/jbNM7G9mdlHSIn0ZjBDOF08vnL66duHMv2vjK2fNWPe5QZwcX9w8INVtBk3w+wl7PTO1X5VShYLnHm6LxlmEu6KlJtKXwQjhfHF24pvhTflfBCdvV8e/mrwnsHKqPq4aX9w8IOk2w/GFsze0aYLl/K/We5OMl4WnovEzjRcVpAiH9AXGI4TzRdXiYtiXB8cKi9Xx+pJi2Fygj6vGFzcPSLzN4LU1dbzwMgrUKME9zwHs+rPxsi/yHp64xm24SOyKlppUnVv0I4TzRfjkAtjwIHy0doE6Dk2FJx5abRQnxhc3D0i8zbBvLcDXM1CgCMjRuGvun2a4LLZ0xb/mPveMYTWxK1pqkm9/SIDUI4QLRc2Uf8DODQXa+NQFBS8YxYnxxc0DEm8z/M/MIl7SfyqNxv38mqJ642V4/MzTUw2rIZ2emiIg1e0P/QjhqPi/Pxft2AQfPjNPE+eLZ2YZx+WtJvEdJNxm+PAVgMgZOLRwYf7ChafEYMOUslicZeunTJkyYcoMw2piV7TUpLr9oR8hnC/aHp4d44tjC4rV8Xe2FMN7+fq4anzxJH7FhNsM52c2tb+0AQWKlGDVpPjLamcdn7thc5xFuCtaalLd/tCPEM4X1TZctE3+X3X80NTdwdL79HHV+OJJnAeJtxnq5hWs8Uv/qRiU+pkYLYPdD05c3xJnkdAVLTWRtz8MRgj/atE05G7aosJnNfHC5bML/nHMIE6OL87OpCligChigChigChigChigChigCjKOKALAfV8hDshxZ4IJFpRq8+HjYq7bBwnaJxhS62/dXTeZcYB/a5UPd++9D9irJZzJlpRq1l5F+IuO2W3v8Ftt9tPGbZUUpiCy4sFqEUb25/fyxBQ3BPrsU+pZrXrNfKfuXFLkQFfpeAydUC+wl//dNw5+LDPN9B++1/s/V7L7XvPSXnpUI4bC5z93gnw9cT+fUedQDsvjtWWFqkA6VaEnHdv633DLoA9g3tfvx5Gc9xd4CwYmDuDP17wejkvXtPn7uaiAb/El3gCIKOWdo5MxWXqgO4Zfejo9Fs9cP/vYe2tEfvluW/X/KmfkhKAPgduxI5zcMvdlQdHDsNHt/B5fabfg8gVIefqtxun/Tj8Xc8nPl3H1aM9KDZseFXV8GHiejmDDlf/qvfqU3N6eBRABi3NWZWKy5QB2X/E24sO2Aqu/qt/chTsHP+Zh/pvVDe9FiBW+m+A7f0ogIgVIYefOc2dOXDZGYAPmhCgqh7fATT1OITXg5x/AjxyM8C3XKMCyKClO95IxWXKgMq5HF7c0wA7uWW8E/xf582Tl+OmP+Ynwm+tyr+SBohYEXL2ATi5M6ERvfI2nsXfQS/diFYctEVYL6ca4Mk/AHhJQAYt5X6cisuUAb0zUJraxE2Sms6fo266FiAweOi69ytogIgVMQAeEMQOP3lHz3cRoBcxoBs34/UEQGM1gAxaunZvKi5TBnT6Mv67rvnOL+DMFW/23cU3vZv/GHKf1zW9p6cbYIcJQNKKEqBKdDhNnIQAHcw5B3A2pyoxIG1Ld1ak4jL1L+nJN31UOeq2aGxkIWwc6LFz1753ZOJVSobZyEc8uOk6buvZ3df1caCmUUwHSLuiBOjg5Zsad15Zgr+kh46oPTxiaCwxIG1Ljz2eisvUAbXMv+bneU2wIdcLHcMfsnN7b+87ulFZvPUXk4X/p2TgVfedumUYahrFdIC0K8qH2AuDel23sh2fBzmmDxiIf+YTAdK2tG9wKi7TfaJo51J9Upp2ddygy2pJwmU3AARl+dqIhYCOX4+brhkjakeat98ZxaZpL1aTcMlud1DEAFHEAFHEAFHEAFHEAFHEAFHEAFHEAFHEAFHEAFHEAFHUtQDhRHrx7RGQpkR9iroUIJxIL749Aikdifo0dSVAQiK9+PYIgPQk6tOUGBDepckduXqFqKUrNFq6RBfRBpYs0wSWE4HXzbhFOZri2yMAlER91/z589cuWbJk8RKsx4VysVg8LkQTLZTnFsurbjUDCO/Sqh0Zd8NG8sc06/o9moC6oziSK6gJRJWAuW7h+cLDc/T2CAAlUR8BWlLj9TqFMRCdwnCIbrGQoqqxEl3CUIpOMepSDbDIB11S15xEgIRdWrUjZwMg4e0RouREfeRMbNElYIwIHfSl3vPCQqnbflD4yH1eoXUh+6XVEZVWlf0kPsSQIXlHbjt58uSmo+J7SH1t7Wr53JoA+QZTQS6/JtCqBCKmhipFfsS3RwCoE/WtAyTvyGbGyrCJb0SyJftmJHNjdyA/4tsj4HRIlahPB6Qa8CSNgLAZtCOb2YNs4itzUZmhPUjK6s87qUrUtw6QbkfGivMdJI3Ag8oMfAclkHWAdDsyFgOEhQ8x7Y6MxQDFt4GVdYDQIOni8OZO4Q1RXq/qHVFuFSC3MEq6Sxg63e9WvTbKTYyQfskAKv8sFvPFsFwtuAiHcRF0CVHNHtQuBL24aA/gIuIQo3zkktuD2CFmwhkDRHHGAJEjT4B+TM5uD4gceQIMxuTs7oBUI0+AwZic3R2QauQJIKZaGxoa1je0tXnbsFwBXLS04CLgEqJeFSB/BAc9HlxEfEINhxDltxOWrg2TBFR22C/I5dfI5UR/UdNSqVvF6dYEfErAE2eYLT0gZSgK82Ny2oxl5u5C19uDlKEo5CnpPoN4/8AVwEVLCy4CLnxzwUh8DY9wW8GHi5CjVbqTFenkHpQFgJShKHRjcsb9DooDKGu+g8RoegBJg0acDuneuMoAkSNP5J3UvXGVAaI4Y4AozhggijMGKIHKu9DPvPjxXFxAZXXBoDuI5fTiwufDhdcZDMYBxNdw4XX8Qk2fI4BLfs5n7o6i5uo5mwFZcohprp6BAVJLe/UMDJBWqqtn6ZpZuOBtbVPL60Z/+VbxH1R65WU42tbm9GvqRJRA2NSTVaRd6MnqePyPkYk51gFK6j0WNnSJLLROPKNX5lN8No/UweM89jSe1CbmxAEUBw9SOgDJ18zR5ubmLZ+KL471tUfV8uF3z9rQu2mFtvE7ZvG8EIhGXQFNnVa/MpnMi2iji/GA3rrEHGsA6d9jgWTyO8imvKspjd9BH5TjQjUUrHWAdNfMWFYCai8UkrXkxBxnQUFB0SGPx+nBcrhw4XbjwuVMAMgpVHGLNRxCFT7oNJNhpr56zh5A1WuUaZyY4yspKVlVK7+QW/voOQEg40fPrgw/es40oOK9Qmk2wyyDh5hOCQGJYBRA4nyaAUXyfkCF+Qyz7gboxFxcmMwwS0BH8nupAYrvjAGiOGOAKM4YIIqzrAZENKgzIE1mEtDWGq/XZdAVgQrIoCuCy2RXBANlL6As2YPQ5yRBVov4RHSfkDSpq6QE5E8svYCoeARE3eEQY4BMOMtKQDblQY8xIB0uRZ0DpNxHpKfgdUtA8n1EEyl43RGQch/RRApedwSk3EeUHyeEP/nkk2ftkYgnguX04yIQiERMAYpEQl5cI+gI4ZLfTkva7gdddEDKK6VMPE4wBSizKXgXHRCS0FNVfpzQ4fP5XjPsq2EKUGb7augAJVbqgJT7iCZS8Ex6Su07SPt4Dss6QNJ9RFMpeBcFUOJ+8xcdkHQf0VQK3sUApHs8h2UhIGNZB0j+WZV+TIVfwrA4YUM/kxGTP6jiupJCPnlS/lnteoAo/eZtylN4M4rXpd78s3mzgEwaUiilcqKIf1alH1Msf4c4YUM/kzGTP6hoXakC6gkZjMmTqR5i27R9NZIE1Pm+GpR+87Zkv4NsyhdROr+DttWHQu4QltOH/iYJCFUJOIJ4A/x2AubPpCn95rMFkIWHWOJ+8wyQzgYWAxTfBhYJKDkfMiCbrXOAjAd5Y4AkxRnkTQ0oSVeit0sDUJxB3hggSQaDvDFAGmkHebswYcKEhdVut9ON5XC5qT6M5HbxNZH47TjSecPsIgPSD/Lm37hxY3FdMOiS+2p0CpDcV8Nluq9GhgDZxKdlnQJkapC3TgFKxyFWbnQtloykCrhmp67FTA3yZhkgfb/5JE1IFXDNpPvNYwv0Qd46hQcrS76kbakcYgmcMUAUZ1kFSPonUxEDxABlDyAyBa/zjtKYgpdtgNgeZMJZam4IMUAUpQgIbSDZE0RD8YB0hrogIF2/eYsBxR8FLx2etH7ogPT95q0FlGAUvHR40vqhA9L3m7cWUIJR8NLhSeuHDkhOdAtWVFSsqceP5EIhTzp8BL3i4z3iQZ15P3pn6fCk9UMHlFS/+U7L/LP5zoyCl4IfOiBl1Dmg91kllYFB3tR+dM4SDzSJ/6b5rQhISfWbJ5UhQJ0ZBU+xkwFA6kS3HQ2CquwNalVVNmgj2sCBQ5rAESVQb/5XLM4oeDvkFg/U4KK2Fhc1B0g79v1C9NARIXhQMFKNi7r9R6RV65M4DyIT3b6sEFVWXqHWysc0gYqt2sD8Yk2gvEyZrkzowcCPJgXvS6XFBU/hYts2XPxtHmnn1dmlpP/lS0gjz89+WV5V8pO2t0P9dSZ1lT9upq6SFt37Ejn3yhhyzjFE9UksXkDO1Q05p9sYA0SIAbISkOs8dZUfdL9rmdH3qnZ835NzHc0hctapehN9uLldt7G0AbpUxQBRlC5A1BdW6kYczJyUDCuybSydAWWRQT2kdAGivbBSP+JgxkRkWBFtY+kMKIsM6mGlCZAql8lABiMOZkxKhhXZNpbWALFIX09QmgCpcpkMpR6vP5NSMqzItrH0BvJbEtRDShMgpRtgPKnH68+wxAwrsm0svQECkL4eUhp/xeQXVhpKPV5/pvTx/PmfExlWeJagoDdAAlJlZklKEyBVLpOhdCMOZkxKhhXZNpbeAAFIXw8pTYBUuUyG0o04mDFJGVbqtrH0BghA+npI6TrEyFwmQ+lGHMyYpAwrTdtYOgMEIH09JHYmTREDRBEDRBEDRBEDRBEDRBEDRFF8QJ8PGxV32ThO0DgyGOFOwAX0BK71t4449WRdCKin5LpW2dDUBXhCXBQf0Kw83c0jWafs9je47Xb7KTLYvvQ/8LtSfqKkkGZMWI+YkutaZUNTF2o5pxDi4p7+jn1KNatdr5GHrV+Ath4ZoLskjuesRTu1P98iG+rY/vxeMqAJ8PXE/n1H8c3kvHtb7xt2AewZ3Pv69TCa4+4CZ8HA3Bn8jsrZ750AOS9e0+fu5qIBv9ygOMMLxA3w++dQjhsLO0fyC3yFv/7puHPwYZ9voP32v9j7vZbb956Tsgm8Hq77G47LGcXZx3An8FRFqTU2NHVrS4tkQOfglrsrD44cxgO6+u3GaT8Of9fziU/XcfXoo4sNG15VNZxfxI3YcQ5yBh2u/lXv1afm9PAoztACcQPoAEb456ziF9wz+tDR6bd64P7fw9pbI/bLc9+u+VM/5XEMWk+o+7N5B0dyI8q5E3jKMhvqugCfyYAgVvpvgO39eGdrAU5zZw5cdgbggybkrKrHdwBNPQ4Bxy+CnH8CPHIzwLdco+KMXyBtQNr6HW8A2H/Eu48O2Aqu/qt/chTsHL9LhPpvVDvDdYeWwnZubYQ7jqcss6GuqwIE4bdW5V+JnO0DcHJnQiN65W08iw/+l25EqwzaAtzHyFk1wJN/APCSztACcQPS1nP5WDmXw4t7GmAnt4w3ipvLm6d2huvePCb/Su5jvi6essyGuq4KUGDw0HXvV/QTWuadQezwk3f0fBc5exE7u3EzcLWSs7EaZ/wCaQPS1q/dC/DOQMnDJm6S5Cx/jtoZrnvFuPcruNoIV4+nLLOhrqsCtKenG2CH4qwS7ccTJyFnB3POAZzNqUrsTNqAtPU7K/hj5DL+q7D5zi/gzBVv9t3FO9vNf0q5z+uc7el5VynsQIA24CnLbKjrqgDVcVvP7r6uj0NydvDyTY07ryzB345DR9QeHjE0ltiZtAG09ZGPeOAxlHs0+aaPKkfdFo2NLISNAz127tr3jky8SklCQ+sJdW96YPt13N4I9zqessyGuq4EqHwz+g4qGXjVfaduGSbv2y8M6nXdynZ8AuKYPmAg/n1N5EzaANr61l9Mhn2D+VjL/Gt+ntcEG3K90DH8ITu39/a+oxtlY3g9oe7PLrv6FHcLXxdPWWZDXVcCNGZwJq7FOm6o04bsXKKnrplRWmw8mJGL1TLd6bAVgNJhY8srGQEUm6a9Sjx+PXZWM0bUjgy0mgkbTna7gyoGiCIGiCIGiCIGiCIGiCIGiCIGiKL/B1UoFjCfAoI8AAAAAElFTkSuQmCC)
:::
:::

::: {#polygenic-model-univariate .section .level2}
Polygenic model (univariate)
----------------------------

`solarPolygenic` function runs the polygenic analysis in `SOLAR` and
returns back the results in R. The basic call looks like the standard
call of `lm` function with the formula interface.

::: {#cb24 .sourceCode}
``` {.sourceCode .r}
M1 <- solarPolygenic(trait1 ~ age + sex, dat30)
M1
#> 
#> Call: solarPolygenic(formula = trait1 ~ age + sex, data = dat30)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1                Individuals:  174 
#>  
#>           H2r is 0.8061621  p = 6.1167535e-10  (Significant) 
#>         H2r Std. Error:  0.1100465 
#>  
#>  
#>  Proportion of Variance Due to All Final Covariates Is 
#>                0.0330070 
#>  
#>  Loglikelihoods and chi's are in trait1/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor, nocovar 
#>  
#>  Residual Kurtosis is -0.3603, within normal range
```
:::

The print method for the returned object of class `solarPolygenic` shows
the output produced by `SOLAR`. The tables of the results are printed by
`summary` method.

::: {#cb25 .sourceCode}
``` {.sourceCode .r}
summary(M1)
#> 
#> Call: solarPolygenic(formula = trait1 ~ age + sex, data = dat30)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1                Individuals:  174 
#>  
#>           H2r is 0.8061621  p = 6.1167535e-10  (Significant) 
#>         H2r Std. Error:  0.1100465 
#>  
#>  
#>  Proportion of Variance Due to All Final Covariates Is 
#>                0.0330070 
#>  
#>  Loglikelihoods and chi's are in trait1/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor, nocovar 
#>  
#>  Residual Kurtosis is -0.3603, within normal range 
#> 
#>  Loglikelihood Table:
#>       model    loglik
#> 1  sporadic -229.3348
#> 2 polygenic -210.8689
#> 
#>  Covariates Table:
#>   covariate      Estimate         SE Chi pval
#> 1      mean  7.7823346631 0.33112422  NA   NA
#> 2       age  0.0004591117 0.01483783  NA   NA
#> 3       sex -0.4559509166 0.31250447  NA   NA
#> 
#>  Variance Components Table:
#>   varcomp       Var        SE         pval
#> 1     h2r 0.8061621 0.1100465 6.116753e-10
#> 2      e2 0.1938379 0.1100465           NA
```
:::

::: {#testing-covariates .section .level3}
### Testing covariates

By default `SOLAR` does not perform the LRT tests for the covariates.
`covtest` argument of `solarPolygenic` function tells to `SOLAR` to make
this job.

::: {#cb26 .sourceCode}
``` {.sourceCode .r}
M2 <- solarPolygenic(trait1 ~ age + sex, dat30, covtest = TRUE)
M2
#> 
#> Call: solarPolygenic(formula = trait1 ~ age + sex, data = dat30, covtest = TRUE)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1                Individuals:  174 
#>  
#>           H2r is 0.8061621  p = 6.1167535e-10  (Significant) 
#>         H2r Std. Error:  0.1100465 
#>  
#>                                       age  p = 0.9753082  (Not Sig., but fixed) 
#>                                       sex  p = 0.1455341  (Not Sig., but fixed) 
#>  
#>  Proportion of Variance Due to All Final Covariates Is 
#>                0.0330070 
#>  
#>  Loglikelihoods and chi's are in trait1/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor, nocovar 
#>  Initial sporadic and polygenic models are s0 and p0 
#>  Constrained covariate models are named no<covariate name> 
#>  
#>  Residual Kurtosis is -0.3603, within normal range
```
:::

The test statistics and p-values are stored in `cf` slot of the returned
object `M2`.

::: {#cb27 .sourceCode}
``` {.sourceCode .r}
M2$cf
#>   covariate      Estimate         SE    Chi      pval
#> 1      mean  7.7823346631 0.33112424     NA        NA
#> 2       age  0.0004591117 0.01483782 0.0010 0.9753082
#> 3       sex -0.4559509166 0.31250456 2.1184 0.1455341
```
:::
:::

::: {#custom-kinship-matrix .section .level3}
### Custom kinship matrix

The custom kinship matrix is passed by `kinship` argument. The following
code shows such an example of `dat50` data set, which phenotype data are
stored in `phenodata` data.frame.

::: {#cb28 .sourceCode}
``` {.sourceCode .r}
M3 <- solarPolygenic(trait ~ 1, phenodata, kinship = kin)
M3
#> 
#> Call: solarPolygenic(formula = trait ~ 1, data = phenodata, kinship = kin)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait                 Individuals:  66 
#>  
#>           H2r is 0.1000000  p = 0.5000000  (Not Significant) 
#>  
#>  
#>  Loglikelihoods and chi's are in trait/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor 
#>  
#>  Residual Kurtosis is -0.6742, within normal range
```
:::

The polygenic and sporadic models have the same LRT statistics, that
means `SOLAR` was not able to make the polygenic analysis in a proper
way.

::: {#cb29 .sourceCode}
``` {.sourceCode .r}
M3$lf
#>       model    loglik
#> 1  sporadic -27.77137
#> 2 polygenic -27.77137
```
:::

This bug seems the case of unrelated individuals or the individuals
withount such ID fields as `fa` and `mo`.
:::
:::

::: {#polygenic-model-bivariate .section .level2}
Polygenic model (bivariate)
---------------------------

The bivariate polygenic analysis is run in the same way as the
univariate one.

::: {#cb30 .sourceCode}
``` {.sourceCode .r}
B1 <- solarPolygenic(trait1 + trait2 ~ 1, dat30)
B1
#> 
#> Call: solarPolygenic(formula = trait1 + trait2 ~ 1, data = dat30)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1 trait2         Individuals:  174 
#>  
#>           H2r(trait1) is 0.8218823   
#>         H2r(trait1) Std. Error:  0.1053258 
#>  
#>           H2r(trait2) is 0.6270026   
#>         H2r(trait2) Std. Error:  0.1158107 
#>  
#>           RhoE is 0.4120487   
#>         RhoE Std. Error:  0.1969379 
#>  
#>           RhoG is 0.9728759 
#>         RhoG Std. Error:  0.0374442 
#>  
#>         Derived Estimate of RhoP is 0.8045957 
#>  
#>  
#>  Loglikelihoods and chi's are in trait1.trait2/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor
```
:::

::: {#testing-correlations .section .level3}
### Testing correlations

The test of significance of environmental and genetic correlation
coefficients is indicated by `polygenic.options` argument of
`solarPolygenic` function.

::: {#cb31 .sourceCode}
``` {.sourceCode .r}
B2 <- solarPolygenic(trait1 + trait2 ~ 1, dat30, polygenic.options = "-testrhoe -testrhog")
B2
#> 
#> Call: solarPolygenic(formula = trait1 + trait2 ~ 1, data = dat30, polygenic.options = "-testrhoe -testrhog")
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1 trait2         Individuals:  174 
#>  
#>           H2r(trait1) is 0.8218823   
#>         H2r(trait1) Std. Error:  0.1053258 
#>  
#>           H2r(trait2) is 0.6270026   
#>         H2r(trait2) Std. Error:  0.1158107 
#>  
#>           RhoE is 0.4120487  p = 0.1424023 
#>         RhoE Std. Error:  0.1969379 
#>  
#>           RhoG is 0.9728759 
#>         RhoG Std. Error:  0.0374442 
#>  
#>  
#>         RhoG different from zero  p = 2.2463156e-09 
#>         RhoG different from 1.0   p = 0.2098041 
#>         Derived Estimate of RhoP is 0.8045957 
#>  
#>  
#>  Loglikelihoods and chi's are in trait1.trait2/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor
```
:::

The results are stored in `vcf` slot of the returned object.

::: {#cb32 .sourceCode}
``` {.sourceCode .r}
B2$vcf
#>       varcomp  Estimate         SE         Z        pvalZ
#> 1 h2r(trait1) 0.8218823 0.10532575  7.803242 5.215258e-03
#> 2  e2(trait1) 0.1781177 0.10532575  1.691112 1.934544e-01
#> 3 h2r(trait2) 0.6270026 0.11581068  5.414031 1.997553e-02
#> 4  e2(trait2) 0.3729974 0.11581068  3.220751 7.271026e-02
#> 5        rhog 0.9728759 0.03744417 25.982039 3.446085e-07
#> 6        rhoe 0.4120487 0.19693793  2.092277 1.480453e-01
```
:::
:::

::: {#trait-specific-covariates .section .level3}
### Trait-specific covariates

It might happen that traits in the multivariate polygenic analysis have
different covariates. Trait-specific covariates are passed to
`solarPolygenic` function according the specification of `SOLAR`.

::: {#cb33 .sourceCode}
``` {.sourceCode .r}
B3 <- solarPolygenic(trait1 + trait2 ~ sex + age(trait2), dat30)
B3
#> 
#> Call: solarPolygenic(formula = trait1 + trait2 ~ sex + age(trait2), 
#>     data = dat30)
#> 
#> File polygenic.out:
#>  Pedigree:    dat.ped 
#>  Phenotypes:  dat.phe 
#>  Trait:       trait1 trait2         Individuals:  174 
#>  
#>           H2r(trait1) is 0.7956974   
#>         H2r(trait1) Std. Error:  0.1129518 
#>  
#>           H2r(trait2) is 0.6029008   
#>         H2r(trait2) Std. Error:  0.1197632 
#>  
#>           RhoE is 0.4655421   
#>         RhoE Std. Error:  0.1736882 
#>  
#>           RhoG is 0.9676250 
#>         RhoG Std. Error:  0.0397582 
#>  
#>         Derived Estimate of RhoP is 0.8027999 
#>  
#>  
#>  Loglikelihoods and chi's are in trait1.trait2/polygenic.logs.out 
#>  Best model is named poly and null0 
#>  Final models are named poly, spor
```
:::

In the results of this analysis, `sex` covariate is common for both
traits, while `age` covariate is applied only to `trait2`. The given
model specification can be ckecked by looking at the model files from
`SOLAR`.

::: {#cb34 .sourceCode}
``` {.sourceCode .r}
tail(B3$solar$files$model$poly.mod, 3)
#> [1] "# mu = \\{t1*(<Mean(trait1)>+<bsex(trait1)>*Female) + t2*(<Mean(trait2)>"
#> [2] "# +<bsex(trait2)>*Female+<bage(trait2)>*(age-39.87931034))\\}"           
#> [3] "loglike set -361.499294"
```
:::
:::
:::

::: {#association-model .section .level2}
Association model
-----------------

`solarAssoc` function runs the association analysis in `SOLAR`, in
particular via
[mga](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#mga)
command. The genotype data are passed to `solarAssoc` function in the
form of SNP genotypes (`snpdata` argument) or covariates (`snpcovdata`
or `genocov.files` arguments). The later case assumes that the user has
converted the genotypes to covariates by some method (outside
`solareclipser` package) before computing the association model.

The function returns an object of class `solarAssoc`. This class has
additional methods such as `annotate`, `plotQQ` and `plotManh` for
exploration of the association results. The main table of the results is
`snpf` slot of the returned object, which class is `data.table` from the
R package `data.table`.

The association analysis is shown with `dat50` data set with 66
individuals and 50 SNP variants. The kinship matrix `kin` of this data
set is also used.

::: {#snp-data-by-snpdata-argument .section .level3}
### SNP data by `snpdata` argument

The option of using `snpdata` is recommended in the case when the small
number of SNPs is to be tested.

::: {#cb35 .sourceCode}
``` {.sourceCode .r}
A1 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, kinship = kin)
#> [1] "snp.geno-list"
A1
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, kinship = kin, 
#>     snpdata = genodata)
#> 
#>  Input SNP data:
#>   *  50 SNP genotypes passed by `snpdata` argument
#> 
#>  Output results of association:
#> 
#>   *  Table of association results (first 5 out of 50 rows):
#>    SNP NAv      chi     pSNP      bSNP   bSNPse   Varexp  est_maf   est_mac
#> 1:  s1  66 0.189994 0.662922  0.160860 0.369045 0.002875 0.053031  7.000026
#> 2: s10  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 3: s11  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 4: s12  66 0.909920 0.340136 -0.196598 0.206100 0.013692 0.204545 27.000006
#> 5: s13  66 0.182890 0.668901  0.118876 0.277971 0.002767 0.106061 13.999986
#>    dosage_sd
#> 1:  0.307915
#> 2:  0.244311
#> 3:  0.244311
#> 4:  0.549856
#> 5:  0.408810
#> 
#>  CPU time on 1 core(s): 00:00:01
```
:::

`summary` method reports the number of significant SNPs based on
Bonferroni correction.

::: {#cb36 .sourceCode}
``` {.sourceCode .r}
summary(A1)
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, kinship = kin, 
#>     snpdata = genodata)
#> 
#> Association model
#>  * Number of SNPs: 50 
#>  * Input format: snpdata 
#>  * Number of significal SNPs: 0 (Bonferroni correction with alpha 0.05)
```
:::

The significance threshold can be passed as `alpha` parameter of
`summary` method.

::: {#cb37 .sourceCode}
``` {.sourceCode .r}
summary(A1, alpha = 1)
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, kinship = kin, 
#>     snpdata = genodata)
#> 
#> Association model
#>  * Number of SNPs: 50 
#>  * Input format: snpdata 
#>  * Number of significal SNPs: 1 (Bonferroni correction with alpha 1)
#>    SNP NAv      chi    pSNP      bSNP   bSNPse   Varexp  est_maf  est_mac
#> 1: s43  66 6.504508 0.01076 -1.186139 0.465081 0.093852 0.030303 3.999996
#>    dosage_sd
#> 1:  0.238606
```
:::
:::

::: {#snp-data-by-snpcovdata-argument .section .level3}
### SNP data by `snpcovdata` argument

SNP covariates are passed to the association model by `snpcovdata`
parameter.

::: {#cb38 .sourceCode}
``` {.sourceCode .r}
A2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata, kinship = kin)
#> [1] "snp.geno-list"
A2
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, kinship = kin, 
#>     snpcovdata = genocovdata)
#> 
#>  Input SNP data:
#>   *  50 SNP covariates passed by `snpcovdata` argument
#> 
#>  Output results of association:
#> 
#>   *  Table of association results (first 5 out of 50 rows):
#>    SNP NAv      chi     pSNP      bSNP   bSNPse   Varexp  est_maf   est_mac
#> 1:  s1  66 0.189994 0.662922  0.160860 0.369045 0.002875 0.053031  7.000026
#> 2: s10  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 3: s11  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 4: s12  66 0.909920 0.340136 -0.196598 0.206100 0.013692 0.204545 27.000006
#> 5: s13  66 0.182890 0.668901  0.118876 0.277971 0.002767 0.106061 13.999986
#>    dosage_sd
#> 1:  0.307915
#> 2:  0.244311
#> 3:  0.244311
#> 4:  0.549856
#> 5:  0.408810
#> 
#>  CPU time on 1 core(s): 00:00:01
```
:::

The values of the SNP covariates can take any numerical values. That
means dosage format of imputed SNP data is allowed.

The previous example is repeated for the SNP covariate data, where
values of `2` are replaced by `1.9`.

::: {#cb39 .sourceCode}
``` {.sourceCode .r}
genocovdata2 <- genocovdata
genocovdata2[genocovdata2 == 2] <- 1.9

A2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata2, kinship = kin)
#> [1] "snp.geno-list"
A2
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, kinship = kin, 
#>     snpcovdata = genocovdata2)
#> 
#>  Input SNP data:
#>   *  50 SNP covariates passed by `snpcovdata` argument
#> 
#>  Output results of association:
#> 
#>   *  Table of association results (first 5 out of 50 rows):
#>    SNP NAv      chi     pSNP      bSNP   bSNPse   Varexp  est_maf   est_mac
#> 1:  s1  66 0.189994 0.662922  0.160860 0.369045 0.002875 0.053031  7.000026
#> 2: s10  66 0.050972 0.821380 -0.110595 0.489857 0.000772 0.014394  1.900008
#> 3: s11  66 0.050972 0.821380 -0.110595 0.489857 0.000772 0.014394  1.900008
#> 4: s12  66 0.927276 0.335572 -0.201592 0.209348 0.013951 0.203031 26.800026
#> 5: s13  66 0.182890 0.668901  0.118876 0.277971 0.002767 0.106061 13.999986
#>    dosage_sd
#> 1:  0.307915
#> 2:  0.232095
#> 3:  0.232095
#> 4:  0.541289
#> 5:  0.408810
#> 
#>  CPU time on 1 core(s): 00:00:01
```
:::
:::

::: {#snp-data-by-genocov.files-single-value .section .level3}
### SNP data by `genocov.files` (single value)

`genocov.files` argument is supported to pass the SNP covariate data in
files of the appropriate `SOLAR` format. This option suits for
large-scale calculations, when the real genome-wide scan is run. Both
arguments `genocov.files` and `snplists.files` are required, and
`snpmap.files` argument is optional.

::: {#cb40 .sourceCode}
``` {.sourceCode .r}
dir <- package.file("extdata", "solarAssoc", package = "solareclipser")
genocov.files <- file.path(dir, "snp.genocov")
snplists.files <- file.path(dir, c("snp.geno-list1", "snp.geno-list2"))

A3 <- solarAssoc(trait ~ 1, phenodata, 
  genocov.files = genocov.files, snplists.files = snplists.files)
A3
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, genocov.files = genocov.files, 
#>     snplists.files = snplists.files)
#> 
#>  Input SNP data:
#>   *  SNP covariates passed in 1 file(s) by `genocov.files` argument
#> 
#>  Output results of association:
#> 
#>   *  Table of association results (first 5 out of 25 rows):
#>    SNP NAv      chi     pSNP      bSNP   bSNPse   Varexp  est_maf   est_mac
#> 1:  s1  66 0.189994 0.662922  0.160860 0.369045 0.002875 0.053031  7.000026
#> 2: s10  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 3: s11  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 4: s12  66 0.909920 0.340136 -0.196598 0.206100 0.013692 0.204545 27.000006
#> 5: s13  66 0.182890 0.668901  0.118876 0.277971 0.002767 0.106061 13.999986
#>    dosage_sd
#> 1:  0.307915
#> 2:  0.244311
#> 3:  0.244311
#> 4:  0.549856
#> 5:  0.408810
#> 
#>  CPU time on 1 core(s): 00:00:00
```
:::
:::

::: {#snp-data-by-genocov.files-many-values .section .level3}
### SNP data by `genocov.files` (many values)

Many values can be passed to `genocov.files` parameter, but this case is
thought for parallel calculations.

::: {#cb41 .sourceCode}
``` {.sourceCode .r}
dir <- package.file("extdata", "solarAssoc", package = "solareclipser")
genocov.files <- file.path(dir, c("snp.genocov1", "snp.genocov2"))
snplists.files <- file.path(dir, c("snp.geno-list1", "snp.geno-list2"))

A4 <- solarAssoc(trait ~ 1, phenodata, 
  genocov.files = genocov.files, snplists.files = snplists.files)
A4
#> 
#> Call: solarAssoc(formula = trait ~ 1, data = phenodata, genocov.files = genocov.files, 
#>     snplists.files = snplists.files)
#> 
#>  Input SNP data:
#>   *  SNP covariates passed in 2 file(s) by `genocov.files` argument and 2 files(s) by `snplists.files` argument
#> 
#>  Output results of association:
#> 
#>   *  Table of association results (first 5 out of 25 rows):
#>    SNP NAv      chi     pSNP      bSNP   bSNPse   Varexp  est_maf   est_mac
#> 1:  s1  66 0.189994 0.662922  0.160860 0.369045 0.002875 0.053031  7.000026
#> 2: s10  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 3: s11  66 0.050972 0.821380 -0.105065 0.465364 0.000772 0.015152  1.999998
#> 4: s12  66 0.909920 0.340136 -0.196598 0.206100 0.013692 0.204545 27.000006
#> 5: s13  66 0.182890 0.668901  0.118876 0.277971 0.002767 0.106061 13.999986
#>    dosage_sd
#> 1:  0.307915
#> 2:  0.244311
#> 3:  0.244311
#> 4:  0.549856
#> 5:  0.408810
#> 
#>  CPU time on 1 core(s), 2 batches: 00:00:00
```
:::
:::

::: {#exploration-of-association-results .section .level3}
### Exploration of association results

QQ and Manhattan plots help to explore the association results visually.

::: {#cb42 .sourceCode}
``` {.sourceCode .r}
A5 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, snpmap = snpdata, kinship = kin)
#> [1] "snp.geno-list"
```
:::

::: {#cb43 .sourceCode}
``` {.sourceCode .r}
plot(A5)
```
:::

![](data:image/png;base64,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)

::: {#cb44 .sourceCode}
``` {.sourceCode .r}
plot(A5, "qq")
```
:::

![](data:image/png;base64,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)
:::

::: {#parallel-computation .section .level3}
### Parallel computation

Parallelization of the association analysis is an easy task with
`solareclipser`. The user only needs to specify how many cores are to be
used for the computation.

::: {#cb45 .sourceCode}
``` {.sourceCode .r}
A6 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, kinship = kin, cores = 2)
#> [1] "snp.geno-list"
```
:::

The gain in CPU time is not obvious, as the number of SNPs (50) is too
small.

::: {#cb46 .sourceCode}
``` {.sourceCode .r}
A1$assoc$tprofile$cputime.sec
#> [1] 1.213
A6$assoc$tprofile$cputime.sec
#> [1] 0.805
```
:::

`dat30` data set contains 100 SNPs (`genocovdat30` matrix). The same
test on speed up by parallel computing looks as following.

::: {#cb47 .sourceCode}
``` {.sourceCode .r}
A7 <- solarAssoc(trait1 ~ 1, dat30, 
  snpcovdata = genocovdat30, snpmap = mapdat30)
#> [1] "snp.geno-list"
```
:::

::: {#cb48 .sourceCode}
``` {.sourceCode .r}
A8 <- solarAssoc(trait1 ~ 1, dat30, 
  snpcovdata = genocovdat30, snpmap = mapdat30,
  cores = 2)
#> [1] "snp.geno-list"
```
:::

::: {#cb49 .sourceCode}
``` {.sourceCode .r}
A7$assoc$tprofile$cputime.sec
#> [1] 1.995
A8$assoc$tprofile$cputime.sec
#> [1] 1.064
```
:::
:::
:::

::: {#linkage-model .section .level2}
Linkage model
-------------

`solarMultipoint` function performs the linkage analysis by calling
`SOLAR` command
[multipoint](http://helix.nih.gov/Documentation/solar-6.6.2-doc/91.appendix_1_text.html#multipoint).
The function looks for multipoint IBD matrices in the directory pointed
out by `mibddir` argument. The analysis can be customized by combination
of options via `multipoint.options` and `multipoint.settings` arguments.

Phenotypes in `dat30` data set and IBD matrices are placed in
[inst/extdata/solarOutput](https://github.com/ugcd/solarius/tree/master/inst/extdata/solarOutput)
directory of the package.

::: {#basic-linkage-model .section .level3}
### Basic linkage model

The basic model requires only the path to a directory with IBD matrices
(`mibddir` argument). By default, all `SOLAR` specific files, which
represent IBD matrices, will be used in the analysis.

::: {#cb50 .sourceCode}
``` {.sourceCode .r}
mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solareclipser")  
list.files(mibddir)
#>  [1] "mibd.2.0.gz" "mibd.2.1.gz" "mibd.2.2.gz" "mibd.2.3.gz" "mibd.2.4.gz"
#>  [6] "mibd.2.5.gz" "mibd.5.0.gz" "mibd.5.1.gz" "mibd.5.2.gz" "mibd.5.3.gz"
#> [11] "mibd.5.4.gz" "mibd.5.5.gz"
```
:::

::: {#cb51 .sourceCode}
``` {.sourceCode .r}
L1 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir)
L1
#> 
#> Call: solarMultipoint(formula = trait1 ~ 1, data = dat30, mibddir = mibddir)
#> 
#>  Input IBD data:
#>   *  directory /tmp/Rtmp8xhLLr/Rinst2e9d41898f786/solareclipser/extdata/solarOutput/solarMibdsCsv
#> 
#>  Output results of association:
#> 
#>   * Table of association results (first 5 out of 12 rows):
#>   chr pos    LOD  Loglike      H2r     H2q1
#> 1   2   0 1.4888 -208.503 0.448464 0.403780
#> 2   2   1 2.4232 -206.352 0.244560 0.614310
#> 3   2   2 3.2759 -204.388 0.123184 0.742126
#> 4   2   3 3.5569 -203.741 0.097187 0.763608
#> 5   2   4 3.1610 -204.653 0.179701 0.668166
#> 
#>  CPU time on 1 core(s): 00:00:01
```
:::

Methods `print`, `summary` and `plot` help to explore the results of
analysis. The data table of results is stored in `lodf` slot of the
returned object.

::: {#cb52 .sourceCode}
``` {.sourceCode .r}
summary(L1)
#> 
#> Call: solarMultipoint(formula = trait1 ~ 1, data = dat30, mibddir = mibddir)
#> 
#> Multipoint model
#>  * Number of used markers: 12 
#>  * Number of passes: 1 
#>  * Maximum LOD score: 3.56 
#>   -- chr: 2 
#>   -- position: 3 cM
```
:::

::: {#cb53 .sourceCode}
``` {.sourceCode .r}
plot(L1)
```
:::

![](data:image/png;base64,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)

The results for chromosomes 2 and 5 are identical, since the IBD
matrices for chromosome 2 are duplicates of the matrices for chromosome
5.
:::

::: {#linkage-options .section .level3}
### Linkage options

The following code shows an example of a custom linkage analysis. In
particular, the linkage is run on chromosome 5 with a 2 cM interval,
Fine mapping is disabled and the second pass is enables via
`multipoint.options` argument.

::: {#cb54 .sourceCode}
``` {.sourceCode .r}
L2 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir, chr = 5, interval = 2, 
  multipoint.settings = "finemap off", multipoint.options = "3")
summary(L2)
#> 
#> Call: solarMultipoint(formula = trait1 ~ 1, data = dat30, mibddir = mibddir, 
#>     chr = 5, interval = 2, multipoint.options = "3", multipoint.settings = "finemap off")
#> 
#> Multipoint model
#>  * Number of used markers: 3 
#>  * Number of passes: 2 
#>  * Maximum LOD score: 3.28 
#>   -- chr: 5 
#>   -- position: 2 cM
```
:::

The LOD scores for the second path (stored in `lodf2`) are plotted with
the code given below.

::: {#cb55 .sourceCode}
``` {.sourceCode .r}
plot(L2, pass = 2)
```
:::

![](data:image/png;base64,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)
:::

::: {#bivariate-linkage .section .level3}
### Bivariate linkage

The bivariate linkage analysis is performed in the same way as the
bivariate polygenic analysis.

::: {#cb56 .sourceCode}
``` {.sourceCode .r}
L3 <- solarMultipoint(trait1 + trait2 ~ 1, dat30, mibddir = mibddir, chr = 5, 
  interval = 2, multipoint.settings = "finemap off")
summary(L3)
#> 
#> Call: solarMultipoint(formula = trait1 + trait2 ~ 1, data = dat30, 
#>     mibddir = mibddir, chr = 5, interval = 2, multipoint.settings = "finemap off")
#> 
#> Multipoint model
#>  * Number of used markers: 3 
#>  * Number of passes: 1 
#>  * Maximum LOD score: 2.47 
#>   -- chr: 5 
#>   -- position: 2 cM
```
:::

::: {#cb57 .sourceCode}
``` {.sourceCode .r}
L4 <- solarMultipoint(trait1 + trait2 ~ age + sex(trait1), dat30, 
  mibddir = mibddir, chr = 5, interval = 2, multipoint.settings = "finemap off")
summary(L4)
#> 
#> Call: solarMultipoint(formula = trait1 + trait2 ~ age + sex(trait1), 
#>     data = dat30, mibddir = mibddir, chr = 5, interval = 2, multipoint.settings = "finemap off")
#> 
#> Multipoint model
#>  * Number of used markers: 3 
#>  * Number of passes: 1 
#>  * Maximum LOD score: 2.68 
#>   -- chr: 5 
#>   -- position: 2 cM
```
:::
:::

::: {#parallel-computation-1 .section .level3}
### Parallel computation

Parallelization of the linkage analysis is again parametrized with
`cores` argument, but the chromosome values are expected to be a vector
of values. That means that each core does the job on a single
chromosome.

::: {#cb58 .sourceCode}
``` {.sourceCode .r}
L5 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir,
  chr = c(2, 5), cores = 2)
```
:::

::: {#cb59 .sourceCode}
``` {.sourceCode .r}
L1$multipoint$tprofile$cputime.sec
#> [1] 1.063
L5$multipoint$tprofile$cputime.sec
#> [1] 0.662
```
:::
:::
:::
:::

::: {#miscellaneous-material .section .level1}
Miscellaneous material
======================

::: {#solar-references .section .level2}
`SOLAR` references
------------------

A list of references given on the official web page of `SOLAR`
[solar-eclipse-genetics.org](https://solar-eclipse-genetics.org):

-   The main reference for SOLAR, including theoretical explanations of
    the variance component linkage method and the approximate multipoint
    IBD calculations in pedigrees, is [\[\@Almasy1998\]]{.citation};
-   Bivariate Quantitative trait linkage is described in the following
    papers: [\[\@Almasy1997\]]{.citation},
    [\[\@Williams1999\]]{.citation};
-   The Liability Threshold model for discrete traits is described in
    the preceding paper as well as the following one:
    [\[\@Duggirala1997\]]{.citation};
-   Gene By Environment Interaction is discussed in:
    [\[\@Towne1997\]]{.citation}, [\[\@Blangero2009\]]{.citation};
-   An examination of LOD adjustment is given in:
    [\[\@Blangero2000\]]{.citation}, [\[\@Blangero2001\]]{.citation};
-   Additional references include: [\[\@Blangero1997\]]{.citation},
    [\[\@Williams2004\]]{.citation}, [\[\@Williams1999b\]]{.citation}.
:::
:::

::: {#r-session-info .section .level1}
R session info
==============

::: {#cb60 .sourceCode}
``` {.sourceCode .r}
sessionInfo()
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Pop!_OS 22.04 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=C              
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gridExtra_2.3       ggplot2_3.4.4       plyr_1.8.9         
#> [4] solareclipser_0.3.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.0  xfun_0.43         bslib_0.5.1       qqman_0.1.9      
#>  [5] lattice_0.20-45   colorspace_2.1-0  vctrs_0.6.4       generics_0.1.3   
#>  [9] htmltools_0.5.8.1 yaml_2.3.7        utf8_1.2.3        rlang_1.1.1      
#> [13] jquerylib_0.1.4   pillar_1.9.0      glue_1.6.2        withr_2.5.1      
#> [17] calibrate_1.7.7   foreach_1.5.2     lifecycle_1.0.3   munsell_0.5.0    
#> [21] gtable_0.3.4      codetools_0.2-18  evaluate_0.22     labeling_0.4.3   
#> [25] knitr_1.46        fastmap_1.1.1     doParallel_1.0.17 parallel_4.1.2   
#> [29] fansi_1.0.5       highr_0.10        Rcpp_1.0.12       scales_1.2.1     
#> [33] cachem_1.0.8      jsonlite_1.8.7    farver_2.1.1      digest_0.6.33    
#> [37] dplyr_1.1.3       grid_4.1.2        quadprog_1.5-8    kinship2_1.9.6.1 
#> [41] cli_3.6.1         tools_4.1.2       magrittr_2.0.3    sass_0.4.7       
#> [45] tibble_3.2.1      crayon_1.5.2      pkgconfig_2.0.3   MASS_7.3-55      
#> [49] Matrix_1.4-0      data.table_1.14.8 rmarkdown_2.25    rstudioapi_0.15.0
#> [53] iterators_1.0.14  R6_2.5.1          compiler_4.1.2
```
:::
:::

::: {#license .section .level1}
License
=======

This document is licensed under the Creative Commons Attribution 4.0
International Public License.

[![Creative Commons
License](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFgAAAAfCAMAAABUFvrSAAAAAXNSR0IB2cksfwAAAARnQU1BAACxjnz7UZMAAAAgY0hSTQAAeiUAAICDAAD5/wAAgOkAAHUwAADqYAAAOpgAABdvkl/FRgAAAW5QTFRFAAAAAAAA////////////7+/v39/f1tXV09bS0tXS0tXR0dTR0dTQ0NTQ0NPPz9PPztLOztHNzdHNzdHMz8/PzdDMzNDMzNDLzM/Ly8/Ly8/Ky87Kys3Jyc3Jyc3IyMzIyMzHx8vHxsrGxsrFxcnFxcnExMnExMjDw8jDxMfDw8fCwsfCwcXAwMXAwMW/wMS/v8S+v8O+vsO+vsK9vcK9vcK8v7+/vMG8vMG7vMC8u8C7u8C6ur+6ur+5ub65ub64uL23t7y2tru1tbq0tLqztLmzs7iysrixsrexsbewsbawsLavsLWvr7Wur7SusLOvrrStrrOtr7KvrbOsrLKrr6+vq7Gqn6OenqCdn5+flpmWk5iTkZSRkZORj4+PiYyJhIaEhIWEgoWCgICAfX98fH98eXx5cHJvcHBwYGBgXV5dUFFQUFBQQ0RDQEBAPj8+NTY1MjMxMDAwKSkpKCkoICAgGxsbEBAQDg4ODQ4N2y3MbAAAAAR0Uk5T/wAKDnDBpeYAAALVSURBVHjatZX9V9JgFMdvNQh1Tme2zU1othSl5WsQZoqIrxmmqaUpqS2JKczibfLfd5+hHCYybB6/B859ftg+5+57Xx54Ag+iR/hPxGPRSXVYkf2SwHGC2C8rQXUiOhdPLK992tze3vn6/wKLOxseCykBkWXoNp+vrYNh+wJKaDQ8S8jJze0dl2TkvgkOSCx9kqsAqpI7prtFeUh9i+SV9eRWY8rpAqAKaWdwDLkyz6TKUFP5kOECg2p4Lr60ut6Q8kERwMhkDIDigRM4OhaUOfocT6auRSKabuLxfORZYGg0GltAM26k/O0Ssl4K5c3C5YEDeDI0wBOuqXmoqvYRXRzh5NDEzHwCU7aDi6DjMwDkQSg6gFVFYpCb91I1efKYMyMqahhT/rhh8yIN2adVMCoL6ebg4RdsCrmYrp182B0Ijs/EF1eTW/XgAtRl4IVCc7Ai0mUweymbPCaU6T5FfRdLrNi9AINC6QA6iQY0B8s9JwCanUvsg2NWDk3Nohd2cIZCAYrEjAPY35UDE42IZA1Dq4Yf+IoJOaafeLG0tuEOLNEV0BEIRMauFXTysZUOcXD0/fziWtKdFUI7AKaYh5I3UtrPY3vslvYtL9oERY1iwyG4oXgALYvH+wAiFGX5XAvWF/i4l41gbDfKod3cgu0DksEBuYsVnt6ShlaQ0GCF00g7F++qaleBZFOhpcbiEXIRIK9n8q2WkL/rHEyq1m67hhWoErabPzj+odpubtam3JNCL24fkNdTcwvLZPLcLPpX1kh7HUbaJRiX0N5tSyjFBoYncPDIEnIFhlOJOUNy3RryWmtT+gn31G+hs37RezQT4A/Nn8K99YvrPLtxNdHPv0MzVeeDHFqTBWav/jJNMXxzLlBXP+oubkgsnbq+/g9pVnLy4Xptwl3094vY00W3+3ztNMP2fb4AB9mMaK2Lo2lJ4HlBmj4iWLcZP7wePwz20T8rQcP0CuFIbQAAAABJRU5ErkJggg==)](http://creativecommons.org/licenses/by/4.0/)
:::

::: {#references .section .level1}
References
==========
:::
