# TODO

- [ ] generate documentation.
- [ ] check out git hooks - https://www.git-scm.com/docs/githooks#_post_checkout 
- [ ] set project path for the bundled solar
- [ ] write a test for solar.lib.R. It fails to load sym. link in project
      root.

## BUGS

The following is what doesn't work following the [tutorial](https://ugcd.github.io/solarius/vignettes/tutorial.html). 🔴 is a breaking bug and 🟡 is a non-breaking bug that may astonish user.


- 🔴 [4.2.1](https://ugcd.github.io/solarius/vignettes/tutorial.html#plot-kinship-matrix)

    ```txt
    > plotKinship2(2*kin)
    Error in h(simpleError(msg, call)) : 
      error in evaluating the argument 'x' in selecting a method for function 'image': object 'kin' not found
    Called from: h(simpleError(msg, call))
    ```

- 🟡 [4.2.3](https://ugcd.github.io/solarius/vignettes/tutorial.html#transform-traits)

    - `library("gridExtra")`

      - wasn't installed, so not in solareclipser dependency list, however after
        install it works. Must be included as dependency for examples, or
        suggested.

- 🟡 [4.3.2](https://ugcd.github.io/solarius/vignettes/tutorial.html#custom-kinship-matrix)

    - When following the tutorial `phenodata` is not loaded at this point and
    gives an error. Load the data first.

- 🟡 [4.5.3](https://ugcd.github.io/solarius/vignettes/tutorial.html#snp-data-by-genocov.files-single-value)

    - change package in `dir <- package.file("extdata", "solarAssoc", package = "solarius")` to `solareclipser`.
      This is just for when using the tutorial. Other references have been
      replaced.

- 🟡 [4.5.4](https://ugcd.github.io/solarius/vignettes/tutorial.html#snp-data-by-genocov.files-many-values)

    - same as 4.5.3

- 🔴 [4.5.5](https://ugcd.github.io/solarius/vignettes/tutorial.html#exploration-of-association-results)

    ```txt
    > plot(A5)
    Error in plotManh(x, ...) : 
      requireNamespace("qqman", quietly = TRUE) is not TRUE
    Called from: plotManh(x, ...)
    ```

    - `qqman` should be a dependency. NOTE: This was added. Needs a test.

- 🟡 [6](https://ugcd.github.io/solarius/vignettes/tutorial.html#r-session-info) has
    a dependency listing that differs from other documentation.


## Other notes

I keep finding errors caused by checks for matrix. See list:

```txt
R/solar.lib.R:  #stopifnot(class(mat) == "matrix")
R/solar.lib.R:  #stopifnot(class(mat) == "matrix")
R/solar.lib.R:  stopifnot(class(map) == "data.frame")
R/solar.lib.R:  stopifnot(class(kf) == "data.frame")
R/solar.lib.R:  stopifnot(class(kmat) == "matrix")
R/solar.lib.R:  stopifnot(class(kf) == "data.frame")
R/solar.lib.R:  stopifnot(class(pf) == "data.frame")
R/multipoint.lib.R:  stopifnot(class(mf) == "data.frame")
R/classSolarMultipoint.R:  stopifnot(class(lodf) == "data.frame")
R/solarBayesAvg.R:  stopifnot(class(data) == "data.frame")
R/solarPolyAssoc.R:  #  stopifnot(class(snpcovdata) == "matrix")
R/plot.lib.R:  stopifnot(class(data) == "data.frame")
R/solarPolygenic.R:  stopifnot(class(data) == "data.frame")
R/solarAssoc.R:    stopifnot(class(mga.files) == "list")
R/solarAssoc.R:  #  stopifnot(class(snpdata) == "matrix")
R/solarAssoc.R:  #  #stopifnot(class(snpcovdata) == "matrix")
R/transforms.lib.R:  stopifnot(class(transforms) == "character")
R/transforms.lib.R:  stopifnot(class(transform) == "character")
R/transforms.lib.R:  stopifnot(class(x) %in% c("integer", "numeric"))
```

## R notes

1. CLI run R: `Rscript main.R`


## References

- https://r6.r-lib.org/articles/Introduction.html
- https://adv-r.hadley.nz/oo-tradeoffs.html
- https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure
- https://debruine.github.io/tutorials/your-first-r-package-with-unit-tests.html#unit-tests



