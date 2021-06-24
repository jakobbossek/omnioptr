
# omnioptr: Interface for the Omni-Optimizer Algorithm

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/jakobbossek/omnioptr/workflows/R-CMD-check/badge.svg)](https://github.com/jakobbossek/omnioptr/actions)
[![codecov](https://codecov.io/gh/jakobbossek/omnioptr/branch/main/graph/badge.svg?token=88YGQRYJ4W)](https://codecov.io/gh/jakobbossek/omnioptr)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/omnioptr)](https://cran.r-project.org/package=omnioptr)
<!-- badges: end -->

## Introduction

Simple interface to the C-implementation of the Omni-optimizer by Deb
and Tiwari (1,2). The algorithm “is designed as a generic
multi-objective, multi-optima optimizer” (2).

1.  Kalyanmoy Deb, Santosh Tiwari: Omni-optimizer: *A generic
    evolutionary algorithm for single and multi-objective optimization*.
    European Journal of Operations Research 185(3): 1062-1087.

2.  Kalyanmoy Deb, Santosh Tiwari: *Omni-optimizer: A Procedure for
    Single and Multi-objective Optimization*. In: Proceedings of the
    Evolutionary Multi-Criterion Conference (EMO) 2005: 47-61.s

## Example

``` r
library(omnioptr)
library(smoof)

# Single-Objective Example (see reference (2), Section 4.2)
# ===

set.seed(1) # reproducibility
fn = smoof::makeHimmelblauFunction()
res = omniopt(fn, 100, 200, var.space.niching = TRUE, delta = 0.001,
  verbose = FALSE envir = environment())

plot(fn)
points(t(res$pareto.set))

# Multi-Objective Example
# ===

set.seed(1) # reproducibility
fn = smoof::makeZDT2Function(dimension = 4L)
res = omniopt(fn, 100, 1000, p.cross = 0.9, verbose = FALSE)
plot(t(res$pareto.front))
pairs(t(res$pareto.set))
```

## Installation Instructions

The package will be available at [CRAN](http://cran.r-project.org) *when
it is done*. If you are interested in trying out and playing around with
the current github developer version use the
[devtools](https://github.com/hadley/devtools) package and type the
following command in R:

``` r
remotes::install_github("jakobbossek/omnioptr")
```

## Getting help

Please address questions and missing features about the *omnioptr* as
weell as annoying bug reports in the [issue
tracker](https://github.com/jakobbossek/omnioptr/issues). Pay attention
to explain your problem as good as possible. At its best you provide an
example, so I can reproduce your problem quickly. Please avoid sending
e-mails.
