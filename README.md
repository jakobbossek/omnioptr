
# omnioptr: Interface to OmniOptimizer

## Introduction

…

## Example

``` r
library(omnioptr)
library(smoof)

fn = smoof::makeDTLZ1Function(dimensions = 2L, n.objectives = 2L)
fn = smoof::addLoggingWrapper(fn, logg.x = TRUE, logg.y = TRUE)

set.seed(52357) # reproducibility
res = omniopt(fn, pop.size = 20L, n.gens = 100L)
print(res)
print(smoof::getLoggedValues(fn, compact = TRUE))
plot(t(res$pareto.front))
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
