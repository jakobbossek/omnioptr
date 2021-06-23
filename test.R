library(devtools)
library(smoof)

load_all(".")

fn = smoof::makeDTLZ1Function(dimensions = 2L, n.objectives = 2L)

res = omniopt(fn, 12, 100, envir = environment())
prin(res)

