library(devtools)
library(smoof)

load_all(".")

fn = smoof::makeZDT2Function(dimensions = 4L)
#fn = smoof::makeMMF2Function()
res = omniopt(fn, 100, 200, seed = 0.1, envir = environment())
print(res)
plot(t(res$pareto.front))

