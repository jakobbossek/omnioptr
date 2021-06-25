library(devtools)
library(smoof)

load_all(".")

set.seed(1)
fn = smoof::makeHimmelblauFunction()
res = omniopt(fn, 100, 1000, var.space.niching = TRUE, frequency = 1, delta = 0.001, envir = environment())
plot(fn)
points(t(res$dec))

for (i in 1:length(res$history)) {
  e = res$history[[i]]
  if (is.null(e))
    next
  plot(fn)#, main = sprintf("Iteration %i", i))
  points(t(e$dec))
  Sys.sleep(0.1)
}
stop()


fn = smoof::makeZDT2Function(dimensions = 4L)
res = omniopt(fn, 100, 1000, var.space.niching = TRUE, delta = 0.001, p.cross = 0.9, envir = environment())

plot(t(res$obj))
pairs(t(res$dec))
