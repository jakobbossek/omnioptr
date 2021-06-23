#' @title OmniOptimizer.
#'
#' @description ...
#'
#' @export
omniopt = function(fn, pop.size, n.gens, envir = environment()) {
  return(.Call("omnioptC",
    fn,
    as.integer(pop.size),
    as.integer(n.gens),
    envir))
}
