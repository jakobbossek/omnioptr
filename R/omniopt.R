#' @title
#' OmniOptimizer.
#'
#' @description Simple interface to the C-implementation of the OmniOptimizer.
#'
#' @param fn [\code{function}]\cr
#'   Multi-objective function of type \code{smoof_function} (see \CRANpkg{smoof}).
#' @param pop.size [\code{integer(1)}]\cr
#'   Population size. Must be a multiple of 4.
#'   The default is 4.
#' @param n.gens [\code{integer(1)}]\cr
#'   The number of generations (stopping condition).
#'   Defaults to 100.
#' @param seed [\code{numeric(1)}]\cr
#'   Single numeric value in \eqn{[0,1]}.
#'   Defaults to a random number within the interval.
#' @param envir [\code{environment}]\cr
#'   This parameter is required for calling R functions from C.
#'   Do not modify unless you know what you are doing!
#' @return [\code{List}] List with the following entries: ...
#'
#' @examples
#' \dontrun{
#' library(smoof)
#'
#' fn = smoof::makeDTLZ1Function(dimensions = 2L, n.objectives = 2L)
#' fn = smoof::addLoggingWrapper(fn, logg.x = TRUE, logg.y = TRUE)
#'
#' set.seed(52357) # reproducibility
#' res = omniopt(fn, pop.size = 20L, n.gens = 1000L)
#' print(res)
#' print(smoof::getLoggedValues(fn, compact = TRUE))
#' plot(t(res$pareto.front))
#' }
#' @export
omniopt = function(fn, pop.size = 4L, n.gens = 100L, seed = runif(1), envir = environment()) {
  n.objectives = smoof::getNumberOfObjectives(fn)
  dimension = smoof::getNumberOfParameters(fn)

  rawres = .Call("omnioptC",
    fn,
    as.integer(n.objectives),
    as.integer(dimension),
    as.integer(pop.size),
    as.integer(n.gens),
    as.double(seed),
    envir)

  names(rawres) = c("pareto.set", "pareto.front")
  rawres$pareto.front = t(matrix(rawres$pareto.front, byrow = TRUE, ncol = n.objectives))
  rawres$pareto.set = t(matrix(rawres$pareto.set, byrow = TRUE, ncol = dimension))
  return(rawres)
}
