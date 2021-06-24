#' @title
#' Omni-Optimizer
#'
#' @description Simple interface to the C-implementation of the Omni-optimizer
#' by Deb and Tiwari [1,2]. The algorithm \dQuote{is designed as a generic
#' multi-objective, multi-optima optimizer} [2].
#'
#' @details The function expects a real-valued benchmark function from package
#' \CRANpkg{smoof}, the population size and the number of generations (the only
#' stopping condition) as mandatory arguments. Besides there are various
#' parameters that can be adjusted (see the referenced papers for an in-depth
#' explanation of the algorithms working principles).
#'
#' The original C-code can be found at the
#' \href{http://www.coin-lab.org/content/source_codes.html}{COIN laboratory website}.
#'
#' @references
#' [1] Kalyanmoy Deb, Santosh Tiwari: Omni-optimizer: A generic evolutionary algorithm
#' for single and multi-objective optimization. European Journal of Operations
#' Research 185(3): 1062-1087.
#'
#' [2] Kalyanmoy Deb, Santosh Tiwari: Omni-optimizer: A Procedure for Single and
#' Multi-objective Optimization. In: Proceedings of the Evolutionary Multi-Criterion
#' Conference (EMO) 2005: 47-61.
#'
#' @keywords optimize
#'
#' @param fn [\code{function}]\cr
#'   Multi-objective function of type \code{smoof_function} (see \CRANpkg{smoof}).
#' @param pop.size [\code{integer(1)}]\cr
#'   Population size. Must be a multiple of 4.
#'   The default is 4.
#' @param n.gens [\code{integer(1)}]\cr
#'   The number of generations (stopping condition).
#'   Defaults to 100.
#' @param p.cross [\code{numeric(1)}]\cr
#'   Probability of crossover (within \eqn{[0.6, 1.0]}).
#'   Defaults to 0.6.
#' @param p.mut [\code{numeric(1)}]\cr
#'   Probablity of mutation (within \eqn{[0,1]}).
#'   Default to \eqn{1/n} where \eqn{n} is the number of decision variables of
#'   \code{fn}.
#' @param eta.cross [\code{integer(1)}]\cr
#'   Value of distribution index for crossover in \eqn{\{5, \ldots, 20\}}.
#'   Default is 20.
#' @param eta.mut [\code{integer(1)}]\cr
#'   Value of distribution index for mutation in \eqn{\{5, \ldots, 50\}}.
#'   Default is 20.
#' @param mate [\code{character(1)}]\cr
#'   Choice for selection restriction. Either \dQuote{normal} for normal selection
#'   or \dQuote{restricted} for restricted selection.
#' @param delta [\code{numeric(1)}]\cr
#'   Value \eqn{\delta \in (0.0,1.0)} for loose domination.
#'   Default is \eqn{0.001}.
#' @param var.space.niching [\code{logical(1)}]\cr
#'   Use variable space niching?
#'   Default is \code{FALSE}.
#'   Note that at least on of \code{var.space.niching} or \code{obj.space.niching}
#'   must be \code{TRUE}.
#' @param obj.space.niching [\code{logical(1)}]\cr
#'   Use objective space niching?
#'   Default is \code{TRUE}.
#'   Note that at least on of \code{var.space.niching} or \code{obj.space.niching}
#'   must be \code{TRUE}.
#' @param init [\code{character(1)}]\cr
#'   This parameter determines how to initialize the  population: \dQuote{random}
#'   (the default) for uniform random generation and \dQuote{lhs} for
#'   Latin-Hypercube-Sampling (LHS).
#' @param frequency [\code{integer(1)}]\cr
#'   Frequency with which the population information is to be stored.
#'   Defaults to 1.
#' @param seed [\code{numeric(1)}]\cr
#'   Single numeric value in \eqn{[0,1]}.
#'   Defaults to a random number within this interval.
#' @param envir [\code{environment}]\cr
#'   This parameter is required for calling R functions from C.
#'   Do not change the default unless you know what you are doing!
#' @return [\code{List}] List with the following entries:
#' \describe{
#'   \item{pareto.set}{Matrix of non-dominated points (each column is a point).}
#'   \item{pareto.front}{Matrix of non-dominated point objective values
#'   (each column is a point).}
#' }
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
omniopt = function(
  fn,
  pop.size = 4L,
  n.gens = 100L,
  p.cross = 0.6,
  p.mut = 1/smoof::getNumberOfParameters(fn),
  eta.cross = 20L,
  eta.mut = 20L,
  mate = "normal",
  delta = 0.001,
  var.space.niching = FALSE,
  obj.space.niching = TRUE,
  init = "random",
  frequency = 1,
  seed = runif(1),
  envir = environment()) {

  checkmate::assert_class(fn, "smoof_function")
  pop.size = checkmate::asInt(pop.size, lower = 4L)
  if (pop.size<4 || (pop.size %% 4) != 0) {
    stop("[omniopt] pop.size needs to be at least 4 and a multiple of 4.\n")
  }
  n.gens = checkmate::asInt(n.gens, lower = 1L)
  checkmate::assert_number(p.cross, lower = 0.6, upper = 1)
  checkmate::assert_number(p.mut, lower = 0, upper = 1)
  eta.cross = checkmate::asInt(eta.cross, lower = 5L, upper = 20L)
  eta.mu = checkmate::asInt(eta.mut, lower = 5L, upper = 50L)

  mate.options = c("normal", "restricted")
  checkmate::assert_choice(mate, mate.options)
  mate = which(mate == mate.options) - 1L # convert to 0 and 1

  checkmate::assert_number(delta, lower = 0, upper = 1)

  checkmate::assert_flag(var.space.niching)
  checkmate::assert_flag(obj.space.niching)

  if (!var.space.niching & !obj.space.niching) {
    stop("[omniopt] At least one of var.space.niching and obj.space.niching must be set to TRUE.")
  }

  init.options = c("random", "lhs")
  checkmate::assert_choice(base::tolower(init), init.options)
  init = which(init == init.options) - 1L # convert to 0 and 1

  frequency = checkmate::asInt(frequency, lower = 1L)
  if (frequency < 1 || frequency > n.gens) {
    stop(sprintf("Parameter frequency should be in the range 1 to n.gens (%i).", n.gens));
  }

  checkmate::assert_number(seed, lower = 0, upper = 1)

  n.objectives = smoof::getNumberOfObjectives(fn)
  dimension = smoof::getNumberOfParameters(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)

  rawres = .Call("omnioptC",
    fn,
    as.integer(n.objectives),
    as.integer(dimension),
    as.double(lower),
    as.double(upper),
    as.integer(pop.size),
    as.integer(n.gens),
    as.double(p.cross),
    as.double(p.mut),
    as.integer(eta.cross),
    as.integer(eta.mut),
    as.integer(mate),
    as.double(delta),
    as.integer(var.space.niching),
    as.integer(obj.space.niching),
    as.integer(init),
    as.integer(frequency),
    as.double(seed),
    envir)

  # convert in "ecr"-style format
  names(rawres) = c("pareto.set", "pareto.front")
  rawres$pareto.front = t(matrix(rawres$pareto.front, byrow = TRUE, ncol = n.objectives))
  rawres$pareto.set = t(matrix(rawres$pareto.set, byrow = TRUE, ncol = dimension))
  return(rawres)
}
