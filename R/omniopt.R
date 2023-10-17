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
#' Optimization Conference (EMO) 2005: 47-61.
#'
#' @keywords optimize
#'
#' @param fn [\code{function}]\cr
#'   Single- or multi-objective function of type \code{smoof_function}
#'   (see \CRANpkg{smoof}) with continuous decision space.
#' @param pop.size [\code{integer(1)}]\cr
#'   Population size. Must be a multiple of 4.
#'   The default is 4.
#' @param n.gens [\code{integer(1)}]\cr
#'   The number of generations (the only stopping condition).
#'   Defaults to 100.
#' @param p.cross [\code{numeric(1)}]\cr
#'   Probability of crossover of real-valued variables (within \eqn{[0.6, 1.0]}).
#'   Defaults to 0.6.
#' @param p.mut [\code{numeric(1)}]\cr
#'   Probablity of mutation of real-valued variables (within \eqn{[0,1]}).
#'   Default to \eqn{1/n} where \eqn{n} is the number of decision variables of
#'   \code{fn}.
#' @param p.cross.bin [\code{numeric(1)}]\cr
#'   Probability of crossover of binary-valued variables (within \eqn{[0.6, 1.0]}).
#'   Defaults to 0.6.
#' @param p.mut.bin [\code{numeric(1)}]\cr
#'   Probablity of mutation of binary-valued variables (within \eqn{[0,1]}).
#'   Default to \eqn{1/n} where \eqn{n} is the number of decision variables of
#'   \code{fn}.
#' @param eta.cross [\code{numeric(1)}]\cr
#'   Value of distribution index for crossover in \eqn{[5, \ldots, 20]}.
#'   Default is 20.
#' @param eta.mut [\code{numeric(1)}]\cr
#'   Value of distribution index for mutation in \eqn{[5, , 50]}.
#'   Default is 20.
#' @param mate [\code{character(1)}]\cr
#'   Choice for selection restriction. Either \dQuote{normal} for normal selection
#'   or \dQuote{restricted} for restricted selection.
#' @param delta [\code{numeric(1)}]\cr
#'   Value \eqn{\delta \in [0.0,1.0]} for loose domination.
#'   Default is \eqn{0.001}.
#' @param var.space.niching [\code{logical(1)}]\cr
#'   Use variable space niching?
#'   Default is \code{FALSE}.
#'   Note that at least one of \code{var.space.niching} or \code{obj.space.niching}
#'   must be \code{TRUE}.
#' @param obj.space.niching [\code{logical(1)}]\cr
#'   Use objective space niching?
#'   Default is \code{TRUE}.
#'   Note that at least one of \code{var.space.niching} or \code{obj.space.niching}
#'   must be \code{TRUE}.
#' @param init [\code{character(1)}]\cr
#'   This parameter determines how to initialize the population: \dQuote{random}
#'   (the default) for uniform random generation and \dQuote{lhs} for
#'   Latin-Hypercube-Sampling (LHS).
#' @param frequency [\code{integer(1)}]\cr
#'   Frequency with which the population information (points and objetive function values)
#'   is to be stored.
#'   Defaults to 10, i.e., information is stored in every 10th generation.
#' @param seed [\code{numeric(1)}]\cr
#'   Seed for random number generaator. Must be a scalar numeric value in \eqn{[0,1]}.
#'   Defaults to a random number within this interval.
#' @param verbose [\code{logical(1)}]\cr
#'   If \code{TRUE} the algorithm is verbose, i.e., it prints some informative
#'   messages in the course of optimization.
#'   Default is \code{TRUE}.
#' @param envir [\code{environment}]\cr
#'   This parameter is required for calling R functions from C.
#'   Do not change the default unless you know what you are doing!
#' @return [\code{List}] List with the following entries:
#' \describe{
#'   \item{dec}{Matrix of points in decision space (each column is a point).}
#'   \item{obj}{Matrix of objective values in objective space
#'   (each column is a point). Note that this is a matrix even in the
#'   single-objective case.}
#'   \item{history}{A list of named lists. The i-th component contains the
#'   \dQuote{dec} and \dQuote{obj} values (see preceding bullet points)
#'   if the population was logged in i-th generation or \code{NULL} if it
#'   was not logged.}
#' }
#'
#' @examples
#' library(smoof)
#'
#' # Single-Objective Example (see reference [2], Section 4.2)
#' # ===
#'
#' set.seed(1) # reproducibility
#' fn = smoof::makeHimmelblauFunction()
#' res = omniopt(fn, 100, 200, var.space.niching = TRUE, delta = 0.001,
#'   verbose = FALSE)
#'
#' \dontrun{
#' plot(fn)
#' points(t(res$dec))
#' }
#'
#' # Multi-Objective Example
#' # ===
#'
#' set.seed(1) # reproducibility
#' fn = smoof::makeZDT2Function(dimension = 4L)
#' res = omniopt(fn, 100, 1000, p.cross = 0.9, verbose = FALSE)
#' \dontrun{
#' plot(t(res$obj))
#' pairs(t(res$dec))
#' }
#' @export
omniopt = function(
  fn,
  pop.size = 4L,
  n.gens = 100L,
  p.cross = 0.6,
  p.mut = 1/smoof::getNumberOfParameters(fn),
  p.cross.bin = 0.6,
  p.mut.bin = 1/smoof::getNumberOfParameters(fn),
  eta.cross = 20L,
  eta.mut = 20L,
  mate = "normal",
  delta = 0.001,
  var.space.niching = FALSE,
  obj.space.niching = TRUE,
  init = "random",
  frequency = 10L,
  seed = runif(1),
  verbose = TRUE,
  envir = environment()) {

  checkmate::assert_class(fn, "smoof_function")
  pop.size = checkmate::asInt(pop.size, lower = 4L)
  if (pop.size < 4 || (pop.size %% 4) != 0) {
    stop("[omniopt] pop.size needs to be at least 4 and a multiple of 4.\n")
  }
  n.gens = checkmate::asInt(n.gens, lower = 1L)
  checkmate::assert_number(p.cross, lower = 0.6, upper = 1)
  checkmate::assert_number(p.mut, lower = 0, upper = 1)
  checkmate::assert_number(eta.cross, lower = 5L, upper = 20L)
  checkmate::assert_number(eta.mut, lower = 5L, upper = 50L)

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

  frequency = checkmate::asInt(frequency, lower = 1L, upper = n.gens)
  if (frequency < 1 || frequency > n.gens) {
    stop(sprintf("Parameter frequency should be in the range 1 to n.gens (%i).", n.gens));
  }

  checkmate::assert_number(seed, lower = 0, upper = 1)
  checkmate::assert_flag(verbose)

  n.objectives = smoof::getNumberOfObjectives(fn)

  par.set = getParamSet(fn)
  types = ParamHelpers::getParamTypeCounts(par.set)
  if ((types$numericvector > 0L) & (types$integervector > 0L)) {
    stop("[omniopt] Functions with both integervector and numericvector parameter sets are
      currently not supported.")
  }

  dimension.real = if (types$numericvector > 0L) smoof::getNumberOfParameters(fn) else 0L
  dimension.bin = if (types$integervector > 0L) smoof::getNumberOfParameters(fn) else 0L

  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)

  rawres = .Call("omnioptC",
    fn,
    as.integer(n.objectives),
    as.integer(dimension.real),
    as.double(lower),
    as.double(upper),
    as.integer(dimension.bin),
    as.integer(rep(1L, dimension.bin)),
    as.double(p.cross.bin),
    as.double(p.mut.bin),
    as.integer(pop.size),
    as.integer(n.gens),
    as.double(p.cross),
    as.double(p.mut),
    as.double(eta.cross),
    as.double(eta.mut),
    as.integer(mate),
    as.double(delta),
    as.integer(var.space.niching),
    as.integer(obj.space.niching),
    as.integer(init),
    as.integer(frequency),
    as.double(seed),
    as.integer(verbose),
    envir)

  # convert in "ecr"-style format
  names(rawres) = c("dec", "obj", "history")
  rawres$obj = t(matrix(rawres$obj, byrow = TRUE, ncol = n.objectives))
  rawres$dec = t(matrix(rawres$dec, byrow = TRUE, ncol = dimension))
  rawres$history = lapply(rawres$history, function(e) {
    if (is.null(e))
      return(NULL)
    list(
      dec = t(matrix(e[[1]], byrow = TRUE, ncol = dimension)),
      obj = t(matrix(e[[2]], byrow = TRUE, ncol = n.objectives))
    )
  })
  return(rawres)
}
