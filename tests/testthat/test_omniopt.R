context("OmniOptimizer C Interface")

test_that("omnioptr::omniopt produces correct outpur", {
  # test functions
  funs = c(
    smoof::makeZDT1Function(dimension = 2L),
    smoof::makeZDT2Function(dimension = 4L)
  )
  pop.size = 100L

  for (f in funs) {
    d = smoof::getNumberOfParameters(f)
    o = smoof::getNumberOfObjectives(f)
    res = omniopt(f, pop.size = pop.size, n.gens = 200L)
    checkmate::expect_list(res)
    checkmate::expect_matrix(res$pareto.set, nrows = d, ncols = pop.size, any.missing = FALSE, all.missing = FALSE)
    checkmate::expect_matrix(res$pareto.front, nrows = o, ncols = pop.size, any.missing = FALSE, all.missing = FALSE)

    # check whether fitness values of individuals match the smoof-output
    expect_true(all(res$pareto.front == apply(res$pareto.set, 2L, f)))

    # check if individuals respect box-constraints
    lower = smoof::getLowerBoxConstraints(f)
    upper = smoof::getUpperBoxConstraints(f)

    expect_true(all(apply(res$pareto.set, 2L, function(e) { all((e >= lower) & (e <= upper))})))
  }
})
