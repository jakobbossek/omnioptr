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

    for (obj.space.niching in c(TRUE, FALSE)) {
      for (var.space.niching in c(TRUE, FALSE)) {
        # at least one of these must be TRUE
        if (!(var.space.niching || obj.space.niching))
          next
        for (mate in c("normal", "restricted")) {
          res = omniopt(f, pop.size = pop.size, n.gens = 100L,
            obj.space.niching = obj.space.niching,
            var.space.niching = var.space.niching,
            mate = mate,
            verbose = FALSE)
          checkmate::expect_list(res)
          checkmate::expect_matrix(res$dec, nrows = d, ncols = pop.size, any.missing = FALSE, all.missing = FALSE)
          checkmate::expect_matrix(res$obj, nrows = o, ncols = pop.size, any.missing = FALSE, all.missing = FALSE)

          # check whether fitness values of individuals match the smoof-output
          expect_true(all(res$obj == apply(res$dec, 2L, f)))

          # check if individuals respect box-constraints
          lower = smoof::getLowerBoxConstraints(f)
          upper = smoof::getUpperBoxConstraints(f)

          expect_true(all(apply(res$dec, 2L, function(e) { all((e >= lower) & (e <= upper))})))
        }
      }
    }
  }
})
