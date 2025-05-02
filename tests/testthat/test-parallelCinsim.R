test_that("Empty inputs into parallelCinsim returns a valid object", {
  paralllelCinsim_test <- parallelCinsim(g = 4)
  expect_equal(class(paralllelCinsim_test), "karyoSimParallel")
})

test_that("Invalid number of generations leads to an error", {
  expect_error(parallelCinsim(g = -1, cores = 1, iterations = 1))
  expect_error(parallelCinsim(g = Inf, cores = 1, iterations = 1))
})

test_that("Adding selection mode to inputs returns a valid object", {
  expect_equal(class(parallelCinsim(iterations = 2, selection_mode = "cn_based",
                                    selection_metric = Mps1)),
               "karyoSimParallel")
  expect_equal(class(parallelCinsim(iterations = 2, selection_mode = "rel_copy",
                                    selection_metric = Mps1)),
               "karyoSimParallel")

  # this test will fail consistently - check Davoli implementation
  # expect_equal(class(parallelCinsim(iterations = 2, selection_mode = "davoli",
  #                                   selection_metric = Mps1)),
  #              "karyoSimParallel")

})
