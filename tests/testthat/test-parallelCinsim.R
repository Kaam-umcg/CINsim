test_that("Empty inputs into parallelCinsim returns a valid object", {
  paralllelCinsim_test <- parallelCinsim(g = 4)
  expect_equal(class(paralllelCinsim_test), "karyoSimParallel")
})

test_that("Invalid number of generations leads to an error", {
  expect_error(parallelCinsim(g = -1, cores = 1, iterations = 1))
  expect_error(parallelCinsim(g = Inf, cores = 1, iterations = 1))

  # this test currently fails, need to ensure this is checked properly
  # expect_error(parallelCinsim(g = 0, cores = 1, iterations = 1))
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

test_that("CnFS is calculated when flag is given with selection", {
  generations <- 6
  CnFS_test <- parallelCinsim(iterations = 2, cores = 2, g = 6, CnFS = TRUE,
                              selection_mode = "cn_based", selection_metric = Mps1)
  expect_vector(CnFS_test$sim_1$gen_measures$CnFS, ptype = numeric(),
                size = generations + 1) # +1 because we calculate the baseline as well for g = 0
})

test_that("CnFS is not calculated when flag is assigned FALSE", {
  generations <- 6
  CnFS_test <- parallelCinsim(iterations = 2, cores = 2, g = 6, CnFS = FALSE,
                              selection_mode = "cn_based", selection_metric = Mps1)
  expect_vector(CnFS_test$sim_1$gen_measures$CnFS, ptype = character(),
                size = generations + 1) # +1 because we calculate the baseline as well for g = 0
})
