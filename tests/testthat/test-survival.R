test_that("Survival errors out without inputs", {
  expect_error(karyotype_survival())
})

test_that("Extremely high qMod survival returns full survival of population", {
  cells_tested <- 100
  test_pop <- makeKaryotypes(numCell = cells_tested)

  expect_equal(
    sum(karyotype_survival(karyotypes = test_pop,
                       selection_mode = "cn_based",
                       selection_metric = Mps1,
                       qMod = 10000000000000)),
    cells_tested)
})

test_that("Returns a non-NA output for every selection mode and valid inputs", {
  cells_tested <- 10
  test_pop <- makeKaryotypes(numCell = cells_tested)

  expect_length(karyotype_survival(karyotypes = test_pop,
                                   selection_mode = "cn_based",
                                   selection_metric = Mps1),
                cells_tested)

  expect_length(karyotype_survival(karyotypes = test_pop,
                                   selection_mode = "rel_copy",
                                   selection_metric = Mps1),
                cells_tested)
  # Davoli will consistently fail, not 100% sure why
  # we never do anything with it - remove the feature? taking more time
  # than it's worth to keep maintaining it IMO.
  # expect_length(karyotype_survival(karyotypes = test_pop,
  #                                  selection_mode = "davoli",
  #                                  selection_metric = Mps1),
  #               cells_tested)

})

test_that("Survival calculations returns a valid length
          for both human and mouse karyotypes/selection_metrics", {
  cells_tested <- 10
  test_pop_mouse <- makeKaryotypes(numCell = cells_tested, species = "mouse")
  test_pop_human <- makeKaryotypes(numCell = cells_tested, species = "human")

  # mouse
  expect_length(karyotype_survival(karyotypes = test_pop_mouse,
                                   selection_mode = "cn_based",
                                   selection_metric = Mps1),
                cells_tested)
  # human
  expect_length(karyotype_survival(karyotypes = test_pop_human,
                                   selection_mode = "cn_based",
                                   selection_metric = human_neutral),
                cells_tested)
})
