test_that("CnFS correctly errors out on wrong inputs", {
  test_karyo <- makeKaryotypes(numCell = 20, species = "mouse")
  expect_error(calc_CnFS(karyotypes = test_karyo))
})

test_that("calc_CnFS returns a correct value for a single euploid cell", {
  test_karyo_mouse <- makeKaryotypes(numCell = 1, species = "mouse")
  test_karyo_human <- makeKaryotypes(numCell = 1, species = "human")
  expect_equal(
    calc_CnFS(
      karyotypes = test_karyo_mouse,
      selection_metric = CINsim::Mps1
    ),
    0.07890889,
    tolerance = 1e-4
  )

  expect_equal(
    calc_CnFS(
      karyotypes = test_karyo_human,
      selection_metric = CINsim::human_neutral
    ),
    0.3636364,
    tolerance = 1e-4
  )
})

test_that("calc_KMS returns a correct value for a single euploid cell", {
  test_karyo_mouse <- makeKaryotypes(numCell = 1, species = "mouse")
  test_karyo_human <- makeKaryotypes(numCell = 1, species = "human")

  expect_equal(
    calc_KMS(
      karyotypes = test_karyo_mouse,
      pop_measures = CINsim::karyotype_measures
    ),
    0.08493877
  )

  # this test needs pop measures to compare against - included those
  # in the fixtures? Can't really test without having a constant
  # to calculate the KMS against

  # expect_equal(
  #   calc_KMS(
  #     karyotypes = test_karyo_human,
  #     pop_measures = CINsim::human_neutral
  #   ),
  # )
})

test_that("calc_CnFS returns a correct value for multiple euploid cells", {
  test_karyo_mouse <- makeKaryotypes(numCell = sample.int(1e4, 1),
                                     species = "mouse")
  test_karyo_human <- makeKaryotypes(numCell = sample.int(1e4, 1),
                                     species = "human")
  expect_equal(
    calc_CnFS(
      karyotypes = test_karyo_mouse,
      selection_metric = CINsim::Mps1
    ),
    0.07890889,
    tolerance = 1e-4
  )

  expect_equal(
    calc_CnFS(
      karyotypes = test_karyo_human,
      selection_metric = CINsim::human_neutral
    ),
    0.3636364,
    tolerance = 1e-4
  )
})

test_that("calc_KMS returns a correct value for multiple euploid cells", {
  test_karyo_mouse <- makeKaryotypes(numCell = sample.int(1e4, 1),
                                     species = "mouse")
  test_karyo_human <- makeKaryotypes(numCell = sample.int(1e4, 1),
                                     species = "human")

  expect_equal(
    calc_KMS(
      karyotypes = test_karyo_mouse,
      pop_measures = CINsim::karyotype_measures
    ),
    0.08493877
  )

  # this test needs pop measures to compare against - included those
  # in the fixtures? Can't really test without having a constant
  # to calculate the KMS against

  # expect_equal(
  #   calc_KMS(
  #     karyotypes = test_karyo_human,
  #     pop_measures = CINsim::human_neutral
  #   ),
  # )
})

