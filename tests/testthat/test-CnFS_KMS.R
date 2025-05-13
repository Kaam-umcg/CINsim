test_that("CnFS correctly errors out on wrong inputs", {
  test_karyo <- makeKaryotypes(numCell = 20, species = "mouse")
  expect_error(calc_CnFS(karyotypes = test_karyo))
})

test_that("calc_CnFS returns a correct value for euploid cells", {
  test_karyos_mouse <- makeKaryotypes(numCell = 1, species = "mouse")
  test_karyos_human <- makeKaryotypes(numCell = 1, species = "human")
  expect_equal(
    calc_CnFS(
      karyotypes = test_karyos_mouse,
      selection_metric = CINsim::Mps1
    ),
    0.07890889
  )

  expect_equal(
    calc_CnFS(
      karyotypes = test_karyos_human,
      selection_metric = CINsim::human_neutral
    ),
    0.3636364
  )
})

test_that("calc_KMS returns a correct value for euploid cells", {

})
