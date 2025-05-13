test_that("Can calculate pop measures for human and mouse single karyotype
          at generation 0", {
  test_karyo_human <- makeKaryotypes(numCell = 1, species = "human")
  test_karyo_mouse <- makeKaryotypes(numCell = 1, species = "mouse")

  expect_equal(dim(population_measures(
    karyotypes = test_karyo_mouse,
    euploid_ref = 2,
    g = 0)),
    c(length(test_karyo_mouse), 7) # should be 7 entries in pop measures
  )
  expect_equal(dim(population_measures(
    karyotypes = test_karyo_human,
    euploid_ref = 2,
    g = 0)),
    c(length(test_karyo_human), 7) # should be 7 entries in pop measures
  )
})

test_that("Can calculate pop measures for human and mouse multiple karyotype
          at generation 0", {
  test_karyo_human <- makeKaryotypes(numCell = sample.int(5000, 1),
                                     species = "human")
  test_karyo_mouse <- makeKaryotypes(numCell = sample.int(5000, 1),
                                     species = "mouse")
  expect_equal(dim(population_measures(
    karyotypes = test_karyo_mouse,
    euploid_ref = 2,
    g = 0)),
    c(dim(test_karyo_mouse)[2], 7)) # should be 7 entries in pop measures
    # and 1 row per chromosome

  expect_equal(dim(population_measures(
    karyotypes = test_karyo_human,
    euploid_ref = 2,
    g = 0)),
    c(dim(test_karyo_human)[2], 7) # should be 7 entries in pop measures
    # and 1 row per chromosome
  )
})
