test_that("Empty inputs returns a Cinsim object", {
  Cinsim_output <- Cinsim()
  expect_equal(class(Cinsim_output), "karyoSim")
})

test_that("Adding selection mode doesn't error out", {
  expect_no_error(Cinsim(selection_mode = "cn_based",
                         selection_metric = CINsim::Mps1))
  expect_no_error(Cinsim(selection_mode = "rel_copy",
                         selection_metric = CINsim::Mps1))

  # Davoli only works with human human genomes
  # or maybe not? Need to investigate a little more
  # before this just works
  # expect_no_error(Cinsim(karyotypes = makeKaryotypes(species = "human",
  #                                                    numCell = 10),
  #                        selection_mode = "davoli",
  #                        selection_metric = CINsim::human_neutral))
})

test_that("CnFS can be calculated when flag is given with selection", {
  generations <- 6
  Cinsim_output <- Cinsim(g = generations, CnFS = TRUE, selection_mode = "cn_based",
                          selection_metric = Mps1)
  expect_equal(class(Cinsim_output), "karyoSim")
  expect_vector(Cinsim_output$gen_measures$CnFS, ptype = numeric(),
                size = generations + 1) # +1 for CnFS of generation 0
})

test_that("Verbosity can be controlled with the designated flag", {
  test_karyos <- makeKaryotypes(numCell = 1, species = "mouse")
  # need to set some other flags to ensure we don't get other messages
  # maybe those should be errors instead?
  expect_no_message(
    Cinsim(karyotype = test_karyos, verbose = 0, CnFS = FALSE)
  )

  expect_message(Cinsim(karyotype = test_karyos, verbose = 1, CnFS = FALSE),
                 regexp = "Processing generation")
  expect_message(Cinsim(karyotype = test_karyos, verbose = 1, CnFS = FALSE),
                 regexp = "Processing generation")


})
