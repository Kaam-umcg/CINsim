# this script creates some objects for use in the testing suite,
# specifically for plotting functions. They are labelled according
# to their class, which is the most relevant thing for error
# reporting and how functions use their control flow

set.seed(42)
simple_cinsim_test_object <- Cinsim_simple()
saveRDS(simple_cinsim_test_object, "tests/testthat/fixtures/simple_cinsim_object.Rds")

set.seed(43)
cinsim_test_object <- Cinsim()
saveRDS(cinsim_test_object, "tests/testthat/fixtures/cinsim_object.Rds")

set.seed(44)
parallel_cinsim_test_object <- parallelCinsim()
saveRDS(parallel_cinsim_test_object, "tests/testthat/fixtures/parallel_cinsim_object.Rds")
