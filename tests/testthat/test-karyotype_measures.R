test_that("Heterogeneity gets 0 in euploid sample", {
  n_chrom <- sample(20:30, 1)
  n_cells <- sample(50:100, 1)

  expect_equal(qHeterogeneity(matrix(rep(2, n_cells * n_chrom), ncol = n_chrom)),
               rep(0, n_chrom))
})

# test_that("Heterogeneity is 0.8 if there are 5 states with 0.2 frequency each", {
#   # add test case here
# })
