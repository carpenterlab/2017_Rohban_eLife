context("Preprocess integrity")

test_that("PCA gives meaningful results", {
  P <- profile.data.dummy()
  P <- pca(P, .80)
  expect_equal(NCOL(feats(P)),7)
})