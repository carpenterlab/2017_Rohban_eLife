context("Loading profile.data")

cf = '../../data/HDAC-Knockdowns-Batch-3/hdac-ko-well-mean-global_norm.yml'
P <- profile.data(cf)

test_that("loading configuration file works", {
  expect_equal(P$cfg$feat_start, 4)
})

