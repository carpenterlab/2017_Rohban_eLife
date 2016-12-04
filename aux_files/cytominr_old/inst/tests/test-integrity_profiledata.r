context("profile.data integrity")
library(stringr)

cf = '../../data/HDAC-Knockdowns-Batch-3/hdac-ko-well-mean-global_norm.yml'
P <- profile.data(cf)
P <- cleanup(P)
test_that("the features returned by feats() are valid", {
  expect_true(all((unique(unlist(llply(names(feats(P)), 
                                       function(x) unlist(str_split(x, "_"))[1]))) %in% 
                                         c("Cytoplasm","Cells","Nuclei" ))))
  })


