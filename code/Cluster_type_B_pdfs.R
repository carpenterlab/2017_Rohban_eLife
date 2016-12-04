rm(list = ls())

library(knitr)

for (i in 1:25) {
  cl.ind <- i
  knit("Cluster_Interpretation.Rmd")
}