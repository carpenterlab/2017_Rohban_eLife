rm(list = ls())

library("knitr")
library(dplyr)

## loading initial analysis
use.cache <- T
workspace.file <- "../results/master/Initial_analysis/Initial_analysis_workspace.RData"

if (use.cache) {
  if (file.exists(workspace.file)) {
    load(workspace.file)
  } else {
    knit("Initial_analysis.Rmd")
    load(workspace.file)
  }
} else {
  knit("Initial_analysis.Rmd")
  load(workspace.file)
}

g1 <- c("BRAF_WT.1", "RAF1_WT.1", "KRAS_WT.1", "CDC42_WT", "AKT3_WT.2", "AKT1_WT.1", "RHOA_WT", "RAC1_WT.1")
g2 <- c("BRAF_V600E", "RAF1_L613V", "KRAS_G12V", "CDC42_Q61L", "AKT3_E17K", "AKT1_E17K", "RHOA_Q63L", "RAC1_Q61L")

x1 <- Pf.trt.trt$data[which(Pf.trt.trt$data$Treatment %in% g1), c("Treatment", "sim_pearson_q50")] %>% dplyr::arrange(Treatment)
x2 <- Pf.trt.trt$data[which(Pf.trt.trt$data$Treatment %in% g2), c("Treatment", "sim_pearson_q50")] %>% dplyr::arrange(Treatment)

t.test(x2$sim_pearson_q50, x1$sim_pearson_q50, paired = T, alternative = "greater")

