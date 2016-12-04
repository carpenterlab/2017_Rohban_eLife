rm(list = ls())

library("knitr")

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

library("reshape2")
library("dplyr")
library("stringr")
library("htmlTable")

Pf.trt <- Pf.collapsed
strong.trt <- Pf.trt.strong.collapsed$data$Treatment
data.annot <- Pf.strong %>% dplyr::filter(Treatment %in% strong.trt) %>% dplyr::select(one_of(c("Treatment", "Pathway", "Gene", "AlleleDesc"))) %>% 
  dplyr::mutate(Pathway = str_replace(Pathway, "Canonical ", "")) %>% 
  dplyr::mutate(AlleleDesc = str_sub(AlleleDesc, 1, 2)) %>%
  dplyr::filter(AlleleDesc == "WT") %>% dplyr::group_by(Gene) %>% dplyr::slice(1) %>%
  dplyr::ungroup()

strong.trt <- data.annot$Treatment
data <- Pf.trt$data %>% dplyr::filter(Treatment %in% strong.trt) %>% dplyr::select(one_of(c("Treatment", Pf.trt$feat_cols))) 
rownames(data) <- data$Treatment
data <- data[,2:NCOL(data)]
cr <- cor(data %>% t)
cr.melt <- cr %>% melt
t1 <- merge(cr.melt, data.annot, by.x = "Var1", by.y = "Treatment")
t2 <- merge(t1, data.annot, by.x = "Var2", by.y = "Treatment")

thr <- 0.4333212
v1 <- t2 %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & 
                       Pathway.x == Pathway.y & (value) >= thr) 
v10 <- t2 %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & 
                              Pathway.x == Pathway.y & (value) < thr)
v2 <- t2 %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & 
                             Pathway.x != Pathway.y & (value) >= thr) 
v20 <- t2 %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & 
                             Pathway.x != Pathway.y & (value) < thr)
contingency.table <- rbind(c(NROW(v1), NROW(v10)), c(NROW(v2), NROW(v20))) %>% t
print(contingency.table)
fisher.test(contingency.table, alternative = "greater") %>% print

v1.subset <- v1 %>% dplyr::select(one_of(c("Var1", "Var2", "Pathway.x"))) 
colnames(v1.subset) <- c("gene 1", "gene 2", "common pathway")
v1.subset %>% htmlTable()
