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

library(dplyr)
library(stringr)

get.interacting.proteins <- function(protein.name) {
  tbl <- read.table(sprintf("http://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=%s&taxId=9606&includeHeader=true&accesskey=ca950b072394ce1897811022f7757222", protein.name), sep="\t", header = FALSE, fill = T)
  tbl <- tbl[, c(8, 9, 12, 13, 14)] 
  colnames(tbl) <- c("Protein.1", "Protein.2", "Method", "Type", "Evidence")
  return(tbl)
}

get.all.interacting.proteins <- function(protein.names) {
  tmp <- lapply(protein.names, get.interacting.proteins)
  tmp.bind <- do.call(rbind, tmp) %>% unique
  tmp.bind <- tmp.bind %>% dplyr::filter(Protein.1 %in% protein.names & Protein.2 %in% protein.names & (Protein.1 %>% as.character()) != (Protein.2 %>% as.character()))
}

Pf <- Pf.trt.strong.collapsed$data

genes <- Pf$Treatment %>% lapply(., function(x) str_split(x, "_")[[1]][1]) %>% unlist %>% unique
pr.pr.interaction <- get.all.interacting.proteins(genes)
saveRDS(pr.pr.interaction, "../results/master/protein_interaction_data/pr.pr.interaction.rds")
write.csv(pr.pr.interaction, "../results/master/protein_interaction_data/pr.pr.interaction.csv")

pr.pr.binary <- outer(rep(0, length(genes)), rep(0, length(genes)), "*")
rownames(pr.pr.binary) <- genes
colnames(pr.pr.binary) <- genes
for (i in 1:NROW(pr.pr.interaction)) {
  g1 <- pr.pr.interaction[i, 1] %>% as.character()
  g2 <- pr.pr.interaction[i, 2] %>% as.character()
  pr.pr.binary[g1, g2] <- 1
  pr.pr.binary[g2, g1] <- 1
}

saveRDS(pr.pr.binary, "../results/master/protein_interaction_data/pr.pr.binary.rds")
