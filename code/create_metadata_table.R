library("knitr")

## loading initial analysis
use.cache <- T
workspace.file <- "../results/master/Initial_analysis/Initial_analysis_workspace.RData"
load(workspace.file)

library("readr")
library("plyr")
library("dplyr")
library("xtable")
library("corrplot")

df <- Pf_org.org$data %>% dplyr::select(one_of(c("Well", "Gene", "Treatment", "PUBLICID"))) %>% unique 

pid <- df$PUBLICID %>% lapply(., function(x) ifelse(x != 0, as.character(x), "NA")) %>% unlist
df$PUBLICID <- pid
df <- cbind(df, data.frame(cell_line = "U2OS"))
colnames(df) <- c("well_position",   "gene_name",       "pert_name",       "broad_sample",    "cell_line")

write.csv(df, "metadata.csv", row.names = F)