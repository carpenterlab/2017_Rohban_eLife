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
library(topGO)

if (!require("org.Hs.eg.db", quietly = T)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}

Pf.trt <- Pf.trt 
Pf.trt.collapsed <- Pf.trt.strong.collapsed 
strongs.indx.General <- Pf.trt.collapsed$data$Treatment %>% unique

min.cluster.size <- 2
cor_type <- "pearson"
clustering_method <- "average"

## used "Dendrogram_cutoff_selection_based_on_Stability.R" to figure out the following threshold
depth.to.cut <- 1 - 0.522

cluster.size.cut.off.cls <- min.cluster.size

df.to.use <- Pf.trt.collapsed
Px <- data.frame(df.to.use$data %>% dplyr::filter(Treatment %in% strongs.indx.General) %>% dplyr::select(starts_with("PC")))
trts.Px <- as.character((df.to.use$data %>% dplyr::filter(Treatment %in% strongs.indx.General))$Treatment)
pr <- permute::shuffle(length(trts.Px))
rownames(Px) <- trts.Px

cr <- cor(t(Px), method = cor_type)

saveRDS(cr, "cr.rds")
dsts <- as.dist(1 - cr)
hcl_0 <- hclust(d=dsts, method = clustering_method)
ct <- cutree(hcl_0, h = depth.to.cut)
ct.cls <- ct

cls <- c()
cls.gene <- c()
cl.num <- max(ct)

for (i in 1:cl.num) {
  trts <- names(which(ct == i))
  if (length(trts) >= min.cluster.size) {
    cls <- c(cls, list(trts))
    gn <- lapply(trts, function(x) return(str_split(x, "_")[[1]][1])) %>% unlist %>% unique
    if (length(gn) >= min.cluster.size) {
      cls.gene <- c(cls.gene, list(gn))
    }
  }
}

genes <- cls.gene %>% unlist %>% unique
cls.gene.int <- rep(0, length(genes))
names(cls.gene.int) <- genes

for (g in genes) {
  ind <- lapply(cls.gene, function(x) (g %in% x)) %>% unlist 
  cls.gene.int[g] <- which(ind)[1]
}

gene.id.map <- Pf.trt$data[,c("TA_GeneID", "Gene")] %>% dplyr::filter(Gene %in% genes) %>% unique
rownames(gene.id.map) <- gene.id.map$Gene

cls.gene.int.new <- cls.gene.int
names(cls.gene.int.new) <- gene.id.map[genes, "TA_GeneID"] %>% as.character()

saveRDS(cls, "cls.rds")
print(cls)

for (i in 1:length(cls)) {
  gns <- lapply(cls[[i]], function(x) (str_split(x, "_")[[1]][1])) %>% unlist %>% unique
  if (length(gns) <= 1) {
    next
  }
  
  v <- 1:NROW(gene.id.map)
  names(v) <- gene.id.map$TA_GeneID
  
  sampleGOdata <- new("topGOdata", description = "Simple session",
                      ontology = "BP", allGenes = v, 
                      geneSel = function(x) (gene.id.map[x, "Gene"] %in% gns), 
                      nodeSize = 2, annot = annFUN.org, mapping = "org.Hs.eg.db")
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  w <- GenTable(sampleGOdata, classicFisher = resultFisher, topNodes = 5, numChar = 100)
  w <- w %>% dplyr::filter(Significant >= 2) %>% dplyr::arrange(classicFisher) %>% dplyr::mutate(sig = paste(" ", Significant, "/", Annotated, sep = ""))
  print(gns)
  print(w)
  print("------------------")
  write.csv(w, sprintf("%d_clust_topGO.csv", i))
}
