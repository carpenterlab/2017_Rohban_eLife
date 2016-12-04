rm(list = ls())
library(dplyr)

## finding out common genes up or downregulate as a result of following genes overexpression
gene.list <- c("TRAF2", "REL")

## threshold for repeatability of a landmark gene significance for each of the gene overexpressions mentioned above
thrs <- c(2, 2) 

## number of top/bottom landmark genes considered (based on z-scores)
rank.thr <- 50

## indicating whether top (positive) and/or bottom (negative) landmark genes should be considered
top <- T
bottom <- T
verbose <- F

t.comm <- c()
b.comm <- c()
itr <- 0

for (over.expressed.gene in gene.list) {
itr <- itr + 1
target.gene <- over.expressed.gene

thr <- thrs[itr]

ds.ta <- readRDS("../input/TA-OE-L1000-B1/TA.OE005_U2OS_72H_ZSPCQNORM_n1092x978.rds")

ds.ta.data <- attr(ds.ta, "mat")
ds.ta.cdesc <- attr(ds.ta, "cdesc")
ds.ta.rdesc <- attr(ds.ta, "rdesc")
gene.matching <- 0

lm.genes <- ds.ta.rdesc$pr_gene_symbol
ta.genes <- ds.ta.cdesc$x_genesymbol_mutation
lm.ta.genes <- intersect(lm.genes, ta.genes)

lm.ta.genes <- 1
verbose <- T
gene <- over.expressed.gene

for (over.expressed.gene in lm.ta.genes) {
  over.expressed.gene <- gene
  target.gene <- gene
  
  i <- which(ds.ta.cdesc$x_genesymbol_mutation %in% over.expressed.gene)
  tot.samp <- length(i)
  top.genes.pool <- c()
  bottom.genes.pool <- c()
  target.expression <- c()
  
  for (j in i) {
    yap.data <- ds.ta.data[,j]  
    sr <- sort(yap.data, decreasing = T)
    top.gene.loc <- names(sr[1:rank.thr])
    sr <- sort(yap.data, decreasing = F)
    bottom.gene.loc <- names(sr[1:rank.thr])
    
    top.genes.pool <- c(top.genes.pool, ds.ta.rdesc$pr_gene_symbol[which(ds.ta.rdesc$id %in% top.gene.loc)])
    bottom.genes.pool <- c(bottom.genes.pool, ds.ta.rdesc$pr_gene_symbol[which(ds.ta.rdesc$id %in% bottom.gene.loc)])
    
    target.expression <- c(target.expression, yap.data[which(ds.ta.rdesc$pr_gene_symbol == target.gene)])
  }
  
  top.genes.pool <- top.genes.pool[which(top.genes.pool != "-666")]
  bottom.genes.pool <- bottom.genes.pool[which(bottom.genes.pool != "-666")]
  
  tp.unique <- unique(top.genes.pool)
  vl <- rep(0, length(tp.unique)) 
  mrna.levels <- rep(0, length(tp.unique)) 
  is <- i
  i <- 1
  
  for (gn in tp.unique) {
    vl[i] <- which(top.genes.pool == gn) %>% length
    mrna.levels[i] <- ds.ta.data[ds.ta.rdesc$pr_gene_symbol == gn, is] %>% mean
    i <- i + 1
  }
  
  if (top) {
    indx <- which(vl >= thr)
    tp.data <- data.frame(gene = tp.unique[indx], repeatability = vl[indx]/tot.samp, average.mRNA.level = mrna.levels[indx]) %>% dplyr::mutate(repeatability = round(repeatability, 2), average.mRNA.level = round(average.mRNA.level, 2)) %>% dplyr::arrange(-repeatability, -average.mRNA.level)
    tp <- tp.data %>% htmlTable::htmlTable(.)
  }
  
  indx <- which(vl >= thr)
  df <- data.frame(gene = tp.unique[indx], repeatability = vl[indx]/tot.samp)
  rep.gene <- df[which(df$gene == over.expressed.gene), 2]
  if (length(rep.gene) != 0 && rep.gene >= 0.2) {
    print(sprintf("%s %f", over.expressed.gene, rep.gene))
    gene.matching <- gene.matching + 1
  }
  
  bt.unique <- unique(bottom.genes.pool)
  vl <- rep(0, length(bt.unique)) 
  i <- 1
  
  for (gn in bt.unique) {
    vl[i] <- which(bottom.genes.pool == gn) %>% length
    mrna.levels[i] <- ds.ta.data[ds.ta.rdesc$pr_gene_symbol == gn, is] %>% mean
    i <- i + 1
  }
  
  if (bottom) {
    indx <- which(vl >= thr)
    b.data <- data.frame(gene = bt.unique[indx], repeatability = vl[indx]/tot.samp, average.mRNA.level = mrna.levels[indx]) %>% dplyr::mutate(repeatability = round(repeatability, 2), average.mRNA.level = round(average.mRNA.level, 2)) %>% dplyr::arrange(-repeatability, average.mRNA.level)
    b <- b.data %>% htmlTable::htmlTable(.)
  }
}

if (top) {
  tp
}

if (bottom) {
  b
}

if (itr == 1) {
  if (top) {
    t.comm <- tp.data$gene  
  }
  
  if (bottom) {
    b.comm <- b.data$gene  
  }
} else {
  if (top) {
    t.comm <- intersect(tp.data$gene, t.comm)  
  }
  if (bottom) {
    b.comm <- intersect(b.data$gene, b.comm)  
  }
}

}

cat("top common")
print(t.comm)
cat("bottom common")
print(b.comm)