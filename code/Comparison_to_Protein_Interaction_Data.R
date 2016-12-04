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
library(reshape2)
library(stringr)
library(htmlTable)
library(ggplot2)

## use gather_protein_interaction_data.R to obtain updated versions of the following PPI data
ppi <- readRDS("../results/master/protein_interaction_data/pr.pr.binary.rds")
ppi.detailed <- readRDS("../results/master/protein_interaction_data/pr.pr.interaction.rds")

Pf <- Pf.trt.strong.collapsed

matching.per.table <- readRDS("../results/master/ORFs_sequence_matching_transcripts_percentage/matching.per.table.rds")
leg.orfs <- matching.per.table %>% dplyr::filter(Nuc..Match.. > 99 & Prot..Match.. > 99) %>% dplyr::select(one_of("Treatment")) %>% as.matrix() %>% as.vector()
Pf$data <- Pf$data %>% dplyr::filter(Treatment %in% leg.orfs)

d.trt <- Pf$data %>% dplyr::select(one_of(Pf$feat_cols))
rownames(d.trt) <- Pf$data$Treatment
cr <- d.trt %>% t %>% cor

new.col <- intersect(colnames(cr), leg.orfs)
cr <- cr[new.col, new.col]
cr.org <- cr

genes <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][1]) %>% unlist
alleles <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][2]) %>% unlist %>% 
  lapply(., function(x) str_split(x, "\\.")[[1]][1]) %>% unlist

u <- outer(genes, genes, function(x1, x2) return(x1 == x2))
v <- outer(alleles, alleles, function(x1, x2) return(x1 == "WT" & x2 == "WT"))
z <- outer(rownames(cr), rownames(cr), function(x, y) return(x < y))

w <- u & v & z
rownames(w) <- rownames(cr)
colnames(w) <- colnames(cr)
cr.melt <- melt(cr)
w.melt <- melt(w)
wt.cr.melt <- cr.melt[which(w.melt$value),]

trt.to.go <- c()
for (i in 1:length(genes)) {
  if (alleles[i] == "WT" && (length(trt.to.go) == 0 || !(genes[i] %in% genes[trt.to.go]))) {
    trt.to.go <- c(trt.to.go, i)
  } 
}
cr <- cr[trt.to.go, trt.to.go]

d1 <- density(cr[upper.tri(cr)], from = -1, to = 1, n = 512) 
d2 <- density(wt.cr.melt$value %>% as.matrix() %>% as.vector(), from = -1, to = 1, n = 512) 
d1$y <- d1$y * 0.5
d2$y <- d2$y * 0.5

d1.x <- cr[upper.tri(cr)]
d2.x <- wt.cr.melt$value %>% as.matrix() %>% as.vector()

D <- data.frame("correlation" = c(d1.x, d2.x), "group" = c(rep("all non-clone pairs", length(d1.x)), rep("wild type clone pairs", length(d2.x))))
g <- ggplot(D, aes_string(x = "correlation", fill = "group", color = "group")) + geom_density(alpha = 0.5) + scale_fill_discrete(name = "") + scale_color_discrete(guide = F) + xlab("Profile Correlation") + ylab("Density")
ggsave("all_corr_vs_WT.png", g, width = 9, heigh = 6)

ind <- which(d1$x > 0 & d1$x < 1)
cand.indx <- which.min(abs(d1$y[ind]-d2$y[ind]))
thr <- min(d1$x[ind[cand.indx]])

consistency.score <- function(cr, ppi, verbose = F, thr = 0.5, ppi.detailed = NULL) {
  cr.melt <- cr %>% melt
  cl <- colnames(cr.melt)
  cl[1] <- "Var1"
  cl[2] <- "Var2"
  colnames(cr.melt) <- cl
  
  cr.melt.subset <- cr.melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & (value) >= thr)
  
  if (verbose) {
    print(min(cr.melt.subset$value))
    print(sprintf("The number of pairs suggested : %d", NROW(cr.melt.subset)))
  }
  genes1 <- lapply(cr.melt.subset[,1] %>% as.matrix %>% as.vector, function(x) str_split(x, "_")[[1]][1]) %>% unlist
  genes2 <- lapply(cr.melt.subset[,2] %>% as.matrix %>% as.vector, function(x) str_split(x, "_")[[1]][1]) %>% unlist
  gene.pairs <- cbind(genes1, genes2) %>% unique
  
  gene.ref <- rownames(ppi)
  gene.pairs <- gene.pairs[which(gene.pairs[,1] %in% gene.ref & gene.pairs[,2] %in% gene.ref),]
  sm <- apply(gene.pairs, 1, function(x) return(max(ppi[x[1], x[2]], ppi[x[1], x[2]]))) %>% sum
  u <- c()
  if (verbose) {
    for (i in 1:NROW(gene.pairs)) {
      if (ppi[gene.pairs[i, 1], gene.pairs[i, 2]] == 1) {
        if (!is.null(ppi.detailed)) {
          v <- ppi.detailed %>% dplyr::filter((Protein.1 == gene.pairs[i, 1] & Protein.2 == gene.pairs[i, 2]) | (Protein.2 == gene.pairs[i, 1] & Protein.1 == gene.pairs[i, 2])) 
          u <- rbind(u, v)
        }
      }
    }
    
    if (!is.null(u)) {
      cut.str <- function(x) {
        x <- as.character(x)
        if (str_length(x) > 100) {
          i <- str_locate(x, "\\)")
          return(str_sub(x, 1, i[1]))
        } else {
          return(x)
        }
      }
      
      u$Evidence <- lapply(u$Evidence, cut.str) %>% unlist
      u %>% htmlTable() %>% cat(., file = "evidence.table.html")
    }
  }
  return(sm)
}


sig <- consistency.score(cr, ppi, T, thr, ppi.detailed = ppi.detailed)

cr.melt <- cr %>% melt
cl <- colnames(cr.melt)
cl[1] <- "Var1"
cl[2] <- "Var2"
colnames(cr.melt) <- cl

cr.melt.subset <- cr.melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & (value) >= thr)

genes <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][1]) %>% unlist
genes.in <- intersect(colnames(ppi), genes)
all.verified <- (ppi[genes.in, genes.in] %>% sum())/2
x2 <- sig
x1 <- NROW(cr.melt.subset) - x2

data <- (rbind(c(x1, x2), c(NROW(cr) * (NROW(cr) - 1)/2 -all.verified-x1, all.verified-x2)))
fisher.test(data, alternative = "less") %>% print
data %>% print
