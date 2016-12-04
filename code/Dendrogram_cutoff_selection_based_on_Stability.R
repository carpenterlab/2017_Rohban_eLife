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
library(forecast)
library(ggplot2)
library(rwantshue)
library(ape)
library(randomcoloR)

set.seed(42)

## step size for stability analysis
scale <- 0.002

Pf.trt <- Pf.trt 
Pf.trt.collapsed <- Pf.trt.strong.collapsed 
Pf.trt.trt <- Pf.trt.trt 
strongs.indx.General <- Pf.trt.collapsed$data$Treatment %>% unique
Pf.trt$data <- Pf.trt$data %>% dplyr::filter(Treatment %in% strongs.indx.General)

Px <- Pf.trt.collapsed$data[,Pf.trt.collapsed$feat_cols]
rownames(Px) <- Pf.trt.collapsed$data$Treatment
cr <- Px %>% t %>% cor
hcl <- hclust((1 - cr) %>% as.dist(), method = "average")
hcl0 <- hcl

find.clusters <- function(ct, min.size = 1) {
  cls <- c()
  for (i in 1:max(ct)) {
    nm <- names(which(ct == i))
    if (length(nm) < min.size) {
      next
    }
    cls <- c(cls, list(nm))
  }
  return(cls)
}

compare.clusters <- function(cls1, cls2) {
  indx <- lapply(cls1, function(x) return(lapply(cls2, function(y) (identical((x %>% sort), (y %>% sort)))) %>% unlist %>% any)) %>% unlist 
  sm <- cls1[which(indx)] %>% unlist %>% length
  return(sm)
}

thrs <- seq(from = 0.4, to = 0.7, by = scale)
stab.score <- c()
num.clusters <- c()
pb <- progress::progress_bar$new(total = length(thrs) - 2)

for (i in 2:(length(thrs) - 1)) {
  ct <- cutree(hcl, h = 1 - thrs[i]);
  ct.p <- cutree(hcl, h = 1 - thrs[i-1]);
  ct.n <- cutree(hcl, h = 1 - thrs[i+1]);
  
  cls <- find.clusters(ct)
  cls.p <- find.clusters(ct.p)
  cls.n <- find.clusters(ct.n)
  v.n <- compare.clusters(cls, cls.n)
  v.p <- compare.clusters(cls, cls.p)
  stab.score <- c(stab.score, (v.n + v.p)/2)
  num.clusters <- c(num.clusters, length(cls))
  pb$tick()
}

stab.score.ma <- ma(stab.score/NROW(cr), order = 1/scale * 0.02)
D <- data.frame(thresholds = thrs[2:(length(thrs) - 1)], stability = stab.score.ma)
D.num <- data.frame(thresholds = thrs[2:(length(thrs) - 1)], num.clusters = num.clusters)
g <- ggplot(D, aes(x = thresholds, y = stability))
g <- g + geom_line() + scale_x_continuous(breaks = seq(0.4, 0.7, 0.05))
g <- g + xlab("Threshold to cut dendrogram") + ylab("Stability score")
ggsave("stability_clustering.pdf", g, width = 8, height = 6)

thr.sel <- D %>% dplyr::arrange(-stability) %>% head(., 1) %>% dplyr::select(thresholds) %>% as.matrix() %>% as.vector()

print(thr.sel)

hcl <- hcl0
ct <- cutree(hcl, h = 1 - thr.sel)
cls <- find.clusters(ct, 2)
print(cls)
saveRDS(cls, "cluster_names.rds")
cls <- find.clusters(ct, 1)

labs <- hcl$labels 
cols <- rep("", length(labs))
fn <- rep(0, length(labs))
cexs <- rep(1.0, length(labs))

ncl <- lapply(cls, function(x) length(x) > 1) %>% unlist %>% sum

colors <- distinctColorPalette(ncl)
i <- 1
for (ci in colors) {
  hi <- rgb2hsv(col2rgb(ci))
  hi[2] <- runif(1, 0.6, 1)
  colors[i] <- hsv(hi[1], hi[2], hi[3], 1)
  i <- i + 1
}

j <- 1
for (i in 1:length(cls)) {
  indx <- which(labs %in% cls[[i]])
  if (length(indx) > 1) {
    cols[indx] <- colors[j]
    fn[indx] <- 2
    cexs[indx] <- 1.0
    j <- j + 1
  } else {
    cols[indx] <- "gray55"
    fn[indx] <- 1
    lb <- hcl$labels[hcl$labels %in% cls[[i]]] 
    new.lb <- lapply(lb, function(x) paste("*", x, "*", collapse = "")) %>% unlist
    hcl$labels[hcl$labels %in% cls[[i]]] <- new.lb
  }
}

print("number of clusters of size greater than 1 : ")
print(j - 1)

graphics.off()
plot.new()
quartz(width = 15, height = 15)
plot(as.phylo(hcl), type = "fan", tip.color = cols, font = fn, cex = cexs)
dev.print(pdf, "dendro.pdf")