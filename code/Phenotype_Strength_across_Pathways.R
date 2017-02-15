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
library(ggplot2)

Pf.trt <- Pf.trt 
Pf.trt.collapsed <- Pf.trt.strong.collapsed 

null.data <- nulls.samp.rep.corr
hit.sel.thr <- quantile(null.data, 0.95) %>% as.matrix() %>% as.vector()

cor_type <- "pearson"
var.name <- "sim_pearson_q50"
strongs.indx <- strongs.indx.General
strongs.indx.General <- Pf.trt$data$Treatment %>% unique

Pf.trt.tmp <- Pf.trt
Pf.trt$data$Pathway <- as.character(Pf.trt$data$Pathway)
Pf.trt$data$Pathway <- lapply(Pf.trt$data$Pathway, function(x) str_replace_all(x, "Canonical ", "")) %>% unlist
Pf.trt$data$Pathway[which(Pf.trt$data$Pathway == "Transcription Factors")] <- "C-EBP alpha"
pthws <- unique(Pf.trt$data$Pathway)
pthws.all <- pthws

n <- 0
p.name <- c()
rep.values <- c()
rep.sd.values <- c()
cnts.n <- c()
ns <- c()
p.type.name <- "All"
trts.all <- c()

for (pth in pthws.all) {
  trts <- unique(Pf.trt$data$Treatment[which(Pf.trt$data$Treatment %in% strongs.indx.General & (Pf.trt$data$Pathway == pth) & !str_detect(Pf.trt$data$Treatment, "mismatch"))])
  trts.tm <- c()
  for (trt in trts) {
    if (!str_detect(str_split(trt, "_")[[1]][2], "WT")) {
      ii <- c()
    } else {
      gn <- str_split(trt, "_")[[1]][1]
      ii <- which(str_detect(trts.tm, paste(gn, "_", sep = "")) & trts.tm != trt)  
    }
    
    if (length(ii) == 0) {
      trts.tm <- c(trts.tm, trt)
    } else {
      fl <- TRUE
      for (jj in ii) {
        if (str_detect(str_split(trts.tm[jj], "_")[[1]][2], "WT")) {
          fl <- FALSE
          break
        }
      }
      
      if (fl) {
        trts.tm <- c(trts.tm, trt)
      }
    }
  }
  
  trts <- trts.tm
  trts.all <- c(trts.all, trts.tm)
  
  if (length(trts) > 0) {
    indx <- which(Pf.trt.trt$data$Treatment %in% trts)
    rep.values <- c(rep.values, mean(Pf.trt.trt$data[[var.name]][indx]))
    vl <- sd(Pf.trt.trt$data[[var.name]][indx])/(length(indx)^0.5)
    rep.sd.values <- c(rep.sd.values, ifelse(is.na(vl), 0, vl))
    ns <- c(ns, length(indx))
    p.name <- c(p.name, pth)
  }  
}

par(mar=c(3, 19, 5, 2))
ord <- order(rep.values, decreasing = F)
D <- data.frame(pathway.name = p.name[ord], rep.cor = rep.values[ord], 
                err = rep.sd.values[ord], cnts = ns[ord])
D <- D %>% dplyr::mutate(pathway.name = paste(pathway.name, " (", cnts, ")", sep = ""))
g <- ggplot(D, aes(reorder(pathway.name, rep.cor), rep.cor))
g <- g + geom_bar(stat="identity", width = 0.4)
g <- g + geom_errorbar(aes(ymin = rep.cor - err, ymax = rep.cor + err), width = 0.35)
g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
               axis.text.y = element_text(size = 14))
g <- g + ylab("Average replicate correlation")
g <- g + xlab("")
g <- g + geom_hline(yintercept = hit.sel.thr, color = "red")
g <- g + coord_flip()
ggsave(sprintf("%s.pdf", p.type.name), g, width = 10, height = 12)


Pf.trt.trt$data$Pathway <- lapply(Pf.trt.trt$data$Pathway, function(x) str_replace_all(x, "Canonical ", "")) %>% unlist
Pf.trt.trt$data$Pathway[which(Pf.trt.trt$data$Pathway == "Transcription Factors")] <- "C-EBP alpha"

Pf.trt.trt$data %>% dplyr::filter(Treatment %in% trts.all) %>% 
  dplyr::group_by(Pathway) %>% 
  dplyr::summarise(strong.phenotype = length(which(Treatment %in% strongs.indx)), total.genes = n()) %>% 
  dplyr::arrange(-strong.phenotype/total.genes) %>%
  htmlTable::htmlTable()
