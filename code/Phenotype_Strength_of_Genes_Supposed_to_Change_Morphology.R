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

var.name <- "sim_pearson_q50"

Pf.trt.collapsed <- Pf.collapsed
matching.per.table <- readRDS("../results/master/ORFs_sequence_matching_transcripts_percentage/matching.per.table.rds")
leg.orfs <- matching.per.table %>% dplyr::filter(Nuc..Match.. > 99 & Prot..Match.. > 99) %>% dplyr::select(one_of("Treatment")) %>% as.matrix() %>% as.vector()

nulls.samp <- nulls.samp.rep.corr
thr.qn <- quantile(nulls.samp, 0.95)

Pf.trt.collapsed$data <- Pf.trt.collapsed$data %>% dplyr::filter(Treatment %in% leg.orfs)
Pf.trt.collapsed$data$Treatment <- Pf.trt.collapsed$data$Treatment %>% as.character()

i1 <- which(str_detect(Pf.trt.trt$data$Treatment, "RAC1"))
i2 <- which(str_detect(Pf.trt.trt$data$Treatment, "KRAS"))
i3 <- which(str_detect(Pf.trt.trt$data$Treatment, "CDC42"))
i4 <- which(str_detect(Pf.trt.trt$data$Treatment, "RHOA"))
i5 <- which(str_detect(Pf.trt.trt$data$Treatment, "PAK1"))
i6 <- which(str_detect(Pf.trt.trt$data$Pathway, "Hippo"))
i.total <- c(i1, i2, i3, i4, i5, i6)
i.total.complement <- setdiff(1:length(Pf.trt.trt$data$Treatment), i.total)

trts <- Pf.trt.trt$data$Treatment[i.total]
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

i.total <- which(Pf.trt.trt$data$Treatment %in% trts.tm)

trts <- Pf.trt.trt$data$Treatment[i.total.complement]
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

i.total.complement <- which(Pf.trt.trt$data$Treatment %in% trts.tm)

rep.of.genes <- c()
mag.of.genes <- c()
mag.all <- c()
rep.of.genes <- c(rep.of.genes, Pf.trt.trt$data[[var.name]][i.total])
Px.trt.collapsed <- Pf.trt.collapsed$data %>% dplyr::select(starts_with("PC"))
j.total <- c()

for (i in i.total) {
  j <- which(Pf.trt.collapsed$data$Treatment == Pf.trt.trt$data$Treatment[i])
  mag.of.genes <- c(mag.of.genes, (sum(Px.trt.collapsed[j, ] ^ 2))^0.5)
  j.total <- c(j.total, j)
}

for (i in i.total.complement) {
  j <- which(Pf.trt.collapsed$data$Treatment == Pf.trt.trt$data$Treatment[i])
  mag.all <- c(mag.all, (sum(Px.trt.collapsed[i, ] ^ 2))^0.5)
}

data <- rbind(c(which(rep.of.genes <= thr.qn) %>% length,  which(rep.of.genes > thr.qn) %>% length), 
              c(which(Pf.trt.trt$data[[var.name]][i.total.complement] <= thr.qn) %>% length,  which(Pf.trt.trt$data[[var.name]][i.total.complement] > thr.qn) %>% length))
fisher.test(data, alternative = "less", simulate.p.value = T) %>% print

tg2 <- sprintf("genes hypothesized to \n change morphology \n n=%d", length(rep.of.genes))
tg1 <- sprintf("other genes \n n=%d", length(Pf.trt.trt$data[[var.name]][i.total.complement]))
D.lab <- data.frame(tag = c(tg1, tg2), val = c(data[2, 2]/(sum(data[2,])), data[1, 2]/(sum(data[1,]))), 
                    ypos = 1.25)

D1 <- data.frame(val = rep.of.genes, tag = tg2)
D2 <- data.frame(val = Pf.trt.trt$data[[var.name]][i.total.complement], tag = tg1)
D <- rbind(D2, D1)
g <- ggplot(D, aes(tag, val, fill = tag))
g <- g + geom_violin(aes_string(colour = "tag", fill = "tag"), width = 0.2, trim = F)
g <- g + geom_hline(yintercept = thr.qn, 
                    linetype = "dashed",
                    color = "black")
g <- g + scale_y_continuous(limits = c(-0.8, 1.4), 
                            breaks = c(-.4,-.2,0,.2,.4,.6,.8, 1))
g <- g + scale_x_discrete(name = "", labels = c(tg1, tg2) )
g <- g + ylab("Median replicate correlation")
g <- g + geom_text(data = D.lab, 
                            aes(tag, 
                                ypos,
                                label = sprintf("%2d%%", round(val*100))),
                            size = 8, 
                            alpha = .7,
                            vjust = -0.95)

cnts <- D %>% dplyr::group_by(tag) %>% dplyr::summarise(cnt = n()) %>% dplyr::select(cnt) %>% as.matrix() %>% as.vector()

x.size <- 14
color.map <- c(rgb(0.4, 0.3, 0.7, 0.3), rgb(0.4, 0.3, 0.7, 1))
g <- g + scale_fill_manual(values = color.map)
g <- g + geom_hline(yintercept = 0, 
                      color = "gray", alpha = 0.8, size = .5)
g <- g + scale_color_manual(values = color.map)
g <- g + theme_bw() + theme(legend.position = "none",
               panel.background = element_blank(),
               panel.border = element_blank(),
               axis.line = element_line(size = .5, color = "gray") ,
               axis.title.y = element_text(size = 20, vjust = 1),
               axis.text.x = element_text(size = x.size),
               axis.text.y = element_text(size = 14),
               legend.text = element_text(size = 14),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_blank(),
               legend.title = element_blank(),
               panel.margin = unit(c(1,1,3,1), "cm")
)
g
ggsave("gene.morph.chng.pdf", g, width = 6, height = 6)
