rm(list = ls())

library("knitr")
library("ggplot2")
library("dplyr")
library("devtools")

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

CYTOMINR_DIR <- "../aux_files/cytominr_old/"
load_all(CYTOMINR_DIR)

cor_method <- "pearson"
cor_tag <- "Pearson correlation"

f.violinplot <- function(D, 
                       D.lab, 
                       qthresh.null, 
                       legend.expanded, 
                       legend.shortened, 
                       ymin, 
                       ymax, 
                       color.map,
                       cor_tag, x.size = 14
                       ) {
  D.lab$qthresh.null <- qthresh.null
  D.lab$ypos <- 1.25 # ymax - .05
  p <- ggplot(D, aes(tag, val, fill = tag)) 
  p <- p + geom_hline(yintercept = 0, 
                      color = "gray", alpha = 0.8, size = .5)

  p <- p + geom_violin(aes_string(colour = "tag", fill = "tag"), width = 0.30, trim = F)

  tag.to.select <- D %>% dplyr::group_by(tag) %>% dplyr::summarise(cnt = n()) %>% dplyr::filter(cnt < 300) %>% dplyr::select(tag) %>% as.matrix() %>% as.vector() %>% as.character()
  tag.not.to.select <- setdiff(D$tag %>% unique, tag.to.select)
  
  D.subset <- D %>% dplyr::filter(tag %in% tag.to.select)
  if (length(tag.not.to.select) != 0) {
    D.subset.sampled <- D %>% dplyr::filter(tag == tag.not.to.select) %>% sample_n(., (log(NROW(D))/log(6000) * 300) %>% round)
    D.subset.sampled <- rbind(D.subset.sampled, D[c(which.min(D$val), which.max(D$val)), ])
  } else {
    D.subset.sampled <- D.subset
  }

  if (!is.na(qthresh.null)) {
    p <- p + geom_hline(yintercept = qthresh.null, 
                        linetype = "dashed",
                        color = "black")
    
  }

  p <- p + ylab(cor_tag)
  p <- p + xlab("")
  if (!any(is.na(D.lab$val))) {
    p <- p + geom_text(data = D.lab, 
                       aes(tag, 
                           ypos,
                           label = sprintf("%2d%%", round(val*100))),
                       size = 8, 
                       alpha = .7,
                       vjust = -.95)
    
  }

  p <- p + scale_fill_manual(values = color.map,
                             labels = legend.expanded)

  p <- p + scale_color_manual(values = color.map,
                             guide = FALSE)

  p <- p + scale_x_discrete(name = "", labels = legend.shortened )
  p <- p + theme_bw()
  p <- p + scale_y_continuous(limits = c(-0.8, 1.4), 
                              breaks = c(-.4,-.2,0,.2,.4,.6,.8, 1))
  p <- p + theme(legend.position = "none",
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
  p
}

# data frame to use for creating the plots
Pf <- Pf.collapsed 

# samples from the null distribution in rep. corr. analysis (taken from initial analysis)
nulls.samp <- nulls.samp.rep.corr

# setting up data needed for the plots
NPCS <- length(Pf$feat_cols)
top.pc.list <- paste("PC", seq(NPCS), sep = "")

Pf.trt.wctrl <- Pf.trt 
Pf.trt.wctrl.collapsed <- Pf.trt.strong.collapsed 

Pf.trt.wctrl$feat_cols <- top.pc.list
Pf.trt.wctrl$data <- 
  Pf.trt.wctrl$data[,c(Pf.trt.wctrl$factor_cols, Pf.trt.wctrl$feat_cols)] 

Pf.trt.wctrl.collapsed$feat_cols <- top.pc.list
Pf.trt.wctrl.collapsed$data <- 
  Pf.trt.wctrl.collapsed$data[,c(Pf.trt.wctrl.collapsed$factor_cols, Pf.trt.wctrl.collapsed$feat_cols)]
  
Pf.trt <- Pf.trt.wctrl 
Pf.trt$data <- 
  subset(Pf.trt.wctrl$data, Type == "Treated") 

Pf.trt.collapsed <- Pf.trt.wctrl.collapsed 
Pf.trt.collapsed$data <- 
  subset(Pf.trt.wctrl.collapsed$data, Type == "Treated") 

cor0 <- cor(t(feats(Pf.trt)), method = cor_method)
cor1 <- cor(t(feats(Pf.trt.collapsed)), method = cor_method)

m0.TreatmentAbbrev <- with(Pf.trt$data, outer(TreatmentAbbrev, TreatmentAbbrev, "=="))
m0.Gene <- with(Pf.trt$data, outer(Gene, Gene, "=="))

m1.Gene <- with(Pf.trt.collapsed$data, outer(Gene, Gene, "=="))
m1.TreatmentAbbrev <- with(Pf.trt.collapsed$data, outer(TreatmentAbbrev, TreatmentAbbrev, "=="))
isWt <- Pf.trt.collapsed$data$AlleleDesc %>% str_detect("WT") 
isMismatch <- Pf.trt.collapsed$data$AlleleDesc %>% str_detect("mismatch") 
isWtNotMismatch <- isWt & !isMismatch
m1.isWT <- outer(isWtNotMismatch, isWtNotMismatch, "&")

# Compute correlation between different-gene pairs
sel <- !m1.Gene 
sel <- sel & upper.tri(sel)
cor.collapsed.diff.gene <- cor1[sel]

# Compute correlation between treatments that are not replicates, not same-gene
sel <- !m0.TreatmentAbbrev & !m0.Gene 
sel <- sel & upper.tri(sel)
cor.diff.pert.gene <- nulls.samp
non.edge <- lapply(Pf.untrt.cor$data$Well, function(x) !((str_sub(x, 1, 1) %in% c("a", "p")) | (str_sub(x, 2, 3) %in% c("01", "24")) ) ) %>% unlist
cor.same.pos.untrt <- Pf.untrt.cor$data$sim_pearson_q50[non.edge]

# Compute correlation between treatments that are treated with the same
# perturbation (replicates)
sel <- m0.TreatmentAbbrev 
sel <- sel & upper.tri(sel)
cor.same.pert <- Pf.trt.trt$data$sim_pearson_q50

# Compute correlation between same-gene, but different-construct pairs, only WT
sel <- m1.Gene & !m1.TreatmentAbbrev & m1.isWT
sel <- sel & upper.tri(sel)
cor.collapsed.same.gene.diff.construct <- cor1[sel]

# 1) Fig 1A : rep. corr. analysis (null = non-replicate corr.),   
# 2) Fig 1B : WT clones vs. non-clones, 
# 3) Supp. Fig. 1 : rep. corr. analysis (null = untreated rep. corr.),

NULLTHRESH <- .95
NULLTHRESH.3 <- .95

qthresh.null.1 <- quantile(cor.diff.pert.gene, NULLTHRESH)[[1]]
xthresh.1 <- 1 - ecdf(cor.same.pert)(qthresh.null.1)

qthresh.null.3 <- quantile(cor.same.pos.untrt, NULLTHRESH.3)[[1]]
xthresh.3 <- 1 - ecdf(cor.same.pert)(qthresh.null.3)

# D1
D1 <- rbind(
  data.frame(tag = "cor.diff.pert.gene", val = cor.diff.pert.gene), 
  data.frame(tag = "cor.same.pert", val = cor.same.pert)
)

D1.lab <- rbind(
  data.frame(tag = "cor.diff.pert.gene", 
             val = 1 - NULLTHRESH),
  data.frame(tag = "cor.same.pert", 
             val = xthresh.1))
             
# D3
D3 <- rbind(
  data.frame(tag = "cor.same.untrt", val = cor.same.pos.untrt), 
  data.frame(tag = "cor.same.pert", val = cor.same.pert)
)

D3.lab <- rbind(
  data.frame(tag = "cor.same.untrt", 
             val = 1 - NULLTHRESH.3),
  data.frame(tag = "cor.same.pert", 
             val = xthresh.3))

# D2
D2 <- rbind(
  data.frame(tag = "cor.collapsed.diff.gene", 
             val = cor.collapsed.diff.gene),
  data.frame(tag = "cor.collapsed.same.gene.diff.construct", 
             val = cor.collapsed.same.gene.diff.construct))

qthresh.null.2 <- quantile(cor.collapsed.diff.gene, .95)[[1]]
xthresh.gene.2 <- 1 - ecdf(cor.collapsed.same.gene.diff.construct)(qthresh.null.2)

D2.lab <- rbind(
  data.frame(tag = "cor.collapsed.diff.gene", 
             val = 1 - NULLTHRESH),
  data.frame(tag = "cor.collapsed.same.gene.diff.construct", 
             val = xthresh.gene.2))


get_title <- function(D, lut) {
  t1 <- unique(D$tag)[1]
  t2 <- unique(D$tag)[2]
  n1 <- sum(D$tag == t1)
  n2 <- sum(D$tag == t2)
  return("")
}

get_n1_n2 <- function(D) {
  t1 <- unique(D$tag)[1]
  t2 <- unique(D$tag)[2]
  n1 <- sum(D$tag == t1)
  n2 <- sum(D$tag == t2)
  return(c(n1, n2))
}

v1 <- get_n1_n2(D1)
v2 <- get_n1_n2(D2)
v3 <- get_n1_n2(D3)

legend.expanded.1 <- 
  c("cor.diff.pert.gene" = 
      "non-replicate pairs",
    "cor.same.pert" = 
      "replicate pairs (summarized)")

legend.shortened.1 <- c("cor.diff.pert.gene" = sprintf("pairs of \n different constructs \n n=%d", v1[1]), 
                        "cor.same.pert" = sprintf("replicates of \n the same construct \n n=%d", v1[2]))

legend.expanded.2 <- 
  c("cor.collapsed.diff.gene" = 
      "different genes",
    "cor.collapsed.same.gene.diff.construct" = 
      "same gene targeted with different constructs")

legend.shortened.2 <- c("cor.collapsed.diff.gene" = sprintf("pairs of \n different genes \n n=%d", v2[1]),
                        "cor.collapsed.same.gene.diff.construct" = sprintf("pairs of constructs \n containing the same gene \n n=%d", v2[2]))

legend.expanded.3 <- 
  c("cor.same.untrt" = 
      "untreated pairs on \n the same position (summarized)",
    "cor.same.pert" = 
      "replicate pairs (summarized)")

legend.shortened.3 <- c("cor.same.untrt" = sprintf("pairs of the untreated \n in the same position \n n=%d", v3[1]),
                        "cor.same.pert" = sprintf("pairs of the \n same construct \n n=%d", v3[2]))

p1 <- f.violinplot(D1, D1.lab, qthresh.null.1, 
                legend.expanded.1, 
                legend.shortened.1, 
                ymin = min(-0.45, D1$val), 
                ymax = max(0.8, D1$val),
                color.map = c(rgb(0.3, 0.7, 0.4, 0.3), rgb(0.3, 0.7, 0.4, 1)),
                cor_tag)
p1 <- p1 + ggtitle(get_title(D1, legend.shortened.1))

p3 <- f.violinplot(D3, D3.lab, qthresh.null.3, 
                legend.expanded.3, 
                legend.shortened.3, 
                ymin = min(-0.45, D3$val), 
                ymax = max(0.8, D3$val),
                color.map = c(rgb(0.3, 0.7, 0.4, 0.3), rgb(0.3, 0.7, 0.4, 1)),
                cor_tag, x.size = 10)
p3 <- p3 + ggtitle(get_title(D3, legend.shortened.3))

p2 <- f.violinplot(D2, D2.lab, qthresh.null.2, 
                legend.expanded.2, 
                legend.shortened.2, 
                ymin = min(-0.45, D2$val), 
                ymax = max(0.8, D2$val),
                color.map = c(rgb(0.7, 0.3, 0.4, 0.3), rgb(0.7, 0.3, 0.4, 1)), cor_tag = "")
p2 <- p2 + ggtitle(get_title(D2, legend.shortened.2))
ggsave(sprintf("Fig1A.pdf"), p1, 
       width = 6, height = 6)
ggsave(sprintf("FigSupp1.pdf"), p3, 
       width = 6, height = 6)
ggsave(sprintf("Fig1B.pdf"), p2, 
       width = 6, height = 6)