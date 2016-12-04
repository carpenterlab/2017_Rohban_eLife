rm(list = ls())

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(htmlTable)
library(gridExtra)
library(ggplot2)
library(grid)

## change 'sm' variable to the content of one of the following files 
## sought upregulation of three YAP/TAZ targets : result_UPTAG_summly_n7503_upregulate_CYR61_CTGF_BIRC5__YAP_TAZ_Trad_Targets.csv
## sought upregulation and downregulation of predicted YAP/TAZ targets : result_UPTAG_summly_n7503_upregulation_of_INPP4B_MAP7_LAMA3_STMN1_TRAM2_downregulation_SPP1_IER3_RAB31_GPR56__YAP_TAZ_predicted_Targets.csv
## sought upregulation and downregulation of predicted TRAF2/REL targets : result_UPTAG_summly_n7503_upregulation_of_NFKBIA_IKBKE_AKAP8_BIRC2_downregulation_of_RPA3__TRAF2_REL_predicted_Targets.csv
sm <- readr::read_csv("../results/master/LINCS_Matches_L1000/result_UPTAG_summly_n7503_upregulation_of_CYR61_CTGF_BIRC5__YAP_TAZ_Trad_Targets.csv")

## a description of the query 
desc <- "YAP/TAZ targets"

## pathway id in KEGG; other alternatives : TNF alpha = hsa04668, NFkB = hsa04064, Cell Cycle = hsa04110, JAK/STAT = hsa04630, MAPK = hsa04010, Hippo = hsa04390
pathway.to.check.for.enrichment <- "hsa04064"

## the mapping used for going from 
gene.mapping <- as.list(org.Hs.egALIAS2EG)

## filtered only for overexpression constructs
sm.filtered <- sm %>% dplyr::filter(pert_type == "trt_oe") %>% arrange(mean_rankpt_4) %>% dplyr::select(one_of(c("pert_iname", "mean_rankpt_4"))) 
gene.list <- sm.filtered[,"pert_iname"] %>% as.matrix %>% as.vector()
gene.id.list <- lapply(gene.mapping[gene.list], function(x) (x[1])) %>% unlist %>% as.matrix() %>% as.vector()

## defining the phenotype in the GSEA plot
phenotype <- -sm.filtered[,"mean_rankpt_4"] %>% as.matrix() %>% as.vector()
dt <- phenotype
names(dt) <- (gene.id.list %>% as.character())

## do the KEGG enrichment analysis
kk2 <- gseKEGG(geneList     = dt,
               organism     = 'hsa',
               nPerm        = 20000,
               pvalueCutoff = 1.0,
               verbose      = FALSE)
print(summary(kk2))
u <- gseaplot(kk2, geneSetID = pathway.to.check.for.enrichment)
u$runningScore$labels$x <- "Position in the Ranked List of ORFs"
u$position$labels$y <- sprintf("matching score to %s", desc)
u$position$data$geneList <- -u$position$data$geneList
g1 <- u$position
g2 <- u$runningScore
g1 <- ggplot_gtable(ggplot_build(g1))
g2 <- ggplot_gtable(ggplot_build(g2)) 
maxWidth <- unit.pmax(g1$widths[2:3], g2$widths[2:3])
g1$widths[2:3] <- maxWidth
g2$widths[2:3] <- maxWidth
quartz(width = 15, height = 10)
grid.arrange(g2, g1, heights = c(1.33, 1))
dev.print(device = pdf, "GSEA_ORF.pdf")