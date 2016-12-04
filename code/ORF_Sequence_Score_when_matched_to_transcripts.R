## note : this code requires to be run on a platform connected to the Broad Institute internal network

rm(list = ls())
library(dplyr)
library(devtools)

CYTOMINR_DIR <- "../aux_files/cytominr_old/"
load_all(CYTOMINR_DIR)

prune.plates.list <- data.frame(Plate=c('P2'))
Pf <- profile.data("../input/profiles/median.plus.mad.robstd.empty.well.profiles.yml") 
fac <- read.csv('../input/augmented_metadata.csv')
names(fac)[names(fac)=="ReferencePathway.Process"] <- "Pathway"
expect_equal(intersect(names(fac), names(Pf$data)), "Well")
new.fac.cols_org <- setdiff(names(fac), names(Pf$data))
Pf$data <- merge(Pf$data, fac, by=c("Well"))
Pf$factor_cols <- c(Pf$factor_cols, new.fac.cols_org)

prune.wells.list <- data.frame(Plate=c('P3', 'P6', 'P2'), 
                               Well=c('b01', 'e17', 'j13'))
Pf <- prune.outliers(Pf, read.outlier.data(Pf))
Pf <- prune.rows(Pf, !(with(Pf$data, paste(Plate, Well, sep="_")) %in% 
                         with(prune.wells.list, 
                              paste(Plate, Well, sep="_"))))
Pf <- prune.rows(Pf, !(Pf$data$Plate %in% 
                                 prune.plates.list$Plate))
Pf <- create.plate.pos(Pf, modify.data=TRUE)

ORFs <- Pf$data %>% dplyr::filter(Type == "Treated") %>% dplyr::select(one_of("Treatment")) %>% as.matrix() %>% as.vector() %>% unique

## removed following ORFs for technical/plate position/QC reasons  within plate issues
ORFs <- ORFs[which(!str_detect(ORFs, "mismatch") & ! ORFs %in% c("CREB1_WT", "MAPK1_WT.1", "CYLD_WT"))]
match.list <- c() 
pb <- progress::progress_bar$new(total = length(ORFs))

for (trt in ORFs) {
  brd.id <- Pf$data %>% dplyr::filter(Treatment == trt) %>% dplyr::select(PUBLICID) %>% unique %>% as.matrix() %>% as.vector()
  v <- read.csv(sprintf("https://iwww.broadinstitute.org/rnai/db/clone/details?view=csv&grid=1&cloneId=%s", brd.id))
  if (!("Taxon" %in% colnames(v)) || is.null(v) || NROW(v) < 1 || NCOL(v) < 2) {
    pb$tick()
    next
  }
  w <- v %>% dplyr::filter(Taxon == "human") %>% dplyr::mutate(Treatment = trt)
  match.list <- c(match.list, list(w[1,]))
  pb$tick()
}

matching.per.table <- do.call(rbind, match.list)
saveRDS(matching.per.table, "../results/master/ORFs_sequence_matching_transcripts_percentage/matching.per.table.rds")
write.csv(matching.per.table, "../results/master/ORFs_sequence_matching_transcripts_percentage/matching.per.table.csv")