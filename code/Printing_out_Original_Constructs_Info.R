## printing out original set of constructs

rm(list=ls())
plot.new()
frame()
library("dplyr")
library("devtools")
library("htmlTable")

CYTOMINR_DIR <- "../aux_files/cytominr_old/"
load_all(CYTOMINR_DIR)

Pf_org <- profile.data("../input/profiles/median.plus.mad.robstd.empty.well.profiles.yml") 
fac <- read.csv('../input/augmented_metadata.csv')
names(fac)[names(fac)=="ReferencePathway.Process"] <- "Pathway"
expect_equal(intersect(names(fac), names(Pf_org$data)), "Well")
new.fac.cols_org <- setdiff(names(fac), names(Pf_org$data))
Pf_org$data <- merge(Pf_org$data, fac, by=c("Well"))
Pf_org$factor_cols <- c(Pf_org$factor_cols, new.fac.cols_org)

constructs <- Pf_org$data %>% dplyr::filter(Type == "Treated") %>%
  dplyr::select(one_of(c("Treatment", "TARGETTRANS", "PUBLICID"))) %>% unique %>% dplyr::arrange(Treatment)
  
constructs %>% htmlTable(., rname = rep("", NROW(constructs)))
