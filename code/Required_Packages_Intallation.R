## installing prereqs 

## recommended versions : 
## R version 3.2.2,  MASS 7.3.45, plyr 1.8.4, dplyr 0.5.0, ggplot2 2.2.0, caret 6.0.73, doMC 1.3.4, 
## yaml 2.1.14, xtable 1.8.2, reshape2 1.4.2, gridExtra 2.2.1, hash 2.2.6, devtools 1.12.0,
## knitr 1.15.1, stringr 1.1.0, magrittr 1.5, testthat 1.0.2, matrixStats 0.51.0, httr 1.2.1,
## XML 3.98.1.5, png 0.1.7, wordcloud 2.5, RCurl 1.95.4.8, corrplot 0.77, htmlTable 1.7, DBI 0.5.1,
## readr 1.0.0, forecast 7.3, ape 3.5, randomcoloR 1.0.0, progress 1.0.2, RMySQL 0.10.9, tiff 0.1.5,
## sampling 2.7, digest 0.6.10, robust 0.4.16, roxygen2 5.0.1, KEGGREST 1.10.1, EBImage 4.12.2,
## topGO 2.22.0, org.Hs.eg.db 3.2.3


## required versions : 
## GOSemSim 1.28.2, DOSE 2.8.3, clusterProfiler 2.4.3

req.packages <- c("MASS", "plyr", "dplyr", "ggplot2", "caret", 
                  "doMC", "yaml", "xtable", "reshape2", 
                  "gridExtra", "hash", "devtools", "knitr", "stringr",
                  "magrittr", "testthat", "matrixStats", "httr", 
                  "XML", "png", "wordcloud", "RCurl", "corrplot", 
                  "htmlTable", "DBI", "readr", "forecast", "ape", 
                  "randomcoloR", "progress", "RMySQL", "tiff", 
                  "sampling", "digest", "robust", "roxygen2")

for (pkg in req.packages) {
  if (!requireNamespace(pkg, quietly = T)) {
    install.packages(pkg, quiet = T)
  }
}

if (requireNamespace("GOSemSim", quietly = T)) {
  if (packageVersion("GOSemSim") != "1.28.2") {
    remove.packages("GOSemSim")
    install.packages("../aux_files/GOSemSim_1.28.2.tar", repos = NULL, type = "source", quiet = T)
  }
} else {
  install.packages("../aux_files/GOSemSim_1.28.2.tar", repos = NULL, type = "source", quiet = T)
}

if (requireNamespace("DOSE", quietly = T)) {
  if (packageVersion("DOSE") != "2.8.3") {
    remove.packages("DOSE")
    install.packages("../aux_files/DOSE_2.8.3.tar", repos = NULL, type = "source", quiet = T)
  } 
} else {
  install.packages("../aux_files/DOSE_2.8.3.tar", repos = NULL, type = "source", quiet = T)
}

if (requireNamespace("clusterProfiler", quietly = T)) {
  if (packageVersion("clusterProfiler") != "2.4.3") {
    remove.packages("clusterProfiler")
    install.packages("../aux_files/clusterProfiler_2.4.3.tar", repos = NULL, type = "source", quiet = T)
  }
} else {
  install.packages("../aux_files/clusterProfiler_2.4.3.tar", repos = NULL, type = "source", quiet = T)
}

bio.cond <- c("KEGGREST", "EBImage", "topGO", "org.Hs.eg.db", "DO.db")

for (pkg in bio.cond) {
  if (!requireNamespace(pkg, quietly = T)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkg)
  }
}

if (!requireNamespace("rwantshue", quietly = T)) {
  devtools::install_github("hoesler/rwantshue")  
}

cat("Testing the installed packages : \n")
flag <- T
for (pkg in req.packages) {
  if (!requireNamespace(pkg, quietly = T)) {
    message(sprintf("Error : %s did not install successfully.", pkg))
    flag <- F
  }
}

if (packageVersion("GOSemSim") != "1.28.2") {
  message("Error : GOSemSim version 1.28.2 is required.")
  flag <- F
}
if (packageVersion("DOSE") != "2.8.3") {
  message("Error : DOSE version 2.8.3 is required.")
  flag <- F
}
if (packageVersion("clusterProfiler") != "2.4.3") {
  message("Error : clusterProfiler version 2.4.3 is required.")
  flag <- F
}

if (!flag) {
  message("Installation failed.")
} else {
  message("Installed successfully.")
}