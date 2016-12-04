## Set-up

- Clone the repository in say ~/tmp.
 
- Create a ~/.Rpackages file with this content:
```r
    list(
     default = function(x) {
       file.path("~/tmp",x)
     }
    )
```

- Start R, and install these packages
```r
    install.packages(c("doMC","yaml","xtable","reshape2","MASS","gridExtra","hash"))
    install.packages(c("devtools","testthat"))
    install.packages(c("knitr","markdown"))
    install.packages("caret", dependencies=TRUE)
```

- Go into the demo directory and run the demo
```r
    setwd("~/tmp/cytominr/demo")
    library(knitr)
    knit2html("hdac-demo.Rmd")
```