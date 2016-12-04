
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  # http://www.r-bloggers.com/memory-management-in-r-a-few-tips-and-tricks/
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

confusionMatrix.loocv <- function(fit) {
  params <- gsub("^\\.", "", names(fit$bestTune))
  selvec <- !vector(length=NROW(fit$pred))
  if (params=="parameter") {
   selvec <- !vector(length=NROW(fit$pred)) 
  } else {
    for (p in params) {
      selvec <- selvec & (fit$pred[p] == fit$bestTune[[paste(".",p,sep="")]])
    }
  }
  
  cfmat <- (melt(table(fit$pred[selvec,1:2])))
  names(cfmat)[names(cfmat)=="pred"] <- "Prediction"
  names(cfmat)[names(cfmat)=="obs"] <- "Reference"
  cfmat
}

plot.confmat <- function (cfmat, print_value=FALSE, rotate_labels=FALSE) {
  expect_true(all(levels(cfmat$Reference) == levels(cfmat$Prediction)))
  cfmat <- ddply(cfmat, .(Reference), transform, 
    Percent = 100 * value / sum(value))
#cfmat <- ddply(cfmat, .(Reference), transform, Percent = value)
  p <- ggplot() + geom_tile(aes(x=Reference, y=Prediction,fill=Percent),
                            data=cfmat, color="black",size=0.1) 
  p <- p + labs(x="Actual",y="Predicted") 
  p <- p + coord_equal()
  p <- p +  geom_text(aes(x=Reference, y=Prediction,
                          label=sprintf("%.0f", Percent)),
  data=subset(cfmat, Percent >= 1), 
                      size=4, colour="black") +
                        scale_fill_gradient(
                          low="white",high="gray",limits=c(0,100))
  
  if (print_value) {
    p <- p +  geom_text(aes(x=Reference, y=Prediction,
                            label=sprintf("%d", value)),
    data=subset(cfmat, Percent >= 1),
                        size=3, colour="blue",vjust=3)   
  }
  
  p <- p + theme_bw()
  p <- p + theme(legend.position = "none",
                 plot.margin = unit(c(0, 0, 2, 0), "cm"),
                 axis.title.x = element_text(vjust = -3)
                 #axis.text.y  = element_text(hjust = 0)
                 )
  
  if (rotate_labels) {
    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  }
  return(p)
}




sintegral = function(x, fx, n.pts = max(256, length(n.x))){
    ## numerically integrates fx over x using Simpsons rule
    ## x - a sequence of x values
    ## fx - the value of the function to be integrated at x
    ##    - or a function
    ## n.pts - the number of points to be used in the integration
    ## ret - if true returns the partial sums of the integration
    ## note this argument is deprecated and ignored

    #message("Note: sintegral's behavior has changed.\nTo get the value of the 
    # integral use sintegral(x,fx)$value.\nTo get the cdf use sintegral(x,fx)$cdf")

    if(class(fx) == "function")
        fx = fx(x)

    n.x = length(x)

    if(n.x != length(fx))
        stop("Unequal input vector lengths")

    if(n.pts < 64)
        n.pts = 64

    ## use linear approximation to get equally spaced x values


    ap = approx(x, fx, n = 2 * n.pts + 1)

    h = diff(ap$x)[1]

    integral = h*(ap$y[2 * (1:n.pts) - 1] +
                  4 * ap$y[2 * (1:n.pts)] +
                  ap$y[2 * (1:n.pts) + 1]) / 3


    invisible(list(value = sum(integral),
                   cdf = list(x = ap$x[2*(1:n.pts)], y = cumsum(integral))))
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

getd <- function(d, k, def) if (k %in% names(d)) d[[k]] else def

convert_str_to_dataframe <- function(s) {
  dfnames <- strsplit(s[[1]], split=",")[[1]]
  dfs <- strsplit(s[[2]], split=";")[[1]]
  df <- ldply(strsplit(dfs, split=';'), 
    function(row) data.frame(t(strsplit(trim(row), split=',')[[1]])))
  names(df) <- dfnames
  df
}

library(robust)
corRobv <- function(x, y) covRob(cbind(x,y), corr=T)$cov[1,2]

dcorm <- function(X) {
  n <- NCOL(X)
  C <- matrix(nrow=n, ncol=n, dimnames=list(names(X), names(X)))
  for (i in seq(n-1)) 
    for (j in seq(i+1,n)) 
      C[i,j] <- C[j,i] <- dcor(X[,i], X[,j])
  
  for (i in seq(n))
    C[i,i] <- 1
  C
}

dauc <- function(sig, nul=NULL)  {
  # Areas under signal curve minus area under null
  if (is.null(nul)) {
    nul.val <- 0.5
    message('Assuming AUC for nul = 0.5')
  } else {
    nul.val <- sintegral(c(0,seq_along(nul)/length(nul)), c(0,nul))$value
  }
  sintegral(c(0,seq_along(sig)/length(sig)), c(0,sig))$value - nul.val
}

plot.liftcurve <- function(sig, 
                           nul, 
                           size = 1, 
                           signal.name = "Signal",
                           null.name = "Null",
                           n.breaks = 5) {
  sig <- unname(sig)
  nul <- unname(nul)
  rs <- c(0,seq_along(sig)/length(sig))*100
  rn <- c(0,seq_along(nul)/length(nul))*100
  s <- c(0,sig) * 100
  n <- c(0,nul) * 100
  dfs <- data.frame(rank=rs, value=s, name=signal.name)
  dfn <- data.frame(rank=rn, value=n, name=null.name)
  df <- rbind(dfs, dfn)
  # Reorder factor so that signal is first
  df$name <- factor(df$name, levels = c(as.character(signal.name), 
                                        as.character(null.name)))


  p <- ggplot(df, aes(rank, value, color=name)) +  geom_line(size=size)
  p <- p + annotate("text", x = 80, y = 5, 
                    label = sprintf('dAUC=%0.2f', dauc(sig, nul)))
  p <- p + xlab('Percentile') + ylab('Percent Positive')
  p <- p + scale_x_continuous(breaks=pretty_breaks(n=n.breaks)) 
  p <- p + scale_y_continuous(breaks=pretty_breaks(n=n.breaks)) 
  p <- p + guides(color=guide_legend(title=NULL, nrow=2)) 
  p + theme_bw()
  
}


robust_linear_norm <- function(v) {
  q <- quantile(v, c(.01,.99))
  u <- (v - q[[1]]) / (q[[2]] - q[[1]])
  return (u)
}

standardization_norm <- function(v) {
  return(scale(v))
}

dummy_norm <- function(v) {
  return(v)
}


corwt1 <- function(X, neg_const = 0.01) {
  X <- as.matrix(X)
  Cm <- cor(t(X))
  Cm[Cm < 0] <- neg_const
  Cm <- colSums(Cm)-1
  Cm <- Cm / sum(Cm)
  colSums(X * Cm)
}


cos.dist <- function(X) {
  n <- nrow(X) 
  cmb <- expand.grid(i=1:n, j=1:n) 
  cos.dist.helper <- function(ix) {
    A = X[ix[1],]
    B = X[ix[2],]
    return (1-( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) ))
  }  
  C <- matrix(aaply(as.matrix(cmb),1,cos.dist.helper,
    .parallel=F, .progress = "text"),n,n)
}



plotmatrix3 <- function (data, lab)
{
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], lab = lab)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  mapping <- aes()
  mapping <- defaults(mapping, aes_string(x = "x", y = "y",color="lab"))
  class(mapping) <- "uneval"
  ggplot(all, mapping) + facet_grid(xvar ~ yvar, scales = "free") + 
  geom_point( na.rm = TRUE)
}



plot.pdf.sig.nul.helper <- function(P.sig, 
                             P.nul,
                             sig.name = NULL,
                             nul.name = NULL,
                             meas.name = "sim_spearman_q50") {
  if (is.null(sig.name) & ("label" %in% names(P.sig)))    
    sig.name = P.sig$label
  else
    sig.name = "Signal"
  
  if (is.null(nul.name) & ("label" %in% names(P.nul)))    
    nul.name = P.nul$label
  else
    nul.name = "Null"
  
  P.sig$data$grouping <- sig.name
  P.nul$data$grouping <- nul.name
  D <- rbind(P.sig$data[,c(meas.name, "grouping")],
             P.nul$data[,c(meas.name, "grouping")])
  D$meas.name <- D[,meas.name]
  
  # Order the levels  
  D$grouping <- factor(D$grouping, levels = c(as.character(sig.name), 
                                              as.character(nul.name)))
  
  
  
  # Get thresholds wrt Null
  yte <- 0.95
  sig <- subset(D, grouping==sig.name)$meas.name
  nul <- subset(D, grouping==nul.name)$meas.name
  xt <- quantile(nul, yte)[[1]]
  yt <- ecdf(sig)(xt)
  
  list(D=D, P.sig=P.sig, P.nul=P.nul, yte=yte, xt=xt, yt=yt, sig=sig, nul=nul,
       sig.name=sig.name, nul.name=nul.name)
  
}

plot.bar.sig.nul <- function(P.sig, 
                             P.nul,
                             sig.name = NULL,
                             nul.name = NULL,
                             meas.name = "sim_spearman_q50",
                             plot.dotplot = FALSE) {
  
  r <- plot.pdf.sig.nul.helper(P.sig, 
                               P.nul,
                               sig.name,
                               nul.name,
                               meas.name)
  p <- ggplot(r$D, aes(x = grouping, y = meas.name, fill = grouping)) + 
    geom_boxplot(outlier.size=0) 
  if (plot.dotplot) {
    p <- p + geom_dotplot(binaxis = "y",
                          stackdir = "center",
                          binwidth = .02,
                          dotsize = .5,
                          stackratio = 1,
                          alpha=0.5)      
  } else {
    p <- p + geom_point(position = position_jitter(width = .1), alpha=.5)
  }
  p <- p + geom_hline(yintercept=r$xt, linetype="dotted")
  p <- p + annotate("text", 
                    y = r$xt, 
                    x = .5,
                    label = sprintf("x = %.2f", r$xt), 
                    size = 4, 
                    alpha = 0.5)
  p <- p + theme_bw()
  p <- p + ylab(meas.name) + xlab("") 
  p <- p + coord_flip()
  p <- p + guides(fill=guide_legend(title=NULL, nrow=2))
  p <- p + theme(axis.text.y=element_blank())
  p <- p + theme(legend.position = "bottom")
  p
  
}
  
  
plot.pdf.sig.nul <- function(P.sig, 
                             P.nul,
                             sig.name = NULL,
                             nul.name = NULL,
                             meas.name = "sim_spearman_q50") {

  # Plot density functions of signal and null   
  r <- plot.pdf.sig.nul.helper(P.sig, 
                               P.nul,
                               sig.name,
                               nul.name,
                               meas.name)
  color.map <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, 
                               h.start = 0, direction = 1)
  

  # Get maximum height of the density curves
  ymax <- max(ddply(r$D, .(grouping), 
                    function(X) {
                      if (NROW(X) > 1) {
                        max(density(X$meas.name)$y)
                      } else {
                        0
                      }
                    })$V1)
  
  p <- ggplot()
  p <- p + stat_density(data=r$D, 
                        aes(meas.name, color=grouping),
                        position="identity", geom="line", size=2)
  p <- p + geom_rug(data=subset(r$D, grouping==r$sig.name), aes(meas.name), 
                    color = color.map(2)[[1]], sides = "b" )
  p <- p + geom_rug(data=subset(r$D, grouping==r$nul.name), aes(meas.name), 
                    color = color.map(2)[[2]], sides = "t" )
  p <- p + geom_vline(xintercept=r$xt, linetype="dotted")
  p <- p + annotate("text", x=r$xt, y=.8*ymax,label=
                sprintf("x = %.2f \nNull = %2d%%\nTreatment = %2d%% (%d/%d)", 
                              r$xt, 
                              100-ceiling(r$yte*100), 
                              100-ceiling(r$yt*100),
                              sum(r$sig > r$xt),
                              length(r$sig)), 
                    size=4, alpha=0.5)

  p <- p + ylab("density") + xlab("spearman")
  p <- p + theme_bw()
  p <- p + guides(color=guide_legend(title=NULL, ncol=2))
  p <- p + theme(legend.position = "bottom")
  p

}


plot.liftcurve.sig.nul <- function(P.sig, 
                                   P.nul,
                                   sig.name = NULL,
                                   nul.name = NULL,
                                   size = 1) {
  # Get the default ggplot2 colormap
  color.map <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, 
                             h.start = 0, direction = 1)

  if (is.null(sig.name) & ("label" %in% names(P.sig)))    
    sig.name = P.sig$label
  else
    sig.name = "Signal"

  if (is.null(nul.name) & ("label" %in% names(P.nul)))    
    nul.name = P.nul$label
  else
    nul.name = "Null"

  p <- plot.liftcurve(P.sig$lc, 
                      P.nul$lc, 
                      signal.name=sig.name, 
                      null.name=nul.name, 
                      size=size) 
  p <- p + theme(legend.position = "bottom")
  p 
  
}

library(scales)
plot.ecdf.rank.sig.nul <- function(P.sig, 
                                   P.nul,
                                   sig.name = NULL,
                                   nul.name = NULL,
                                   size = 1,
                                   n.breaks = 5) {
  if (is.null(sig.name) & ("label" %in% names(P.sig)))    
    sig.name = P.sig$label
  else
    sig.name = "Signal"

  if (is.null(nul.name) & ("label" %in% names(P.nul)))    
    nul.name = P.nul$label
  else
    nul.name = "Null"

  f <- function(P, tag) {
    v <- apply(P$facbin2, 1, function(x) {
      if(sum(x)==0)
        NULL
      else
        which(x==1)[1]/(length(x)+1)
    })
    expect_false(any(is.null(v)))
    data.frame(grouping = rep(tag, length(v)), value = v)
  }

  D <- rbind(f(P.nul, nul.name), f(P.sig, sig.name))
  # Order the levels  
  D$grouping <- factor(D$grouping, levels = c(as.character(sig.name), 
                                              as.character(nul.name)))
  

  p <- ggplot(D, aes(value, color=grouping)) + stat_ecdf(size=size)
  p <- p + ylab("Percentage that find a match") 
  p <- p + xlab("Percentile for finding first match")
  p <- p + scale_x_continuous(labels = percent, limits=c(-.1, 1.1), 
                              breaks=pretty_breaks(n=n.breaks))
  p <- p + scale_y_continuous(labels = percent, limits=c(-.1, 1.1), 
                              breaks=pretty_breaks(n=n.breaks))
  p <- p + theme_bw()
  p <- p + guides(color=guide_legend(title=NULL, nrow=2))
  p <- p + theme(legend.position = "bottom")
  p

}


summarize.sig.nul <- function(P.sig, 
                              P.nul,
                              cor.feat.name = 'sim_spearman_q50',
                              q.null = 0.95) {
  # Summarize differences between signal and nul
  # Assumes that lift-curves and correlations have been computed
  
  # Summarize lift-curve statistics
  dauc.val <- dauc(P.sig$lc$lc, P.nul$lc$lc)
  matdim.sig <- P.sig$lc$matdim
  matdim.nul <- P.nul$lc$matdim
  
  # Summarize correlation statistics
  cor.sig <- P.sig$cor$data[,cor.feat.name]
  cor.nul <- P.nul$cor$data[,cor.feat.name]
  
  cor.sig.med <- median(cor.sig)
  cor.nul.med <- median(cor.nul)
  t.nul <- quantile(cor.nul, q.null)[[1]]
  pc.above.nul <- 1-ecdf(cor.sig)(t.nul)
  
  n.set.sig.tot <- length(cor.sig)
  n.set.sig.pos <- sum(cor.sig > t.nul)
  
  data.frame(dauc.val, 
             n.sig = as.integer(matdim.sig[2]+1), 
             n.nul = as.integer(matdim.nul[2]+1),
             cor.sig.med,
             cor.nul.med,
             pc.above.nul,
             n.set.sig.pos = as.integer(n.set.sig.pos),
             n.set.sig.tot = as.integer(n.set.sig.tot))
}


generate.null <- function(D, permute.cols, constrain.groups = F)  {
  # Create null distribution by permuting columns

  permute.cols <- as.character(permute.cols)
  
  # If there are multiple columns, they will be permuted together
  # Merge into one column
  D$data$permute.col <- as.vector(apply(D$data[permute.cols], 1, paste,
                                        collapse="---"))

  D.d <- D$data
  
  if (!constrain.groups) {
    
    new.group.ids <- sample(D.d[, 'permute.col'])
    
  } else {
    
    # Create null distribution by permuting a group column
    # making sure that no pair of elements that are in the same group after
    # permutation were originally in the same group before permutation
    # Requires that all groups are of the same size
    tbl <- table(D.d[,"permute.col"])
    n <- unique(tbl) # number of elements in each group 
    expect_true(length(n)==1)  # all groups should be of the same size  
    k <- length(tbl) # number of groups
    
    # shuffle the rows (no permutation yet)
    D.d <- D.d[sample(NROW(D.d)),]
    
    # order the levels for group variable randomly
    D.d <- D.d[order(D.d[,'permute.col']), ]
    D.d$order.col <- rep(sample(k), each=n)
    D.d <- D.d[order(D.d$order.col), ]
    D.d$order.col <- NULL
    
    # permute the groups
    idx <- ( (rep(seq(n)-1, k) + rep(seq(k)-1, each=n)) %% k ) + 1
    new.group.ids <- unique(D.d[,'permute.col'])[idx]
  }
  
  D.d[,'permute.col'] <- new.group.ids
  # Assign values to the original set of columns
  for (i in seq_along(permute.cols)) {
    ic <- as.character(permute.cols[i])
    vc <- as.vector(sapply(D.d[,'permute.col'], 
                           function(s) str_split(s, '---')[[1]][i]))
    D.d[,ic] <- vc
  }
  D.d[,'permute.col'] <- NULL
  D$data <- D.d  
  D
}

handle.null.na.empty <- function(x) {
  if (length(x)==0)
    x <- NULL
  if (!is.null(x)) {
    if (any(is.na(x))) {
      if (length(x)==1)
        x <- NULL
    }
  }
  if (!is.null(x)) 
    x <- as.character(x)  
  x
}

compute.lift.cor.helper <- function(Px,
                                    grouping.cols, 
                                    exclude.cols = NULL,
                                    select.cols = NULL,
                                    distfunc = "spearman",
                                    summarize.quantiles = TRUE) {
  
  grouping.cols <- handle.null.na.empty(grouping.cols)
  expect_false(is.null(grouping.cols))
  exclude.cols <- handle.null.na.empty(exclude.cols)
  select.cols <- handle.null.na.empty(select.cols)

  # Compute lift curve and correlation  
  
  Pf.lc <- liftcurve(Px,
                     fac_col = grouping.cols,
                     fac_col_exclude = exclude.cols,
                     select.cols = select.cols,
                     distfunc = distfunc)
  
  Pf.cor <- qualitymeas(Px,
                        metacols = grouping.cols,
                        exclude.cols = exclude.cols,
                        select.cols = select.cols,
                        cmpfunc = distfunc,
                        summarize.quantiles = summarize.quantiles)
  
  list(lc = Pf.lc, cor = Pf.cor)
  
}

tableGrob.format.columns <- function(tbl) {
  format.tbl <- rbind()
  for(cn in names(tbl)) {
    tcn <- tbl[[cn]]
    expect_true(is.double(tcn) | is.integer(tcn))
    if(is.double(tcn)) {      
      fmt.str <- "%.2f"
    } else if(is.integer(tcn)) {
      fmt.str <- "%d"
    }
    
    format.tbl.i <- matrix(sprintf(fmt.str, tcn))
    row.names(format.tbl.i) <- c(cn)
    format.tbl <- rbind(format.tbl, format.tbl.i)
  }
  tableGrob(format.tbl, core.just="right")
  
}

combine.qualitymeas <- function(P1, P2) {
  if (is.null(P1))
    return(P2)
  if (is.null(P2))
    return(P1)
  expect_true(all(P1$feat_cols==P2$feat_cols))
  expect_true(all(P1$factor_cols==P2$factor_cols))
  expect_true(NCOL(P1$data)==NCOL(P2$data))
  P1$data <- rbind(P1$data, P2$data)
  P1
}

combine.liftcurve <- function(P1, P2) {
  if (is.null(P1))
    return(P2)
  if (is.null(P2))
    return(P1)

  expect_true(all(NCOL(P1$facbin1)==NCOL(P2$facbin1)))
  expect_true(all(NCOL(P1$facbin2)==NCOL(P2$facbin2)))
  expect_true(all(NCOL(P1$facbin3)==NCOL(P2$facbin3)))
  expect_true(all(names(P1$facp)==names(P2$facp)))
  
  P1$facbin1 <- rbind(P1$facbin1, P2$facbin1)
  P1$facbin2 <- rbind(P1$facbin2, P2$facbin2)
  P1$facbin3 <- rbind(P1$facbin3, P2$facbin3)
  P1$lc <- colMeans(P1$facbin3)
  P1$matdim <- dim(P1$facbin3)
  P1
}

apply.pca <- function (cfg.in.fname, pca.model, nzv.cols) {
  
  sha.model <- str_sub(digest(list(pca.model, nzv.cols)), 1,8)
  
  # Load profile.data 
  Pf <- profile.data(cfg.in.fname) 
  
  # Handle Plate metadata
  Pf$data$Plate <- Pf$data$Plate_
  Pf$data$Plate_ <- NULL
  Pf$factor_cols <-  Pf$factor_cols[Pf$factor_cols != "Plate_"]
  
  # remove features previously identified as having near-zero variance
  Pf <- prune.feats.nzv.apply(Pf, nzv.cols)
  
  # apply PCA based on previous learned model
  Pf <- pca.apply(Pf, pca.model)
  
  # Save feats and factors to variables
  Pf.factors <- factors(Pf)
  Pf.feats <- feats(Pf)
  cfg <- yaml.load_file(cfg.in.fname)
  names(Pf.factors) <- str_join(paste(cfg$mapping$dbname, 
                                      cfg$mapping$metadata_tag, sep="."),
                                Pf$factor_cols)
  
  # Construct output filenames
  csv.out.fname <- gsub("\\.yml", sprintf("_pca_%s\\.csv", sha.model) , cfg.in.fname)
  cfg.out.fname <- gsub("\\.yml", sprintf("_pca_%s\\.yml", sha.model) , cfg.in.fname)
  
  # Write csv file
  write.csv(cbind(Pf.factors, Pf.feats), file=csv.out.fname, row.names=F, quote=F)
  
  # Modify cfg file
  cfg$exclude_cols <- NULL
  cfg$mapping$explicit_mappings <- NULL
  cfg$create_treatment_by_combining <- NULL
  cfg$auto_treatment_synonym <- NULL
  cfg$cwd <- NULL
  cfg$feat_start <- as.integer(length(Pf$factor_cols)+1)
  
  # Write cfg file
  file.conn <- file(cfg.out.fname)
  writeLines(as.yaml(cfg), file.conn)
  close(file.conn)
}


#### Old functions for compatibility 


# Plot plates
plot_plates <- function(Xd)
{
  p <- ggplot(Xd, aes(WellCol, WellRow,fill=not_outlier)) + geom_tile() + 
  facet_grid(Plate ~ .) +coord_equal() + xlab("") + ylab("")
  p <- p + geom_text(aes(x=WellCol,y=WellRow, label=RNAi_), 
    size=3, colour="black")
  p <- p + theme_bw()
  return(p)
}

# Plot plates
plot_plates_treatment <- function(Xd)
{
  Xd$RNAi_ <- droplevels(Xd$RNAi_)
  p <- ggplot(Xd, aes(WellCol, WellRow,fill=RNAi_)) + geom_tile() + 
  facet_grid(Plate ~ .) +coord_equal() + xlab("") + ylab("")
  # p <- p + geom_text(aes(x=WellCol,y=WellRow, label=not_outlier), 
  # data=subset(Xd, not_outlier == FALSE), size=3, colour="black")
  p <- p + geom_tile(aes(WellCol,y=WellRow),
    data=subset(Xd, not_outlier == FALSE), 
    color="black",size=0, fill="black", alpha=0.3)
  p <- p + theme_bw()
  return(p)
}

plot_plates_label <- function(Xd, lab)
{
  Xd$lab <- lab
  p <- ggplot(Xd, aes(WellCol, WellRow,fill=lab)) + geom_tile() + 
  facet_grid(Plate ~ .) +coord_equal() + xlab("") + ylab("")
  #  p <- p + geom_text(aes(x=WellCol,y=WellRow, label=not_outlier), 
  # data=subset(Xd, not_outlier == FALSE), size=3, colour="black")
  p <- p + geom_tile(aes(WellCol,y=WellRow),
    data=subset(Xd, not_outlier == FALSE), 
    color="black",size=1, fill="black", alpha=0)
  p <- p + theme_bw()
  return(p)
}


