count.features.in.categories <- function(feat.names) {
  catx <- c("DNA", "RNA", "Mito", "ER", "AGP")
  caty <- c("Texture", "Intensity", "RadialDistribution")
  catz <- c("Nuclei", "Cytoplasm", "Cells")
  nrm <- data.frame(matrix(rep(0.0, length(catx) * length(caty)), c(length(caty), length(catx))))
  nrm.v <- data.frame(matrix(rep(0.0, 1 * length(catz)), c(1, length(catz))))
  colnames(nrm) <- catx
  rownames(nrm) <- caty
  colnames(nrm.v) <- catz
  rownames(nrm.v) <- "AreaShape"
  for (cx in catx) {
    for (cy in caty) {
      nrm[cy, cx] <- nrm[cy, cx] + length(which(str_detect(feat.names, cy) & str_detect(feat.names, cx)))
    }
  }
  for (cz in catz) {
    nrm.v["AreaShape", cz] <- nrm.v["AreaShape", cz] + length(which(str_detect(feat.names, cz) & str_detect(feat.names, "AreaShape")))
  }
  
  lst.output <- list(normal = nrm, area.shape = nrm.v)
}

print.latex.preamble <- function(file.name, confidential.watermark, tag.pdf) {
  sink(file = file.name, append = FALSE)
  cat("\\documentclass{article}\n")
  cat("\\usepackage{graphicx}\n")
  cat("\\usepackage{array}\n")
  cat("\\usepackage{longtable}\n")
  cat("\\usepackage{geometry}\n")
  cat("\\usepackage{subfig}\n")
  cat("\\usepackage{caption}\n")
  cat("\\usepackage{color}\n")
  cat("\\usepackage{xcolor}\n")
  cat("\\usepackage{pbox}\n")
  cat("\\usepackage{multirow}\n")
  cat("\\usepackage{hyperref}\n")
  cat("\\usepackage{xcolor,colortbl}\n")
  cat("\\hypersetup{colorlinks=false, linkbordercolor=red, pdfborderstyle={/S/U/W 1}}")
  cat("\\captionsetup[subfloat]{labelformat=empty}\n")
  cat("\\geometry{b0paper, margin=0.5in}\n")
  if (confidential.watermark) {
    cat("\\usepackage[printwatermark]{xwatermark}")
    cat("\\usepackage{tikz}")
    cat("\\usepackage{lipsum}")
    cat("\\newsavebox\\mybox")
    cat("\\savebox\\mybox{\\tikz[color=red,opacity=0.1]\\node{\\Huge \\bf CONFIDENTIAL};} \\newwatermark*[ allpages, angle=45, scale=3, xpos=-0, ypos=0]{\\usebox\\mybox}")
    cat("\\newsavebox\\myboxa")
    cat("\\savebox\\myboxa{\\tikz[color=red,opacity=0.1]\\node{\\Huge \\bf CONFIDENTIAL};} \\newwatermark*[ allpages, angle=45, scale=3, xpos=-0, ypos=450]{\\usebox\\mybox}")
    cat("\\newsavebox\\myboxb")
    cat("\\savebox\\myboxb{\\tikz[color=red,opacity=0.1]\\node{\\Huge \\bf CONFIDENTIAL};} \\newwatermark*[ allpages, angle=45, scale=3, xpos=-0, ypos=-450]{\\usebox\\mybox}")
  }
  cat("\\begin{document}\n {\\huge \n")
  cat(sprintf("\\begin{center} {\\bf \\Huge \\color{red} %s} \\end{center} \n", tag.pdf))
  sink()
}

grid.distinctive.feature.summary <- function(Px.zscored, feat.names, metadata.treatment, treatments, significance.ratio) {
  no.features.in.each.category <- count.features.in.categories(feat.names)
  catx <- c("DNA", "RNA", "Mito", "ER", "AGP")
  caty <- c("Texture", "Intensity", "RadialDistribution")
  catz <- c("Nuclei", "Cytoplasm", "Cells")
  ic <- data.frame(matrix(rep(0, length(catx) * length(caty)), c(length(caty), length(catx))))
  colnames(ic) <- catx
  rownames(ic) <- caty
  feat.name.in.org <- feat.names
  m.pos <- data.frame(score = rep(0, length(feat.name.in.org)))
  rownames(m.pos) <- feat.name.in.org
  lstx.tmp <- c()
  
  for (treatment in treatments) {
    trts.indx <- which(metadata.treatment == treatment)
    for (trts.i in trts.indx) {
      match.strong <- which((Px.zscored[trts.i, ] > significance.ratio) | (Px.zscored[trts.i, ] < -significance.ratio))
      match.pos <- which((Px.zscored[trts.i, ] > significance.ratio))
      match.neg <- which((Px.zscored[trts.i, ] < -significance.ratio))
      mt.tmp <- data.frame(score = rep(0, length(feat.name.in.org)))
      mt.tmp[match.pos,] <- as.vector(as.matrix(Px.zscored[trts.i, match.pos]))
      mt.tmp[match.neg,] <- as.vector(as.matrix(Px.zscored[trts.i, match.neg]))
      lstx.tmp <- cbind(lstx.tmp, as.vector(as.matrix(mt.tmp)))
    }
  }
  m.pos$score <- apply(lstx.tmp, 1, median)
  m.pos[which(is.infinite(m.pos[,"score"])), "score"] <- 0
  
  for (cx in catx) {
    for (cy in caty) {
      cx.cy.indx <- feat.name.in.org[which(str_detect(feat.name.in.org, cx) & str_detect(feat.name.in.org, cy))]
      if (length(cx.cy.indx) > 0) {
        u.tmp <- m.pos[cx.cy.indx, "score"]
        u.tmp <- u.tmp[which(!is.infinite(abs(u.tmp)))]
        ic[cy, cx] <- max(abs(u.tmp))
      }
    }
  }
  
  ic[no.features.in.each.category$normal == 0] <- NA
  
  ic.v <- data.frame(matrix(rep(0.0, 1 * length(catz)), c(1, length(catz))))
  colnames(ic.v) <- catz
  rownames(ic.v) <- "AreaShape"
  for (cz in catz) {
    cz.indx <- feat.name.in.org[which(str_detect(feat.name.in.org, "AreaShape") & str_detect(feat.name.in.org, cz))]
    u.tmp <- m.pos[cz.indx, "score"]
    u.tmp <- u.tmp[which(!is.infinite(abs(u.tmp)))]
    ic.v["AreaShape", cz] <- max(abs(u.tmp))
  }
  
  ic.v[no.features.in.each.category$area.shape == 0] <- NA
  
  return(list(normal = ic, area.shape = ic.v, feat.score = m.pos))
}

feats.mds.cluster <- function(z.score, cex.min.thr, cex.max.thr, cex.unit, max.num.of.feats, file.name, cls.fitd.org, third.color = F, second.z.score = NULL, inacitvity.thr = 0.3, reverse.match = F) {
  m.pos <- z.score
  m.pos.clust <- second.z.score
  cex.unit.mds.main <- cex.unit * 8/as.vector(quantile(abs(as.matrix(m.pos$score[m.pos$score != 0])), ((length(which(m.pos$score != 0)) - min(5, length(which(m.pos$score != 0)) - 1))/length(which(m.pos$score != 0)))))
  mds.thr.main <- max(quantile(abs(m.pos$score), 1 - max.num.of.feats/length(m.pos$score)), cex.min.thr / cex.unit.mds.main )
  keep.indx.f.cl <- which(abs(m.pos$score) > mds.thr.main)
  match.ta <- keep.indx.f.cl
  cex.vec <- abs(m.pos$score) * cex.unit.mds.main + 0.02
  cex.vec[cex.vec > cex.max.thr] <- cex.max.thr
  col.vec <- ifelse(m.pos$score > 0, "blue", "red")
  if (third.color) {
    if (reverse.match) {
      opposite.z.score <- which(-m.pos$score * m.pos.clust$score < inacitvity.thr)
      col.vec[opposite.z.score] <- "black"
    } else {
      opposite.z.score <- which(m.pos$score * m.pos.clust$score < inacitvity.thr)
      col.vec[opposite.z.score] <- "black"
    }
  }
  lbls.fea <- rownames(m.pos)
  if (length(keep.indx.f.cl) >= 2) {
    textplot(cls.fitd.org[keep.indx.f.cl, 1], cls.fitd.org[keep.indx.f.cl, 2], lbls.fea[keep.indx.f.cl], xlim = c(min(cls.fitd.org[keep.indx.f.cl, 1]), max(cls.fitd.org[keep.indx.f.cl, 1]))*1.1, ylim = c(min(cls.fitd.org[keep.indx.f.cl, 2]), max(cls.fitd.org[keep.indx.f.cl, 2]))*1.1, cex = cex.vec[keep.indx.f.cl], show.lines = F, col = col.vec[keep.indx.f.cl], font = 2, xaxt='n', yaxt='n', xlab = "", ylab = "")
  } else if (length(keep.indx.f.cl) == 1) {
    plot.new()
    plot(0, font = 2, xaxt="n", yaxt="n", xlab = "", ylab = "", col = "white", xlim = c(min(cls.fitd.org[keep.indx.f.cl, 1]), max(cls.fitd.org[keep.indx.f.cl, 1]))*1.1, ylim = c(min(cls.fitd.org[keep.indx.f.cl, 2]), max(cls.fitd.org[keep.indx.f.cl, 2]))*1.1)
    text(cls.fitd.org[keep.indx.f.cl, 1], cls.fitd.org[keep.indx.f.cl, 2], lbls.fea[keep.indx.f.cl], xlim = c(min(cls.fitd.org[keep.indx.f.cl, 1]), max(cls.fitd.org[keep.indx.f.cl, 1]))*1.1, ylim = c(min(cls.fitd.org[keep.indx.f.cl, 2]), max(cls.fitd.org[keep.indx.f.cl, 2]))*1.1, cex = cex.vec[keep.indx.f.cl], col = col.vec[keep.indx.f.cl], font = 2, xaxt="n", yaxt="n", xlab = "", ylab = "")
  } else {
    plot.new()   
  }
  feat.mds.file <- file.name
  dev.print(device = pdf, feat.mds.file)
}

grid.common.distinctive.feature.summary <- function(Px.compound, Px.orf, compound.metadata, orf.metadata, compound.id, orf.ids, sig.ratio, reverse.match) {
  catx <- c("DNA", "RNA", "Mito", "ER", "AGP")
  caty <- c("Texture", "Intensity", "RadialDistribution")
  catz <- c("Nuclei", "Cytoplasm", "Cells")
  
  no.features.in.each.category <- count.features.in.categories(colnames(Px.compound))
  u.comp <- grid.distinctive.feature.summary(Px.compound, colnames(Px.compound), compound.metadata, compound.id, sig.ratio) 
  u.orf <- grid.distinctive.feature.summary(Px.orf, colnames(Px.orf), orf.metadata, orf.ids, sig.ratio) 
  m.pos <- u.comp$feat.score
  m.pos.orf <- u.orf$feat.score
  feat.names <- rownames(m.pos)
  if (reverse.match) {
    m.pos.orf <- m.pos.orf * -1
  }
  
  ic <- data.frame(matrix(rep(0, length(catx) * length(caty)), c(length(caty), length(catx))))
  colnames(ic) <- catx
  rownames(ic) <- caty
  feats.matching <- feat.names[which((m.pos[feat.names, 1] >= sig.ratio & m.pos.orf[feat.names, 1] >= sig.ratio) | (m.pos[feat.names, 1] <= -sig.ratio & m.pos.orf[feat.names, 1] <= -sig.ratio))]
  
  for (cx in catx) {
    for (cy in caty) {
      cx.cy.indx <- feats.matching[str_detect(feats.matching, cx) & str_detect(feats.matching, cy)]
      avls <- abs(as.vector(m.pos[cx.cy.indx, 1]))
      if (length(avls) > 0) {
        ic[cy, cx] <- max(avls)
      }
    }
  }
  
  nrm <- no.features.in.each.category$normal
  nrm.v <- no.features.in.each.category$area.shape
  
  ic[nrm == 0] <- NA
  ic.v <- data.frame(matrix(rep(0, 1 * length(catz)), c(1, length(catz))))
  colnames(ic.v) <- catz
  rownames(ic.v) <- "AreaShape"
  
  for (cz in catz) {
    tmp.mt <- c()
    cz.indx <- feats.matching[which(str_detect(feats.matching, "AreaShape") & str_detect(feats.matching, cz))]
    avls <- abs(as.vector(m.pos[cz.indx, 1]))
    if (length(avls) > 0) {
      ic.v["AreaShape", cz] <- max(avls)
    }
  }
  ic.v[nrm.v == 0] <- NA
  
  return(list(normal = ic, area.shape = ic.v, feat.score = m.pos))
}
