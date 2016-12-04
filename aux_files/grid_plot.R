my.image <- function(xrange, yrange, figData,  zlim, col, na.color='black', outside.color='black', axes = FALSE, xlab, ylab) {
  newz.na <- zlim[2]+(zlim[2]-zlim[1])/length(col) # new z for NA
  newz.outside <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # new z for values outside zlim
  
  figData[which(is.na(figData) | is.nan(figData))] <- newz.na # we affect newz.outside
  figData[which(figData<zlim[1] | figData>zlim[2])] <- newz.outside # same for newz.na
  
  zlim[2] <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # we finally extend the z limits to include the two new values 
  col <- c(col, na.color, outside.color) # we construct the new color range by including: na.color and outside.color
  image(xrange, yrange, figData,  zlim=zlim, col=col, axes=axes, xlab = xlab, ylab = ylab) # we finally call image(...)
}

myImagePlot <- function(x, min.col, max.col, no.color.scale, ...) {
  if (max.col == 0) {
    max.col <- 1
  }
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(1.0,0.3,length=256),  # Red
                    seq(1.0,0.3,length=256),  # Green
                    seq(1.0,0.6,length=256))  # Blue
  ColorLevels <- seq(min.col, max.col, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,14,2.5,2))
  my.image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
           ylab="", axes=FALSE, zlim=c(min.col,max.col))
  if( !is.null(title) ){
    title(main=title, cex.main = 0.6)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.5)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1.5)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  if (!no.color.scale) {
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
  }
  layout(1)
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}
