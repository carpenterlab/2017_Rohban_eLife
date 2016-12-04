library(plyr)
library(stringr)
library(dplyr)
library(EBImage)
library(tiff)
library(png)
library(DBI)
library(RMySQL)

src_mysql_custom <- function(dbname, host = NULL, port = 0L, user = "root", password = "", ...) {
  con <- DBI::dbConnect(RMySQL::MySQL(), dbname = dbname, host = host,
                        port = port, username = user, password = password, ...)
  info <- DBI::dbGetInfo(con)
  dplyr::src_sql("mysql", con, info = info)
}

## path example /cbnt/cbimageX/HCS/xiaoyunwu/taoe005-u2os-72h-cp-a-au00044858/au00044882/2013-08-09/41757/taoe005-u2os-72h-cp-a-au00044858_p05_s6_w290ae2577-e7e4-41a2-8061-2483f755e5bc.tif
get.image.from.server.tiff <- function(path) {
  #s <- sprintf("ssh mrohban@pirita.broadinstitute.org \"cat %s\" > tmp.tif", path)
  #system(s)
  img.pl <- ""
  fl.name <- ""
  
  for (pl in c("41744", "41754", "41755", "41756", "41757")) {
    v <- stringr::str_split(path, sprintf("/%s/", pl))
    if (length(v[[1]]) > 1) {
      fl.name <- v[[1]][2]
      img.pl <- pl
      break
    }
  }
  
  s <- sprintf("../images/%s/%s", img.pl, fl.name)
  x <- readTIFF(s)
  
  #system("rm tmp.tif")
  return(x)
}

get.image.from.server.png <- function(path) {
  s <- sprintf("ssh mrohban@pirita.broadinstitute.org \"cat %s\" > tmp.png", path)
  system(s)
  x <- readPNG("tmp.png")
  system("rm tmp.png")
  return(x)
}

get.image.path <- function(image.number, channel, dbname, host, user, password, port) {
  channel <- str_replace_all(channel, "ER", "ERSyto")
  channel <- str_replace_all(channel, "RNA", "ERSytoBleed")
  channel <- str_replace_all(channel, "DNA", "Hoechst")
  channel <- str_replace_all(channel, "AGP", "PhGolgi")
  
  db <- src_mysql_custom(dbname = dbname, host = host, user = user, password = password, port = port)

  qu <- sprintf("select Image_URL_Orig%s from SIGMA2_Pilot_2013_10_11_Analysis_Per_Image where ImageNumber = %d", channel, image.number)
  img.db <- 
    db %>%
    dplyr::tbl(dplyr::sql(qu), quiet = T) %>%
    dplyr::collect(.,quiet = T)
  
  path <- (img.db %>% as.matrix() %>% as.vector())
  path <- str_replace(path, "file:", "")
  DBI::dbDisconnect(db$con)
  return(path)
}

get.outline.path <- function(image.number, dbname, host, user, password, port) {
  db <- src_mysql_custom(dbname = dbname, host = host, user = user, password = password, port = port)
  
  qu <- sprintf("select Image_PathName_CellOutlines, Image_FileName_CellOutlines from SIGMA2_Pilot_2013_10_11_Analysis_Per_Image where ImageNumber = %d", image.number)
  img.db <- 
    db %>%
    dplyr::tbl(dplyr::sql(qu), quiet = T) %>%
    dplyr::collect(.,quiet = T)
  path <- (img.db %>% as.matrix() %>% as.vector())
  DBI::dbDisconnect(db$con)
  return(str_join(path, collapse = "/"))
}

get.image <- function(image.number, channel, dbname, host, user, password, port) {
  path <- get.image.path(image.number = image.number, channel = channel, dbname, host, user, password, port)
  return(get.image.from.server.tiff(path = path))
}

get.image.by.plate.well <- function(plate, well, channel, dbname, host, user, password, port) {
  db <- src_mysql_custom(dbname = dbname, host = host, user = user, password = password, port = port)
  
  qu <- sprintf("select ImageNumber from SIGMA2_Pilot_2013_10_11_Analysis_Per_Image where Image_Metadata_Plate = '%s' and Image_Metadata_Well = '%s' and Image_Metadata_Site = 5", plate, well)
  img.db <- 
    db %>%
    dplyr::tbl(dplyr::sql(qu), quiet = T) %>%
    dplyr::collect(.,quiet = T)
  im.number <- (img.db %>% as.matrix() %>% as.vector() %>% head(., 1))
  DBI::dbDisconnect(db$con)
  return(get.image(image.number = im.number, channel = channel, dbname, host, user, password, port))
}

get.image.by.plate.treatment <- function(plate, treatment, channel, dbname, host, user, password, port) {
  db <- src_mysql_custom(dbname = dbname, host = host, user = user, password = password, port = port)
  
  qu <- sprintf("select ImageNumber from SIGMA2_Pilot_2013_10_11_Analysis_Per_Image where Image_Metadata_Plate = '%s' and Image_Metadata_GeneSymbol = '%s' and Image_Metadata_AlleleDesc = '%s' and Image_Metadata_Site = 5", plate, str_split(treatment, "_")[[1]][1], str_split(treatment, "_")[[1]][2])
  img.db <- 
    db %>%
    dplyr::tbl(dplyr::sql(qu), quiet = T) %>%
    dplyr::collect(.,quiet = T)
  im.number <- (img.db %>% as.matrix() %>% as.vector() %>% head(., 1))
  DBI::dbDisconnect(db$con)
  return(get.image(image.number = im.number, channel = channel, dbname, host, user, password, port))
}

get.outline <- function(image.number, dbname, host, user, password, port) {
  path <- get.outline.path(image.number, dbname, host, user, password, port)
  I <- get.image.from.server.png(path)
  return(I)
}

which.nz <- function(I, I.explore) {
  lst <- which(I == 1 & I.explore == 0, arr.ind = T)
  return(lst)
}

fill.image <- function(I.mask, loc) {
  i <- loc[2]  
  j <- loc[1]

  I <- I.mask * 0
  n <- NROW(I)
  m <- NCOL(I)
  I.explore <- I.mask * 0
  I[i, j] <- 1

  for (iter in 1:max(m, n)) {
    lst <- which.nz(I, I.explore)  
    I.last <- I
    
    for (k in 1:NROW(lst)) {
      i <- lst[k, 1]
      j <- lst[k, 2]
      I.explore[i, j] <- 1
      
      if (I.mask[i, j] != 1) {
        if (i + 1 <= n) {
          I[i+1, j] <- 1  
        }
        if (i - 1 >= 1) {
          I[i-1, j] <- 1  
        }
        if (j + 1 <= m) {
          I[i, j+1] <- 1  
        }
        if (j - 1 >= 1) {
          I[i, j-1] <- 1  
        }
      }
    }
    if (all(I.last == I)) {
      break
    }
  }
  return(I)
}

get.cell.mask <- function(image.number, object.number, dbname, host, user, password, port) {
  return(fill.image(get.outline(image.number, dbname, host, user, password, port), get.cell.location(image.number, object.number, dbname, host, user, password, port)))
}

get.cell.location <- function(ImageNumber, ObjectNumber, dbname, host, user, password, port) {
  db <- src_mysql_custom(dbname = dbname, host = host, user = user, password = password, port = port)
  
  qu <- sprintf("select O.Nuclei_Location_Center_X, O.Nuclei_Location_Center_Y from SIGMA2_Pilot_2013_10_11_Analysis_Per_Image as I
                inner join SIGMA2_Pilot_2013_10_11_Analysis_Per_Nuclei as O on I.ImageNumber = O.ImageNumber where I.ImageNumber = %d and O.Nuclei_Number_Object_Number = %d", ImageNumber, ObjectNumber)
  img.db <- 
    db %>%
    dplyr::tbl(dplyr::sql(qu), quiet = T) %>%
    dplyr::collect(.,quiet = T) 
  
  DBI::dbDisconnect(db$con)
  return(img.db %>% as.matrix() %>% as.vector())
}

get.cell.image <- function(ImageNumber, ObjectNumber, WindowSize, dbname, host, user, password, port) {
  chnls <- c("ER", "Mito", "DNA", "AGP", "RNA")
  Ix <- list(ER = c(), Mito = c(), DNA = c(), AGP = c(), RNA = c())
  for (channel in chnls) {
    I <- get.image(image.number = ImageNumber, channel = channel, dbname, host, user, password, port)
    loc <- get.cell.location(ImageNumber = ImageNumber, ObjectNumber = ObjectNumber, dbname, host, user, password, port)
    n <- NROW(I)
    m <- NCOL(I)
    
    xlim <- c(loc[1]-WindowSize/2, loc[1]+WindowSize/2) %>% round
    ylim <- c(loc[2]-WindowSize/2, loc[2]+WindowSize/2) %>% round
    xlim[2] <- min(xlim[2], m)
    xlim[1] <- max(xlim[1], 0)
    ylim[2] <- min(ylim[2], n)
    ylim[1] <- max(ylim[1], 0)
    #print(xlim)
    #print(ylim)
    Is <- I[ylim[1]:ylim[2], xlim[1]:xlim[2]]
    Ix[[channel]] <- Is
  }
  return(Ix)
}

get.object.image <- function(ImageNumber, ObjectNumber, dbname, host, user, password, port) {
  options(warn=-1)
  I.list <- list(ER = c(), RNA = c(), DNA = c(), Mito = c(), AGP = c(), mask = c(), outline = c())
  for (channel in c("ER", "RNA", "DNA", "Mito", "AGP")) {
    I1 <- get.image(ImageNumber, channel, dbname, host, user, password, port)
    I.list[[channel]] <- I1
  }
    
  I.mask <- get.cell.mask(ImageNumber, ObjectNumber, dbname, host, user, password, port)
  I.outline <- get.outline(ImageNumber, dbname, host, user, password, port)
  loc <- get.cell.location(ImageNumber = ImageNumber, ObjectNumber = ObjectNumber, dbname, host, user, password, port)
  n <- NROW(I.list$ER)
  m <- NCOL(I.list$ER)
  for (channel in c("ER", "RNA", "DNA", "Mito", "AGP")) {
    I <- I.list[[channel]]
    I <- I * I.mask
    sx <- apply(I, 2, sum)
    sy <- apply(I, 1, sum)
    x1 <- min(which(sx > 0))
    x2 <- max(which(sx > 0))
    y1 <- min(which(sy > 0))
    y2 <- max(which(sy > 0))
    xlim <- c(x1, x2) %>% round
    ylim <- c(y1, y2) %>% round
    xlim[2] <- min(xlim[2], m)
    xlim[1] <- max(xlim[1], 0)
    ylim[2] <- min(ylim[2], n)
    ylim[1] <- max(ylim[1], 0)
    
    Is <- I[ylim[1]:ylim[2], xlim[1]:xlim[2]]
    I.list[[channel]] <- Is
  }
  I.list$mask <- I.mask[ylim[1]:ylim[2], xlim[1]:xlim[2]]
  I.list$outline <- I.outline[ylim[1]:ylim[2], xlim[1]:xlim[2]] * I.list$mask

  return(I.list)
  options(warn=0)
}

readCellImage <- function(plate, well, channel, dbname, host, user, password, port) {
  x <- get.image.by.plate.well(plate = plate, well = well, channel = channel, dbname = dbname, host = host, user = user, password = password, port = port)
  return(x)
}

readCellImageGene <- function(plate, gene, channel, Pf, dbname, host, user, password, port) {
  x <- get.image.by.plate.treatment(plate = plate, treatment = gene, channel = channel, dbname = dbname, host = host, user = user, password = password, port = port)
  return(x)
}


