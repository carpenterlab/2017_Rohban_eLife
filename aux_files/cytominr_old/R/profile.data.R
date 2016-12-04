library(yaml)
library(plyr)
library(caret)
library(testthat)
library(sampling)
library(stringr)
library(testthat)
library(digest)
library(gridExtra)

#-------------------------------------------------------------------------------
path                  <- function(obj, ...) UseMethod("path", obj)
prune.feats.nzv       <- function(obj, ...) UseMethod("prune.feats.nzv", obj)
prune.feats.nzv.apply <- function(obj, ...) UseMethod("prune.feats.nzv.apply", obj)
rename.feats          <- function(obj, ...) UseMethod("rename.feats", obj)
log.feats             <- function(obj, ...) UseMethod("log.feats", obj)
prune.feats           <- function(obj, ...) UseMethod("prune.feats", obj)
strata1               <- function(obj, ...) UseMethod("strata1", obj)
create.plate.pos      <- function(obj, ...) UseMethod("create.plate.pos", obj)
replace.feats         <- function(obj, ...) UseMethod("replace.feats", obj)
prune.rows            <- function(obj, ...) UseMethod("prune.rows", obj)
prune.ctrl.wells      <- function(obj, ...) UseMethod("prune.ctrl.wells", obj)
prune.outliers        <- function(obj, ...) UseMethod("prune.outliers", obj)
score.cells           <- function(obj, ...) UseMethod("score.cells", obj)
factors               <- function(obj, ...) UseMethod("factors", obj)
feats                 <- function(obj, ...) UseMethod("feats", obj)
plate.orig            <- function(obj, ...) UseMethod("plate.orig", obj)
plate                 <- function(obj, ...) UseMethod("plate", obj)
well                  <- function(obj, ...) UseMethod("well", obj)
pca                   <- function(obj, ...) UseMethod("pca", obj)
pca.apply             <- function(obj, ...) UseMethod("pca.apply", obj)
wtavg                 <- function(obj, ...) UseMethod("wtavg", obj)
liftcurve             <- function(obj, ...) UseMethod("liftcurve", obj)
read.outlier.data     <- function(obj, ...) UseMethod("read.outlier.data", obj)
plot.plates.treatment <- function(obj, ...) UseMethod("plot.plates.treatment", obj)
rbind2                <- function(obj, ...) UseMethod("rbind2", obj)
merge2                <- function(obj, ...) UseMethod("merge2", obj)
loo.liftcurve         <- function(obj, ...) UseMethod("loo.liftcurve", obj)
pick.best.subset      <- function(obj, ...) UseMethod("pick.best.subset", obj)
remove.singletons     <- function(obj, ...) UseMethod("remove.singletons", obj)
randperm              <- function(obj, ...) UseMethod("randperm", obj)
sigstrength           <- function(obj, ...) UseMethod("sigstrength", obj)
selfsim               <- function(obj, ...) UseMethod("selfsim", obj)
qualitymeas           <- function(obj, ...) UseMethod("qualitymeas", obj)
apply.well.position   <- function(obj, ...) UseMethod("apply.well.position", obj)
medpolish.plate       <- function(obj, ...) UseMethod("medpolish.plate", obj)
dist1                 <- function(obj, ...) UseMethod("dist1", obj)
prune.to.identical.freq.for.labels <- 
  function(obj, ...) UseMethod("prune.to.identical.freq.for.labels", obj)
#-------------------------------------------------------------------------------

factors.profile.data    <- function(obj, ...) obj$data[obj$factor_cols]
feats.profile.data      <- function(obj, ...) obj$data[obj$feat_cols]

plate.orig.profile.data <- function(obj, ...) factor(obj$data$Plate,
                                                     levels=obj$cfg$plate_names_abbrev,
                                                     labels=obj$cfg$plate_names)
plate.profile.data      <- function(obj, ...) obj$data$Plate
well.profile.data       <- function(obj, ...) obj$data$Well
path.profile.data       <- function(obj, f,...) file.path(obj$cfg$cwd,f)

randperm.profile.data   <- function(obj, each_col=FALSE, ...) {

  ## Reps can cause problems so this has been removed
  ## if (nreps > 1) {
  ##   obj$data <- do.call("rbind", replicate(nreps, obj$data, simplify = FALSE))
  ## }
  if (!each_col) {
    obj$data[,obj$factor_cols] <- obj$data[sample(NROW(obj$data)),
                                           obj$factor_cols]
  } else {
    for (i in obj$factor_cols)
      obj$data[,i] <- obj$data[sample(NROW(obj$data)),i]
  }
  obj
}

profile.data.dummy <- function(n=1000,d=10) {
  data <- data.frame(data.frame(ceiling(2*runif(n))),
                     as.data.frame(matrix(rnorm(n*d),n,d) ))
  obj <- list(data=data, cfg=NULL)
  class(obj) <- "profile.data"
  obj$factor_cols <- 1
  obj$feat_cols <- 2:NCOL(data)
  obj  
}

profile.data.init <- function() {
  obj <- list(data=NULL, cfg=NULL)
  class(obj) <- "profile.data"
  obj
}

handle_version_pre1 <- function(obj) {
  obj$data$RNAi <- factor(obj$data$RNAi)
  obj$data$Plate <- factor(obj$data$Plate)
  
# Remove No-treatment
  if (!is.null(obj$cfg$exclude_nt) && obj$cfg$exclude_nt)
    obj$data <- obj$data[obj$data$RNAi != obj$cfg$no_treat_name,]
  
  obj$data <- droplevels(obj$data)
  
#------- Format the treatment names
  if (!is.null(obj$cfg$control_name)) {
    ordering <- c(obj$cfg$treatment_names, obj$cfg$control_name)
    ordering_abbrev <- c(obj$cfg$treatment_names_abbrev,
                         obj$cfg$control_name_abbrev)
    ordering_grouped_abbrev <- c(obj$cfg$treatment_names_grouped_abbrev, 
                                 obj$cfg$control_name_abbrev)
  } else {
    ordering <- obj$cfg$treatment_names
    ordering_abbrev <- obj$cfg$treatment_names_abbrev
    ordering_grouped_abbrev <- obj$cfg$treatment_names_grouped_abbrev
  }
  
  obj$data$RNAi <- factor(obj$data$RNAi, levels=ordering)
  obj$data$RNAi_ <- obj$data$RNAi
  levels(obj$data$RNAi_) <- ordering_abbrev
  obj$data$RNAi_g_ <- obj$data$RNAi
  levels(obj$data$RNAi_g_) <- ordering_grouped_abbrev
  
#------- Format the Plate names
  obj$data$Plate <- factor(obj$data$Plate, levels=obj$cfg$plate_names)
  obj$data$Plate_ <- obj$data$Plate
  levels(obj$data$Plate) <- obj$cfg$plate_names_abbrev
  
#------- Set column names  
  obj$factor_cols <- c("Plate", "Plate_", "Well", "RNAi", "RNAi_", "RNAi_g_")
  obj
}

handle_version_pre2 <- function(obj) {
#------- Rename meta-data column names
  dbname <- obj$cfg$mapping$dbname
  metadata_tag <- obj$cfg$mapping$metadata_tag

#	browser() 
  
  metadata_names <- gsub(paste(dbname, '.', sep=''), '', 
                         gsub(metadata_tag, '',
                              names(obj$data)[grep(metadata_tag,
                                                   names(obj$data))]))
  metadata_names_processed <- c()
  for (s2 in metadata_names) {
    expm <- obj$cfg$mapping$explicit_mappings
    s1 <- if (!is.null(expm[[s2]])) expm[[s2]] else s2
    s3 <- paste(dbname, paste(metadata_tag, s2, sep=''), sep=".")
    names(obj$data)[names(obj$data)==s3] <- s1
    s4 <- s1
    
    lastchar <- str_sub(s1,-1,-1)
# convert to factor unless there is * or | in the name
    cast_as_factor <- (!is.null(obj$cfg$mapping$cast_as_factor)
                       && obj$cfg$mapping$cast_as_factor)
    if (!(lastchar %in% c('*', '|')) && cast_as_factor) {
      obj$data[,s1] <- factor(obj$data[,s1])
    }
    
# convert to numeric if there is | in the name
    if (lastchar == '|') {
      obj$data[,s1] <- as.numeric(as.character(obj$data[,s1]))
# Trim the name
      s4 <- str_sub(s1, 1,nchar(s1)-1)
      names(obj$data)[names(obj$data)==s1] <- s4	
    }

# convert to character if there is * in the name
    if (lastchar == '*') {
      obj$data[,s1] <- as.character(obj$data[,s1])				
# Trim the name
      s4 <- str_sub(s1, 1,nchar(s1)-1)
      names(obj$data)[names(obj$data)==s1] <- s4	
    }
    
# Append X to header names starting with a digit
    if (str_sub(s4, 1,1) %in% as.character(seq(10)-1)) {      
      names(obj$data)[names(obj$data)==s4] <- paste('X', s4, sep='')
      s4 <- paste('X', s4, sep='')
    }
    
    metadata_names_processed <- c(metadata_names_processed, s4)
    
  }
  
  obj$factor_cols <- metadata_names_processed
  
#------- Format plate names
  if (!is.null(obj$cfg$plate_names)) {
    obj$data$Plate                <- factor(obj$data$Plate,
                                            levels=obj$cfg$plate_names)
    obj$data$Plate_               <- obj$data$Plate
    obj$factor_cols               <- c(obj$factor_cols, 'Plate_')
    levels(obj$data$Plate)        <- obj$cfg$plate_names_abbrev
  } else {
    if (!("Plate_" %in% obj$factor_cols))
      obj$data$Plate_             <- obj$data$Plate    
  }
  
#------- Remove No-treatment
  if (!is.null(obj$cfg$exclude_nt) && obj$cfg$exclude_nt) {
    obj$data <- obj$data[obj$data$Treatment != obj$cfg$no_treat_name,]      
    obj$data <- droplevels(obj$data)
  }

  
#------- Create a treatment column by concatenating two columns
  if (!is.null(obj$cfg$create_treatment_by_combining)) {
    cols <- str_split(obj$cfg$create_treatment_by_combining, ',')[[1]]
    obj$data$Treatment <- paste(obj$data[,cols[1]], obj$data[,cols[2]], sep='_')
    obj$factor_cols <- c(obj$factor_cols, 'Treatment')
    
  }

    
#------- Format treatment names
  if (!is.null(obj$cfg$auto_treatment_synonym)) {
    syn <- obj$cfg$auto_treatment_synonym
# TODO: Handle the case that Type does not exist
    df <- unique(obj$data[,c(syn, 'Treatment', 'Type')])
    names(df)[names(df)==syn] <- 'syn'
    df <- df[with(df, order(syn, Treatment)), ]
    df0 <- ddply(df, .(syn), index_factor_levels)
    df <- df0
    df <- df[with(df, order(Type, syn, Treatment)), ]	  
    names(df)[names(df)=='Treatment'] <- 'name_Treatment'
    names(df)[names(df)=='abbrev'] <- 'name_TreatmentAbbrev'
    names(df)[names(df)=='syn'] <- sprintf('name_%s',
                       obj$cfg$auto_treatment_synonym)
    df$Type <- NULL
    obj$cfg$tab_treatment_synonym <- df	  
  }
  
  if (!is.null(obj$cfg$treatment_synonyms)) {
    df <- convert_str_to_dataframe(obj$cfg$treatment_synonym)
    obj$cfg$tab_treatment_synonym <- df
  }
  
  if (!is.null(obj$cfg$auto_treatment_synonym) ||
      !is.null(obj$cfg$treatment_synonyms)) {
    obj$data$Treatment <- factor(obj$data$Treatment, levels = df$name_Treatment)
    df <- obj$cfg$tab_treatment_synonym
    sfnames <- names(df)[!(names(df) %in% c('name_Treatment'))]
    for (sfname in sfnames) 
      obj$data[,gsub('name_', '', sfname)] <-
          create_synonym_factor(obj$data$Treatment, df[,sfname])
    
# do this because the syn is already present in the table
    sfnames <- setdiff(sfnames, c(sprintf('name_%s',
                                          obj$cfg$auto_treatment_synonym)))
    
    obj$factor_cols <- c(obj$factor_cols, gsub('name_', '', sfnames))
    
  }

  
#------- Set column names  
  obj
  
}

index_factor_levels <- function(dd, group='syn', abbrev='abbrev') {
  dd[,abbrev] <- paste(dd[,group], seq(NROW(dd)), sep='.')
  dd
}

create_synonym_factor <- function(src_factor, synonym_levels) {
  fac <- src_factor
  levels(fac) <- synonym_levels
  fac
}

handle_non_microplate <- function(obj) {

  obj$factor_cols <- getd(obj$cfg, 'factor_cols',
                          names(obj$data)[1:(obj$cfg$feat_start-1)])
  
  if (getd(obj$cfg, 'flag_dont_cast_factors', F)) {
    for (f in obj$factor_cols) {
      obj$data[,f] <- factor(obj$data[,f])
    }        
  }
  obj	
}

profile.data <- function(cf) {
  'Construct a profile.data object by reading a yml config file'

# ------------ Setup --------------

# Load the config file
  cfg <- yaml.load_file(cf)
  cfg$cwd <- dirname(cf)
  if (!is.null(cfg$profile_dir)) 
    cfg$cwd <- file.path(cfg$cwd, cfg$profile_dir)
  
  
# Get name of data file (both csv and rda)
  cf_name <- sapply(strsplit(basename(cf),"\\."), 
                    function(x) paste(x[1:(length(x)-1)], collapse=".")) 
  if(is.null(cfg$profile_file))
    cfg$profile_file <- paste(cf_name, "csv", sep=".")
  
  if(is.null(cfg$profile_file_binary))
    cfg$profile_file_binary <- paste(cf_name, "rda", sep=".")
  
  frda <- file.path(cfg$cwd, cfg$profile_file_binary)
  fcsv <- file.path(cfg$cwd, cfg$profile_file)
  
  
# Read binary file, or if it doesn't exist, then csv. 
# Save binary file if it doesn't exist
  
  if (file.exists(frda)) {
    data <- readRDS(frda)
  } else {
    data <- data.frame(read.csv(fcsv, header=TRUE))
    saveRDS(data, file=frda)
  }
# TODO: Need to save a checksum of the csv file to make sure that the version
#  of the rda and the csv match. md5 on MAC or md5sum on Linux would work well
#  for this
  
  
# Create profile.data object
  obj <- list(data=data, cfg=cfg)
  class(obj) <- "profile.data"
  
  
  
# ------------ Get features --------------
  allfeat_cols <- c(obj$cfg$feat_start:length(names(obj$data)))
  exclude_cols <- NULL
  if (!is.null(obj$cfg$exclude_cols)) {
    for (f in obj$cfg$exclude_cols)
      exclude_cols <- c(exclude_cols, grep(f, names(obj$data)))  
  }
  obj$feat_cols <- names(obj$data)[setdiff(allfeat_cols, exclude_cols)]
  
  
#------- Get version of the cfg file
  version <- getd(obj$cfg, 'version', 0)
  obj$cfg$treatment_tag <- if (version < 1) "RNAi_" else "Treatment_"
  
  
# ------------ Get factors --------------
  if (obj$cfg$flag_microplate_dataset) {
    if (version < 1) {
      obj <- handle_version_pre1(obj)
    } else if (version < 2) {
      obj <- handle_version_pre2(obj)
    }
  } else {
    obj <- handle_non_microplate(obj)
  }
  
# ------------ Keep only factors and features --------------
  
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)]
  obj
}

prune.feats.profile.data <- function(obj, feat_cols, keep = FALSE, ...) { 
  'Remove features'
  if (!keep)
    obj$feat_cols <- setdiff(obj$feat_cols, feat_cols)
  else
    obj$feat_cols <- feat_cols
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)] 
  obj 
}



rename.feats.profile.data <- function(obj, new_names,  ...) { 
  'Rename features'

  X <- feats(obj)
  names(X) <- new_names
  obj$feat_cols <- names(X)
  obj$data <- cbind(X, factors(obj))
  obj 

}



prune.feats.nzv.profile.data <- function(obj, get.nzv=FALSE, ...) { 
  'Remove features that have near-zero variance but do this on a per plate basis'
  print('Warning: Removing nzv features by computing nzv on a per-plate basis.')
  plates <- plate(obj) # HACK : This should not be hard-coded
  X <- feats(obj)
  nzv_cols <- list()
  for (plate in levels(plates)) {
    r <- plates == plate
    s <- foreach(i=1:NCOL(X), .combine=cbind) %dopar% if(length(nearZeroVar(X[r,i])) > 0) i
    nzv_cols <- union(nzv_cols, unlist(s))
  } 
# remove list formatting
  nzv_cols <- names(X)[unlist(nzv_cols)] 
# get new set of feats_cols
  obj$feat_cols <- setdiff(obj$feat_cols, nzv_cols) 
# re-generate data frame based on new set of feat_cols
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)] 
  if (get.nzv) {
    nzv_cols
  } else {
    obj 
  }
}

prune.feats.nzv.apply.profile.data <- function(obj, nzv_cols, ...) { 
  if (length(nzv_cols) > 0) {
# remove list formatting
    nzv_cols <- feats(obj)[unlist(nzv_cols)] 
# get new set of feats_cols
    obj$feat_cols <- setdiff(obj$feat_cols, nzv_cols) 
# re-generate data frame based on new set of feat_cols
    obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)] 
  }
  obj
}



strata1.profile.data <- function(obj, stratanames, size, method, ...) {
  # Sort obj by stratanames as required by strata
  obj$data <- obj$data[eval(parse(text=sprintf('with(obj$data, order(%s))',
                                      paste(stratanames, collapse=",")))),]

  st <- strata(obj$data, 
               stratanames=stratanames, 
               size=size, 
               method=method)
  obj$data <- droplevels(getdata(obj$data,st))
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)]
  obj
}

create.plate.pos.profile.data <- function(obj, 
                                          numrows = 16, 
                                          numcols = 24, 
                                          modify.data = FALSE,
                                          ...) {
  'Create Row and Col attributes for plotting plate'

  tmp <- obj$data[,obj$factor_cols]
  
  tmp$WellRow <- sapply(tmp$Well, 
                        function(x) toupper(substr(x,1,1)))
  
  tmp$WellCol <- sapply(tmp$Well,
                        function(x) as.integer(substr(x,2,3) ))
  
  tmp$WellCol <- ordered(tmp$WellCol,
                         levels=1:numcols)
  
  levels(tmp$WellCol) <- 1:numcols
  
    tmp$WellRow <- ordered(tmp$WellRow,
                           levels=rev(LETTERS[seq( from = 1,
                                                   to = numrows )]))
  
  
  if(modify.data) {
    obj$data <- merge(obj$data, tmp)
  } else {
    obj$data_plate_plot <- tmp
  }
  
  obj
}

replace.feats.profile.data <- function(obj, X, ...) {
  'Replace the features with X. Typical use involves first making
      a copy of this object before replacing features'
  obj$data <- data.frame(X,factors(obj))
  obj$feat_cols <- names(X)
  obj
}


log.feats.profile.data <- function(obj, replace=F...) {
  'Compute log of each feature and add it to the set of features'
  X <- feats(obj)
  Xlog <- colwise(function(x) log(x-min(x)+1))(X)
  names(Xlog) <- sapply(names(Xlog), function(x) paste( x, 'log', sep="_"))
  
  if(replace) {
    obj$data <- data.frame(Xlog, factors(obj))
    obj$feat_cols <- c(names(Xlog))
  } else {
    obj$data <- data.frame(X, Xlog, factors(obj))
    obj$feat_cols <- c(names(X), names(Xlog))
  }
  
  obj  
}

prune.rows.profile.data <- function(obj, prune_vec, ...) {
  'Prune rows using prune_vec. Typical use involves first making
      a copy of this object before deleting rows'
  expect_equal(typeof(prune_vec), "logical")
  obj$data <- obj$data[prune_vec,]
  obj$data <- droplevels(obj$data)
  obj
}

prune.ctrl.wells.profile.data <- function(obj, ...) {
  prune_vec <- factors(obj)[,obj$cfg$treatment_tag] != obj$cfg$control_name_abbrev
  prune.rows(obj, prune_vec)
}

prune.outliers.profile.data <- function(obj, not_outlier, ...) {
  'Prune outliers and add a column to data_plate_plate to indicate 
  the outlier well. This code is potentially buggy because it assumes
  that data and data_plate_plot have the same rows until this point,
  that is, no call to pruneRows has been made before this'
  obj <- prune.rows(obj, not_outlier)
  if (!is.null(obj$data_plate_plot))
    obj$data_plate_plot$not_outlier <- not_outlier
  obj
}

pca.profile.data <- function(obj, percent_var, scale=FALSE, debug=FALSE,
                             get.proj=FALSE, ...) {
  v <- prcomp(feats(obj), scale.=scale)
  u <- v$x
  npcs <- which(cumsum(v$sdev^2)/sum(v$sdev^2) > percent_var)[1]-1
  if (get.proj) {
    model=list(v=v, npcs=npcs)
  } else {
    replace.feats(obj, as.data.frame(u[,1:npcs]))
  }
}

pca.apply.profile.data <- function(obj, model, ...) {
  replace.feats(obj, as.data.frame(
    predict(model$v, feats(obj))[,1:model$npcs])
    )
}


remove.singletons.profile.data <- function(obj, fac2_col, pick_best=FALSE, ...) {

  fac <- factors(obj)
  tbl <- table(fac[,fac2_col])
  singletons <- names(tbl[tbl==1])

  print('Removing singletons:')
  print(singletons)

  obj <- prune.rows(obj, !(fac[,fac2_col] %in% singletons))

  obj
}



pick.best.subset.profile.data <- function(obj, fac1_col, fac2_col, fac3_col,
                                          n = 2, ...) {

  obj <- remove.singletons(obj, fac2_col)

  fac <- factors(obj)

  # check if rows are unique wrt fac_col1
  expect_true(setequal(unique(fac[, fac1_col]), fac[, fac1_col]))

  data <- foreach (i = unique(fac[,fac2_col]), .combine=rbind) %dopar%  {
    P0 <- prune.rows(obj, fac[,fac2_col] == i)
    P0 <- prune.rows(P0, !duplicated(factors(P0)[, fac3_col]))

    if (NROW(P0$data) <= n) {
      P0$data
    } else {
      Dx <- as.matrix(dist(feats(P0), diag = T, upper = T))
      diag(Dx) <- Inf
      minij <- unname(which(Dx == min(Dx), arr.ind = TRUE)[1,])
      if (n > 2) {
        max_flag <- FALSE
        for (j in seq(n-2)) {
          cs <- colSums(Dx[minij,])
          expect_true(!is.infinite(min(cs)))
          minij <- c(minij, unname(which.min(cs)))
        }
      }
      P0$data <- P0$data[minij,]
      P0$data
    }
  }
 

  obj$data <- data
  obj
}

liftcurve.helper <- function(Dx, 
                             facbin0, 
                             fac, 
                             exclude_first_col=FALSE,
                             compute_closest=FALSE) {
  # Dx     : n x n   matrix where (i,j) is distance between v_i and v_j
  # facbin0: n x n   matrix where (i,j) is the True if i and j belong to the same class
  
  n <- NROW(facbin0)
  
  # Sort rows of facbin0 based on distance matrix
  facbin1 <- laply(seq(n), function(i) 1.0*facbin0[i,order(Dx[i,])],
                   .parallel=T)
  if (exclude_first_col) {
    facbin1 <- facbin1[,2:NCOL(facbin1)]
  }

  if (compute_closest) {
    # Find rank of closest matching vector and normalize to one
    return (aaply(facbin1, 1, function(x) match(T,x)/length(x),  .parallel=T,
                  .progress = "text"))
  }
  
  # Generate cumulative sum along each row of facbin1
  facbin2 <- aaply(facbin1, 1, function(x) cumsum(x), 
                   .parallel=T)
  
  # Exclude cases where no matches exist
  flag <- facbin2[,NCOL(facbin2)] !=0
  facbin2 <- facbin2[flag, ]
  # Make a copy of fac with pruned rows based on the above
  facp <- fac[flag,]

  # Normalize rows of facbin2
  facbin3 <- aaply(facbin2, 1, function(x) x / tail(x,1), 
                   .parallel=T)  

  list(facbin1=facbin1,
       facbin2=facbin2,
       facbin3=facbin3,
       facp=facp,
       lc=colMeans(facbin3),
       matdim=dim(facbin3))

}

  
liftcurve.profile.data <- function(obj, 
                                   fac_col, 
                                   fac_col_exclude=NULL,
                                   select.cols = NULL,
                                   distfunc='euc', 
                                   compute_closest = FALSE,
                                   ...) {
  
  expect_false(is.null(fac_col))
  expect_true(all(fac_col %in% obj$factor_cols))
  
  if(!is.null(fac_col_exclude)) {
    expect_true(all(fac_col_exclude %in% obj$factor_cols))
  }

  if(!is.null(select.cols)) {
    expect_true(select.cols %in% obj$factor_cols)
    expect_true(length(select.cols)==1)
  }

  X <- feats(obj)
  if (distfunc=='euc') {
    Dx <- as.matrix(dist(X, diag = T, upper = T))
  } else if (distfunc=="cos") {
    message(sprintf('Using %s dist.', distfunc))
    Dx <- cos.dist(X)    
  } else if (distfunc %in% c("pearson", "kendall","spearman")) {
    Dx <- 1-cor(t(X), method=distfunc)
  }
  
  fac <- factors(obj)
  # convert to a string that can be easily compared
  fac.l <- apply(fac[fac_col], 1, function(x) digest(x))
  facbin0 <- outer(fac.l, fac.l, "==")
  
  if (!is.null(fac_col_exclude)) {
    fac.l <- apply(fac[fac_col_exclude], 1, function(x) digest(x))
    facbin0a <- outer(fac.l, fac.l, "!=")
    facbin0 <- facbin0 & facbin0a
  }

  if(!is.null(select.cols)) {
    prune.vec <- which(obj$data[,select.cols])
    facbin0 <- facbin0[prune.vec,]
    Dx <- Dx[prune.vec,]  
  }
  
  if(!any(facbin0)) {
    message('facbin0 is all FALSE. Lift curve cannot be generated')
    return(NULL)
  }

  # Exclude the first column because it will always be True
  # (except if it has been zero'ed out based on the fac_col_exclude column,
  # in which case it is still ok to remove the column)
  liftcurve.helper(Dx, 
                   facbin0, 
                   fac, 
                   exclude_first_col = TRUE,
                   compute_closest = compute_closest)

}
         



loo.liftcurve.profile.data <- function(obj, fac1_col, fac2_col,
                                       pick_best=FALSE, type='eqwt', ...) {
  # Loop over each hairpin
  # Leave out the hairpin and compute CGS for gene
  # Compute CGS for all other genes
  # Compute distance between held out hp and all genes

  # Assume n genes and m hairpins
  # facbin is an m x n matrix, where 1 indicates the gene corresponding to the hairpin
  # Dx     is an m x n matrix, where (i,j) is the distance of the i'th hairpin from the th j'th gene

  #remove rows corresponding to singleton fac2_cols
  obj <- remove.singletons(obj, fac2_col)
  fac <- factors(obj)

  # check if rows are unique wrt fac1_col
  if(length(unique(fac[, fac1_col])) != length(fac[, fac1_col])) {
    dupvec <- duplicated(fac[, fac1_col])
    message(sprintf('Rows are not unique wrt %s. Pruning by arbitrarily removing all but one of the repetitions. %d repetitions removed.', fac1_col, sum(dupvec)))
    obj <- prune.rows(obj, !dupvec)
    fac <- factors(obj)    
  }
  

  Dx <- NULL
  facbin0 <- NULL
  
  # Precompute P0. One row at a time will be replaced in the loop
  P0 <- wtavg(obj, metacols = c(fac2_col), type=type)
  
  # check if rows are unique wrt fac_col2
  expect_true(setequal(P0$data[,fac2_col], unique(P0$data[,fac2_col])))
  
  # loop over each row
  #pb <- txtProgressBar(min = 1, max = NROW(fac), style = 3)

  d_fb <- foreach(i = seq(NROW(fac)), .combine=rbind) %dopar% {
    
    #setTxtProgressBar(pb, i)

    # Prune to only hps of the same gene as that of the i'th row hp
    # and exclude the i'th row hp
    P1 <- prune.rows(obj,
                     (fac[,fac2_col] == fac[i,fac2_col]) &
                     (seq(NROW(fac)) != i))

    if (!pick_best) {
      # compute gene signature
      P1 <- wtavg(P1, metacols = c(fac2_col), type = type)
    } else {
      v1 <- feats(obj)[i,]
      d <- unname(aaply(as.matrix(feats(P1)), 1,
                        function(x) sum((x-v1)^2),
                        .parallel=T))
      P1$data <- P1$data[which.min(d),
                         c(fac2_col, obj$feat_cols)]
    }

    # verify that only one row remains after computing gene signature
    expect_equal(NROW(P1$data), 1)
    
    # Replace row of P0 that corresponds to the gene of the i'th hp
    P2 <- P0
    r <- which(P2$data[,fac2_col]  == fac[i, fac2_col])
    # verify that exactly one row corresponds to the selected gene
    expect_equal(length(r), 1)
    # verify that headers are the same
    expect_true(all(names(P2$data) == names(P1$data)))
    # replace row
    P2$data[r,] <- P1$data

    
    # compute distance between all gene signatures and hp
    v1 <- feats(obj)[i,]
    d <- unname(aaply(as.matrix(feats(P2)), 1, function(x) sum((x-v1)^2),
                      .parallel=T))

    
    # create indicator vector indicating the gene of the hp
    fb <- as.character(P2$data[, fac2_col]) ==
        as.character(obj$data[i, fac2_col]) # don't use factors()
    # verify that exactly one indicator is True
    expect_true(sum(fb) == 1)

    
    # concatenate vars
    c(d, fb)
    
  }
  #close(pb)

  # Split d__fb into Dx and facbin0
  Dx <- d_fb[,1:(NCOL(d_fb)/2)]
  facbin0 <- d_fb[,(1+NCOL(d_fb)/2):NCOL(d_fb)]

  liftcurve.helper(Dx, facbin0, fac)
}


wtavg.profile.data <- function(obj, metacols, type = 'eqwt', ...) {
  
  # Usage: e.g.
  # Combine replicates by taking (weighted) averages
  # Pa <- wtavg(Pf, metacols = c('Gene','Treatment', 'TreatmentAbbrev',
  #                              'TargetSeq', 'Type',
  #                              'X7merSeed_11to17', 'X6merSeed_12to17'), type='eqwt')

  if (type == 'eqwt') {
    obj$data <- ddply(obj$data, metacols,
                      function(X) colMeans(X[,obj$feat_cols]))    
  } else if (type == 'corwt') {
    obj$data <- ddply(obj$data, metacols,
                      function(X) corwt1(X[,obj$feat_cols]))    
  }
  
  obj$factor_cols <- metacols
  obj
}


read.outlier.data.profile.data <- function(obj, ...) {
  version <- getd(obj$cfg, 'version', 0)
  
  fname <- path(obj, obj$cfg$qc_file)
  qv <- read.csv(fname)

  if (version < 1) {
    qv <- paste(qv$Plate, qv$Well, sep=":")
    ov <- paste(plate.orig(obj), well(obj), sep=":")
    not_outlier <- !(ov %in% qv)
  } else if (version < 2) {
# HACK - only QCFlag_isBlurry is used for throwing out outliers
    qv$not_outlier <-  qv$QCFlag_isBlurry == 0
    ov <- join(data.frame(Plate = plate.orig(obj), Well = well(obj)), qv)
    not_outlier <- ov$not_outlier
  }
  not_outlier  
}

score.cells.profile.data <- function(obj, rules_filename, ...) {
  s <- readLines(file(rules_filename))
  r = str_match(s, "(IF \\()(.*) (>) (.*), \\[(.*)\\], \\[(.*)\\]\\)")
  feat_name <- r[,3]
  thres <- as.numeric(r[,5])
  
  K <- length(strsplit(r[1,6], ",")[[1]])  #number of classes
  M <- length(s) # number of rules
  
  cond_pos <- cond_neg <- matrix(nrow=M, ncol=K)
  
  for(i in 1:M) {
    for(j in 1:K) {
      cond_pos[i,j] <- (as.numeric(strsplit(r[i,6], ",")[[1]])[j])
      cond_neg[i,j] <- (as.numeric(strsplit(r[i,7], ",")[[1]])[j])
    }
  }
  
  X <- feats(obj)  
  N <- NROW(X)
  scores <- matrix(0.0, nrow=N, ncol=K)
  for(i in 1:N) {
    for(j in 1:M){    
      if (X[i,feat_name[j]] > thres[j]) v <- cond_pos[j,] else v <- cond_neg[j,]
      for(l in 1:K)
        scores[i,l] <- scores[i,l] + v[l]
    }
  }
  
  classp <- matrix(nrow=N, ncol=1)
  
  for (i in 1:N)
    classp[i] <- which.max(scores[i,])
  
  classp
  
}

plot.plates.treatment.profile.data <- function(obj,
                                               numrows = 16,
                                               numcols = 24, 
                                               overlay.outlier = TRUE,
                                               fill.col = NULL,
                                               fill.as.factor = TRUE,
                                               label.col = NULL, 
                                               use.data.for.plotting = FALSE,
                                               ...) {
  
  #expect_true(obj$cfg$version > 0)
  expect_false(overlay.outlier)
  
  if (use.data.for.plotting) {
    data <- obj$data
  } else {
    data <- obj$data_plate_plot    
  }

  all.wells <- ldply(rev(LETTERS[seq( from = 1, to = 16 )]),
                     function(x) ldply(1:numcols,
                                       function(y) ldply(unique(data$Plate),
                                                         function(z) 
                                                           data.frame(
                                                             WellRow=x,
                                                             WellCol=y,
                                                             Plate=z))))

  data <- merge(all.wells, data, all = TRUE)
  
  if (!is.null(label.col))
    data$label.col <- data[,label.col]
  
  if (!is.null(fill.col)) {
    data$fill.col <- data[,fill.col]
    if (fill.as.factor) {
        data$fill.col <- as.factor(data$fill.col)
    }
  }
  
  p <- ggplot()
  
  if (!is.null(fill.col)) {
    p <- p + geom_tile(data=data, aes(as.integer(WellCol), 
                                      as.integer(WellRow),
                                      fill=fill.col))
  } else {
    p <- p + geom_tile(data=data, aes(as.integer(WellCol), 
                                      as.integer(WellRow)),
                       fill="white")
  }
  
  if (!is.null(label.col))
    p <- p + geom_text(data=na.omit(data), 
                       aes(as.integer(WellCol), 
                           as.integer(WellRow),
                           label=label.col), size=2)
  
  p <- (p 
        + facet_wrap( ~Plate) + coord_equal() + xlab("") + ylab("")
        + coord_fixed(ratio=1,
                      xlim = c(0.5,numcols+0.5),
                      ylim = c(0.5,numrows+0.5))
        + scale_x_discrete(breaks=1:numcols,
                           labels=1:numcols) 
        + scale_y_discrete(breaks=1:numrows,
                           labels=rev(LETTERS[seq( from = 1, to = numrows )])))


  p <- p + guides(fill=guide_legend(title=fill.col, nrow = 12, byrow = TRUE))
  p <- p + theme_bw()  
  p
}

rbind2.profile.data <- function(obj, P) {
  if (is.null(obj))
    return(P)
  if (is.null(P))
    return(obj)
  
# ensure that the column names are the sames
  expect_true(setequal(names(obj$data), names(P$data)))                       
# Remove attributes that may no longer be consistent across datasets
  obj$cfg <- NULL 
  obj$data_plate_plot <- NULL
  obj$data <- rbind(obj$data, P$data)
  obj
}

merge2.profile.data <- function(obj, P, discard.cfg=TRUE) {
  if (is.null(obj))
    return(P)
  if (is.null(P))
    return(obj)

  merge.cols <- intersect(obj$factor_cols, P$factor_cols)
  expect_true(length(merge.cols) > 0)

  obj$data <- merge(obj$data, P$data, by = merge.cols)
  
  # Remove attributes that may no longer be consistent across datasets
  if (discard.cfg)
      obj$cfg <- NULL 
  obj$data_plate_plot <- NULL
  obj
}

qualitymeas.profile.data <- function(obj, 
                                     metacols,
                                     exclude.cols=NULL,
                                     select.cols=NULL,
                                     cmpfunc='spearman',
                                     ctrl=NULL, collapse.ctrl=FALSE,
                                     collapse.trt=FALSE,
                                     subsample.ctrl=NULL,
                                     summarize.quantiles=TRUE,
                                     omit.na=TRUE,
                                     ...) {

  obj$factor_cols <- c(as.character(metacols), exclude.cols, select.cols)
  obj$data <- obj$data[, c(obj$feat_cols, obj$factor_cols)]

  if (!is.null(ctrl)) {
    type <- "str"

    obj_ctrl <- obj
    obj_ctrl$data <- merge(obj_ctrl$data, ctrl)
    X_ctrl <- feats(obj_ctrl)

    if (!is.na(subsample.ctrl)) {
        X_ctrl <- X_ctrl[sample(NROW(X_ctrl), subsample.ctrl),]
    }

    if (collapse.ctrl) {
        X_ctrl <- colMeans(X_ctrl)
    }
    
    if (cmpfunc == "euc") {
        h <- function(x, y) sqrt(sum((x-y)**2))
    } else if (cmpfunc %in% c("pearson", "kendall","spearman")) {
        h <- function(x, y)  1-cor(x, y, method=cmpfunc)
    }

    expect_true(collapse.ctrl | !collapse.trt)
    
    # E will just be ignored when computing signal strength
    if (collapse.trt) {
        g <- function(X, E) h(colMeans(X), X_ctrl)
    } else {
        if (collapse.ctrl) {
            g <- function(X, E) as.vector(aaply(as.matrix(X), 1,
                                             function(x) h(x,X_ctrl) ))
        } else {
            g <- function(X, E) as.vector(aaply(as.matrix(X), 1,
                                             function(x)
                                             aaply(as.matrix(X_ctrl),
                                                   function(y) h(x,y) )))
        }
    }

    
    
  } else {
    type <- "sim"
  
    cor.upper <- function(M, E) {
      cm <- cor(M, method=cmpfunc)
      M <- upper.tri(cm) 
      # Remove elements that have the same values in E
      if (NCOL(E) > 0) {
          if (NCOL(E) > 1)
            E <- apply(E, 1, function(Ei) paste(Ei, collapse="."))
          E <- outer(E, E, FUN = "!=")
          M <- M & E
      }
      cm[M]
    }
  
    g <- function(X, E) cor.upper(t(X), E)
    
  }

  if (summarize.quantiles) {
    f <- function(X, E) {
      v <- g(X, E)
      qn <- c(50,75)
      D <- data.frame(t(quantile(v, qn/100, names=F)))
      names(D) <- paste(type, "_", cmpfunc, "_q", qn, sep="")
      D
    }    
  } else {
    f <- function(X, E) data.frame(t(g(X, E)))      
  }

  if(!is.null(select.cols)) {
    message("Handling of select.cols in quality.meas gives a non-intuitive result. Read code for specifics.")
    expect_true(length(select.cols)==1)
    # Prune the entire data down to only those rows where select.cols is true
    obj$data$select.cols <- obj$data[, select.cols]
    obj$data <- droplevels(subset(obj$data, select.cols))
  }

  obj$data <- ddply(obj$data, 
                    metacols, 
                    function(D) f(D[,obj$feat_cols],
                                  D[,exclude.cols]), .parallel=F
                    )
  
  if(omit.na & any(is.na(obj$data))) {
    message("Omitting NA's!")
    obj$data <- na.omit(obj$data)
  }
  

  obj$feat_cols <- setdiff(names(obj$data), metacols)
  obj$factor_cols <- metacols
  obj
}


apply.well.position.profile.data <- function(obj, ...) {
  message('Stub!')  
  obj
    
}


medpolish.plate.profile.data <- function(obj, ...) {
  
  for (plate in unique(obj$data$Plate)) {
    X0 <- subset(obj$data, Plate==plate)
    X1 <- 
      foreach (feat = obj$feat_cols, .combine = 'cbind') %dopar% {
        X <- X0[, c("WellRow", "WellCol", feat)]
        Xc <- melt(medpolish(acast(X, WellRow ~ WellCol,  value.var=feat), 
                             na.rm=T, trace.iter=F)$residuals,
                   varnames=c("WellRow", "WellCol"),
                   value.name=feat)
        Xc$WellRow <- ordered(Xc$WellRow, levels=
                                rev(LETTERS[seq(from=1, 
                                                to=nlevels(X$WellRow))]))
        Xc <- na.omit(Xc)
        Xc <- merge(X[,c("WellRow", "WellCol")], Xc, sort=F)      
        expect_true(all(obj$data[obj$data$Plate==plate, "WellRow"] == 
                          Xc$WellRow))
        expect_true(all(obj$data[obj$data$Plate==plate, "WellCol"] == 
                          Xc$WellCol))
        Xc[, feat]
      }
    colnames(X1) <- obj$feat_cols
    obj$data[obj$data$Plate==plate, obj$feat_cols] <- X1
    
  }
  obj
  
}


dist1.profile.data <- function(obj, dist.func = "spearman.dist",
                                     ...) {

  X <- feats(obj)
  if (dist.func == "spearman.dist") {
    f <- function(X) (1 - cor(X, method="spearman"))
  }
  
  cm <- as.data.frame(f(t(X)))
  obj$data <- cbind(factors(obj), cm)
  obj$feat_cols <- names(cm)
  obj 
  
}



prune.to.identical.freq.for.labels.profile.data <- function(obj, 
                                                            label.col,
                                                            ...) {
  # Prune the profile.data so that the number of rows with the same value in
  # label.col is the identical across all unique values of label.col
  # This pruning is typically done before generate.null with the 
  # constrain.groups set to TRUE. This is because the generate.null function
  # (with the constrain.groups=T) requires that each group has the same number
  # of labels
  
  obj$data <- droplevels(obj$data)
  label.col <- as.character(label.col)
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))] 
  }
  
  tbl <- table(obj$data[, label.col])
  
  if(any(tbl!=Mode(tbl))) {
    obj$data <- droplevels(obj$data[obj$data[, label.col] %in% 
                                      names(tbl[tbl==Mode(tbl)]),])
    message(sprintf("Discarding rows with %s=%s", label.col, 
                    paste(names(tbl[tbl!=Mode(tbl)]), collapse=",")))
  }
  
  obj
}

