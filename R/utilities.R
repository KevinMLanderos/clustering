#!/usr/bin/R

##############################
# Clustering handy functions #
##############################

# This functions are designed to be help with object's handling

# Find 10X directory (barcodes, gene_names, counts), Seurat object, CSV file, TXT, RDS
get_source_data <- function(
  xpath,
  pj = "10XData",
  metadata = NULL,
  merge_counts = FALSE,
  verbose = FALSE
) {
  if(is.null(metadata)) metadata <- "no_file.csv"
  if(length(xpath) > 1) merge_counts = TRUE
  if(verbose) cat("======================= Expression matrix:\n")
  if(grepl("h5$", xpath) && !requireNamespace("hdf5r", quietly = TRUE))
    xpath <- gsub(".outs\\/.*", "", xpath)
  library_names = basename(gsub("outs\\/.*", "", xpath))
  if(all(dir.exists(xpath)))
    xpath <- sapply(xpath, dir.find_sc, pj = pj, verbose = verbose)
  if(verbose) cat(xpath, sep = "\n"); xpath_bk = xpath
  if(!isTRUE(merge_counts)){
    if(dir.exists(xpath)){
      if(verbose) cat(" - 10X directory:", xpath, " ")
      found_edata <- Seurat::Read10X(data.dir = xpath)
    }else if(file.exists(xpath)){
      if(verbose) cat(" - File...");
      found_edata <- readfile(xpath, verbose = verbose)
      otype <- casefold(class(found_edata)); if(verbose) cat("  *", otype);
      if(otype == "seurat"){
        if(verbose) cat(" v", as.character(found_edata@version), "\n", sep = "");
        if(verbose) cat("  Taking @meta.data and @raw.data (or [['RNA']]@counts)\n");
        found_mdata <- found_edata@meta.data
        tvar <- grepl("^2.|^1.", found_edata@version)
        found_edata <- if(tvar) found_edata@raw.data else found_edata[["RNA"]]@counts
      }else if(otype == "saver"){
        found_edata <- found_edata$estimate
      }else if(otype == "singlecellexperiment"){
        print(found_edata)
        found_mdata <- data.frame(SingleCellExperiment::colData(found_edata))
        tmp <- intersect(c("X", "counts"), names(found_edata@assays@data))
        if(length(tmp) == 0) tmp = 1
        cat("  Taking assay", tmp, "\n")
        found_edata <- found_edata@assays@data[[tmp[1]]]
      }
      if(verbose) cat("\n"); xpath <- dirname(xpath)
    }; if(verbose) cat("\n")
  }else{
    if(verbose){ cat("Input for merging:\n"); print(xpath) }
    lnames <- if(length(xpath) == 1 || any(grep("csv$", xpath))){
      fname = if(dir.exists(xpath)){
        list.files(xpath, pattern = "aggregation", full.names = TRUE)
      }else{ xpath }
      y <- data.table::rbindlist(lapply(fname, readfile, stringsAsFactors = FALSE))
      z <- apply(X = y, MARGIN = 2, FUN = function(x) sum(file.exists(x)) )
      y[, max(which(z == max(z)))[1]]
    }else{ xpath }
    lnames <- unname(sapply(unique(lnames), function(x){
      if(!dir.exists(x)) return(paste0(dirname(x), "/filtered_feature_bc_matrix"))
      dir.find_sc(x)
    }))
    names(lnames) <- basename(dirname(dirname(lnames)))
    fnames <- lnames #gsub("beegfs/ciro", "BioAdHoc/Groups/vd-vijay/cramirez", lnames)
    if(verbose) cat(paste0(fnames, collapse = "\n"), "\n")
    fnames <- fnames[dir.exists(fnames)]
    if(verbose) cat("Merging libraries:", length(fnames), "/", length(lnames), "\n")
    found_edata <- lapply(names(fnames), function(x){
      if(verbose) cat("."); y <- Seurat::Read10X(fnames[x])
      if(all(grepl("\\-1", colnames(y)))) colnames(y) <- gsub("\\-.*", "", colnames(y))
      colnames(y) <- sub(paste0(x, "_"), "", paste0(colnames(y), "-", which(names(lnames) == x)))
      y
    }); if(verbose) cat(" ")
    allfeatures <- unique(unlist(lapply(found_edata, rownames)))
    if(verbose) cat("All features:", show_commas(allfeatures), "\n")
    tvar <- rowSums(sapply(found_edata, function(x) allfeatures %in% rownames(x) ))
    allfeatures <- allfeatures[tvar == length(found_edata)]
    if(verbose) cat("Overlaping features:", show_commas(allfeatures), "\n")
    found_edata <- lapply(found_edata, function(x){
      x[allfeatures, ]
    })
    found_edata <- do.call("cbind", found_edata)
  }; gc()
  # if found_mdata exists and metadata file doesn't
  tvar <- if(is.character(metadata)) xpath_bk == metadata && !file.exists(metadata) else TRUE
  if(exists("found_mdata") && tvar){
    return(list(edata = found_edata, source = xpath, mdata = found_mdata))
  }

  if(length(found_edata) > 1 && is.null(dim(found_edata))){
    if(verbose) cat("Keeping only expression data\n")
    tvar <- grep("express", names(found_edata), ignore.case = TRUE)
    found_edata <- found_edata[[tvar]]
  }
  if(is.data.frame(found_edata))
    found_edata <- Matrix::Matrix(as.matrix(found_edata), sparse = TRUE)
  str(found_edata, max.level = 2); gc()
  tvar <- paste0(paste0(".*", library_names, "_"), collapse = "|")
  colnames(found_edata) <- gsub(tvar, "", colnames(found_edata))
  colnames(found_edata) <- gsub(".*filtered_feature_bc_matrix_", "", colnames(found_edata))

  if(verbose) cat("======================= Annotation/metadata\n")
  # This file must be a table with the following structure
  # Library, Condition_1, Condition_2, etc, or it can be the annotation/meta.data
  if(verbose && is.character(metadata)) cat("File:", metadata[1], "\n")
  grpsamples <- if(!is.null(dim(metadata))){
    if(verbose) cat("* Given!\n"); metadata
  }else if(file.exists(metadata[1])){
    if(verbose) cat("* Reading!\n")
    readfile(metadata[1], row.names = 1, stringsAsFactors = FALSE,
      check.names = FALSE, header = TRUE, verbose = verbose)
  }else{
    if(verbose) cat("* Making table of library names!\n")
    data.frame(library = library_names, stringsAsFactors = FALSE)
  };
  if(!is.null(grpsamples$cellname)) rownames(grpsamples) <- grpsamples$cellname
  # find which column contains barcodes
  tvar <- apply(grpsamples, 2, function(x) sum(x %in% colnames(found_edata)) )
  # then check if they're in the rownames
  tvar["in_rownames"] <- sum(rownames(grpsamples) %in% colnames(found_edata))
  if(tvar["in_rownames"]>2 && verbose) cat("Sample names in rows\n")
  if(verbose) str(grpsamples)
  if(any(tvar > 2)){
    # check if any of the columns in the table has the cell names (barcodes)
    if(verbose) cat("- provided from file\n"); # is it the annotation itself?!
    found_mdata <- grpsamples; rm(grpsamples)
    # if the barcodes are not in the rownames, then find the column with them
    tmp <- head(names(tvar[which(tvar == max(tvar))]), 1)
    if(!tvar["in_rownames"]) rownames(found_mdata) <- found_mdata[, tmp]
  }else if(exists("found_mdata")){
    if(verbose) cat("- given in expression input\n");
    print(str(found_mdata))
  }else if(verbose){
    cat("Matrix:\n"); print(head(colnames(found_edata)))
    cat("Table:\n"); print(head(rownames(grpsamples)))
  }

  if(!exists("found_mdata")){
    if(verbose) cat(" - building table...\n");
    found_mdata <- data.frame(
      row.names = colnames(found_edata), stringsAsFactors = FALSE,
      check.names = FALSE)
    # You do need to find aggregation_csv to know the order of the libraries
    aggname <- if(any(dir.exists(xpath))) "aggregation.*csv" else basename(xpath)[1]
    path_i <- if(any(dir.exists(xpath))) xpath[1] else dirname(xpath)[1]
    aggregated <- file.find(
      name = aggname, path = path_i, read = TRUE,
      verbose = verbose, stringsAsFactors = FALSE)
    if(is.null(aggregated)) aggregated <- data.frame(library = library_names)
    # getting cells GEM
    cell_gem <- sub(".*\\-", "", colnames(found_edata), perl = TRUE)
    # and vector named as the library name
    cell_lib <- rownames(aggregated); names(cell_lib) <- aggregated[, 1]
    if(verbose){
      cat("**\nAssigning names per library\n", show_commas(cell_lib), "\n")
      cat(show_commas(names(cell_lib)), "\n***\n")
    }

    grpsamples$its_rnames <- rownames(grpsamples)
    tvar <- sapply(grpsamples, function(x) sum(x %in% names(cell_lib)) )
    libname <- names(which(tvar == max(tvar))); if(verbose) str(libname)
    tvar <- grpsamples[, tail(libname, 1), drop = FALSE]
    tvar <- !sapply(tvar, function(x) any(duplicated(x)) )
    libname <- libname[tvar][1]
    if(is.na(libname)) libname <- colnames(grpsamples)[1]
    if(verbose) cat("Libraries column:", libname, "\n")
    # taking the library names as rownames for ordering
    rownames(grpsamples) <- grpsamples[, libname]

    if(!all(names(cell_lib) %in% rownames(grpsamples)) && verbose){
      tvar <- !names(cell_lib) %in% rownames(grpsamples)
      cat("Missing:", show_commas(names(cell_lib)[tvar]), "\n")
      tvar <- !rownames(grpsamples) %in% names(cell_lib)
      cat("Given:", show_commas(grpsamples[tvar, libname]), "\n")
    }
    # ordering accorgin to GEM from aggregation file
    grpsamples <- grpsamples[names(cell_lib), , drop = FALSE]
    grpsamples_selected <- grpsamples
    rownames(grpsamples) <- NULL;
    colnames(grpsamples) <- gsub("its_rnames", "Library", colnames(grpsamples))
    tvar <- !colnames(grpsamples_selected) %in% "its_rnames"
    grpsamples_selected <- grpsamples_selected[, tvar, drop = FALSE]
    # witching to GEM as names in the named vector
    cell_lib <- rownames(grpsamples_selected)
    names(cell_lib) <- rownames(grpsamples_selected) <- rownames(aggregated)
    group_names <- grpsamples_selected

    print(data.frame(GEM = rownames(aggregated), Library = cell_lib), row.names = F)
    if(verbose){ cat('\nAssigning group names\n'); print(group_names) }
    # don't add origlib if the library names are already in one of the columns
    tvar <- !any(sapply(group_names, function(x) all(cell_lib %in% x) ))
    if(tvar) found_mdata <- cbind(found_mdata, origlib = cell_lib[cell_gem])
    found_mdata <- cbind(found_mdata, group_names[cell_gem, , drop = FALSE])
  };

  return(list(edata = found_edata, source = xpath, mdata = found_mdata))
}

dir.find_sc <- function(x, pj = '10XData', verbose = FALSE){
  mypatterns <- c(
    'outs', 'filtered_feature_bc_matrix', 'filtered_gene_bc_matrices_mex',
    pj, 'hg19')
  dir.find(x, pattern = mypatterns, verbose = verbose)
}
# set directory with names
dir.find <- function(x, pattern = NULL, verbose = FALSE){
  pattern <- unique(unlist(strsplit(pattern, "/")))
  if(!grepl("\\/$", x)) x <- paste0(x, "/")
  if(verbose) cat("Setting folder with root:", x, "\n")
  if(verbose > 1) cat("Patterns:", pattern, sep = "\n")
  n <- length(pattern)
  for(i in 1:n){
    tmp <- list.files(x)
    tmp <- tmp[tmp %in% pattern]
    if(length(tmp) == 1){
      if(verbose) cat(rep(' ', i - 1), "|--", tmp, "\n")
      x <- paste0(x, tmp, '/')
    }
  }
  return(x)
}

ident_info_template = list(names = c("0" = "celltype1"), colours = NULL, order = NULL)
ident_order = function(ident){
  tmp = if(!is.null(ident$order)) ident$order else names(ident[[1]])
  if(is.null(tmp)) return(ident)
  tvar <- names(ident) != "order"
  ident[tvar] <- lapply(ident[tvar], function(x) x[as.character(ident$order)] )
  ident
}
ident_interactions = function(
  ident,
  cluster_colname = "cluster",
  interactions = NULL,
  interactions_sep = ": "
) {
  ident2add = ident[!names(ident) %in% c("order", "colours", "colors", "group_colours")]
  if(is.null(interactions)){
    interactions <- paste0(rep(cluster_colname, length(ident2add)), interactions_sep, names(ident2add))
  }
  if(is.null(names(interactions)))
    names(interactions) <- gsub("\\.{2,}", ".", make.names(interactions))
  return(interactions)
}

# group_colours = list(column = list(col = "#BEBEBE", cluster = 0:2))
ident_colours_infer <- function(
  mdata,
  ident,
  cluster_colname = "cluster",
  verbose = TRUE
) {
  ident2add = ident[which(!names(ident) %in% c("colours", "colors", "order"))]
  y <- lapply(X = setNames(nm = names(ident2add)), FUN = function(i){
    if(verbose) cat(i, "\n")
    tmp <- table(mdata[, c(cluster_colname, i)])
    if(all(colSums(tmp > 0) == 1)) return(NULL)
    tvar <- apply(X = tmp, MARGIN = 1, function(x) colnames(tmp)[which.max(x)] )
    lapply(X = unique(tvar), FUN = function(x){
      list(cluster = names(tvar[tvar == x]), col = ident$colo[[x]])
    })
  }); y <- y[sapply(y, length) > 0]; if(verbose) str(y); unlist(y, recursive = FALSE)
}

ident_colours <- function(
  ident,
  mdata = NULL,
  cluster_colname = "cluster",
  fun = "lighten",
  verbose = TRUE
) {
  ident_i = which(names(ident) %in% c("colours", "colors"))
  col_final = NULL
  if(verbose) cat("Defining identities' colours\n")
  tmp <- !all(levels(mdata[, cluster_colname]) %in% names(ident$colo))
  if(is.null(ident$group_colo) && tmp){
    if(verbose) cat("Calculate grouped\n")
    ident$group_colours = ident_colours_infer(
      mdata = mdata, ident = ident, verbose = verbose > 1)
    col_final <- c(col_final, ident$colo)
  }
  if(!is.null(ident$group_colours)){
    if(verbose) cat("Using list (complex assignments)\n")
    for(i in seq(length(ident$group_colours))){
      if(verbose) cat(" -", i)
      col2set = ident$group_colours[[i]]$col
      col_f = ident$group_colours[[i]]$cluster
      fun_ii <- ident$group_colours[[i]]$fun
      fun_ii <- if(is.null(fun_ii)) fun else fun_ii
      fun_i = switch(
        EXPR = fun_ii,
        lighten = colorspace::lighten,
        darken = colorspace::darken,
        saturation = shades::saturation)
      col_i = if(length(col2set) > 1){
        if(verbose) cat(" > 1"); colorRampPalette(col2set)(length(col_f))
      }else if(length(col_f) <= 6){
        if(verbose) cat(" >", fun_ii)
        fun_i(col2set, seq(0, 1, by = 0.2)[1:length(col_f)])
      }else{
        if(verbose) cat(" >", fun_ii, length(col_f))
        colorRampPalette(fun_i(col2set, 0.7))(length(col_f))
      }; if(verbose) cat("\n")
      tmp <- setNames(col_i, col_f); if(verbose > 1) str(tmp)
      col_final <- c(col_final, tmp[setdiff(names(tmp), names(col_final))])
    }
  }else{
    if(verbose) cat("New assignments\n")
    ident2add = ident[!names(ident) %in% c("order", "colours", "colors", "group_colours")]
    col_final = v2cols(unique(unlist(ident2add)), if(length(ident_i)) ident[[ident_i]])
  }
  if(!is.null(mdata)){
    if(verbose) cat("From metadata\n")
    tvar <- as.character(mdata[, cluster_colname])
    mdata$colours = v2cols(tvar, col_final)[tvar]
    cnames = c(cluster_colname, names(ident), names(ident_interactions(ident)))
    for(i in cnames[cnames %in% colnames(mdata) & !cnames %in% "colours"]){
      if(verbose) cat(" -", i, "\n")
      colours_list <- base::split(x = as.character(mdata[, "colours"]), f = mdata[, i])
      colours_list <- lapply(colours_list, unique)
      tmp <- sapply(colours_list, function(x) tail(x, 1) )
      col_final <- c(col_final, tmp[setdiff(names(tmp), names(col_final))])
    }
  }
  if(length(ident_i)){
    ident[[ident_i]] <- c(col_final, ident[[ident_i]][!names(ident[[ident_i]]) %in% names(col_final)])
  }else{ ident$colours = col_final }
  return(ident)
}

ident_set_names = function(
  object,
  ident,
  cluster_column,
  cluster_colname = "cluster",
  interactions = NULL,
  interactions_sep = ": ",
  verbose = TRUE
){
  if(verbose) cat("---- Identitites ----\n")
  mdata = if(casefold(class(object)) == "seurat") object@meta.data else object
  if(verbose) cat("Dimentions:", dim(mdata), "\n")
  if(is.null(names(ident))) names(ident) <- paste0("names", seq(length(ident)))

  mdata[, cluster_column] <- droplevels(factor(mdata[, cluster_column]))
  clust_names = levels(mdata[, cluster_column])
  if(is.null(ident$order)) ident$order = clust_names
  ident$order = ident$order[ident$order %in% clust_names]
  ident <- ident_order(ident = ident)
  ident2add = ident[!names(ident) %in% c("order", "colours", "colors", "group_colours")]

  for(i in seq(length(ident2add)))
    if(is.null(names(ident2add[[i]]))) names(ident2add[[i]]) <- clust_names
  interactions = ident_interactions(
    ident = ident,
    cluster_colname = cluster_colname,
    interactions = interactions,
    interactions_sep = interactions_sep
  )

  tmp = names(ident2add[[1]])
  if(!is.null(tmp)){ # $order is transmited to the first list
    if(verbose) cat("Creating main column:", cluster_colname, "\n")
    tvar <- levels(mdata[, cluster_column])
    if(verbose) cat(" Re-level:", tvar, "\n")
    if(verbose && any(!tvar %in% tmp))
      cat("WARNING missing:", tvar[!tvar %in% tmp], "\n")
    mdata = mdata[as.character(mdata[, cluster_column]) %in% tmp, ]
    mdata[, cluster_colname] = factor(mdata[, cluster_column], levels = tmp)
    if(verbose) cat("  Levels:", levels(mdata[, cluster_colname]), "\n")
  }
  for(i in names(ident2add)){
    if(verbose) cat("Creating column:", i, "\n")
    if(verbose) cat(" Levels:", unique(ident2add[[i]]), "\n")
    mdata[, i] = factor(ident[[i]][as.character(mdata[, cluster_column])],
      levels = unique(ident2add[[i]]))
  }
  for(i in names(interactions)){
    i_j = unlist(strsplit(interactions[[i]], interactions_sep))
    if(verbose) cat("Interaction:", i, "\n")
    if(verbose) cat(" Merging", i_j, "\n")
    if(all(i_j %in% colnames(mdata))){
      tvar <- interaction(mdata[, i_j], sep = interactions_sep, lex.order = TRUE)
      mdata[, i] <- factor(tvar, levels(tvar)[levels(tvar) %in% as.character(tvar)])
    }
  }
  if(verbose) cat("Dimentions:", dim(mdata), "\n")
  if(casefold(class(object)) == "seurat") object@meta.data = mdata else object = mdata
  if(verbose) cat("---- ----------- ----\n")
  return(object)
}
ident_tag = function(
  mdata,
  tag = "orig~",
  pattern = NA,
  exclude = NULL,
  verbose = FALSE
) {
  tvar <- if(!is.na(pattern)){
    grep(pattern, colnames(mdata), value = TRUE)
  }
  tvar <- filters_columns(
    mdata, tvar,
    types = c("character", "factor", "logical"),
    verbose = verbose,
    exclude = paste0(c(exclude, gsub("~", ".", tag)), collapse = "|")
  )
  if(length(tvar) > 0){
    tvar <- colnames(mdata) %in% tvar
    if(verbose) cat("Tagging names:", colnames(mdata)[tvar], "\n")
    if(grepl("~$", tag)){
      if(verbose) cat("Prefix\n")
      cols <- paste0(gsub("~.*", ".", tag), colnames(mdata)[tvar])
      colnames(mdata)[tvar] <- cols
    }
    if(grepl("^~", tag)){
      if(verbose) cat("Suffix\n")
      cols <- paste0(colnames(mdata)[tvar], gsub(".*~", ".", tag))
      colnames(mdata)[tvar] <- cols
    }
  }
  return(mdata)
}
ident_dimensions = function(
  x, prefix = "dim_"
){
  tvar <- unique(gsub("_1|_2", ".*", grep(prefix, x, value = TRUE)))
  redux <- lapply(tvar, function(y) grep(y, x, value = TRUE) )
  names(redux) <- gsub(paste0(".*", prefix, "|\\.\\*"), "", tvar)
  return(redux)
}

ident_centers = function(
  mdata, x = "x", y = "y", identity = 1, fun = median
) {
  mdata$Identity = mdata[, identity]
  mdata$x1234 <- mdata[, x]; mdata$y1234 <- mdata[, y]
  centers = x %>% group_by(Identity) %>%
    summarize(x = fun(x = x1234), y = fun(x = y1234))
  mdata[, !colnames(mdata) %in% c("x1234", "y1234")]
}

# Create a list of the filters used for QC
filters_qc_limits <- function(f){ # accepts a list(), if $file is given, takes that
  if(is.null(f$file)) f$file = "none"
  sections = c("low", "high", "use")
  ftab <- if(file.exists(f$file)){ # QC filters first
    y <- readfile(f, stringsAsFactors = FALSE, row.names = 1)
    if(any(sections %in% rownames(y))) y else t(y)
  }else{
    tvar <- names(f)[!names(f) %in% c("file", "subset", "nSamples_expressed")]
    sapply(f[names(f) %in% tvar], function(x) as.numeric(unlist(x)[1:3]) )
  }
  ftab = data.frame(ftab, row.names = sections, stringsAsFactors = FALSE)
  ftab[is.na(ftab)] <- -1 # if not specified, don't apply it
  used <- unlist(ftab["use", ]); used[is.na(used)] <- ifelse(all(is.na(used)), 1, 0)
  named_boundary <- function(x, y, r){ x <- x[y]; x[is.na(x)] <- r; names(x) <- y; x }
  list(
    low = named_boundary(x = unlist(ftab["low", ]), y = colnames(ftab), -Inf),
    high = named_boundary(x = unlist(ftab["high", ]), y = colnames(ftab), Inf),
    apply = used > 0
  )
}

# add gene set percentages
add_pcts <- function(
  mdata,
  edata,
  feature_pattern = list(
    pctMitochondrial = c('^mt-', '^m-', '^hg19_mt-', '^mm10_mt-'),
    pctRibosomal = paste0(c('^rps', '^rpl'), "[0-9]"),
    pctHSprefix = c("^hsp[0-9]"),
    pctHeatShock = paste0("^hsp", letters[1:5])
  ),
  verbose = FALSE
){
  if(verbose) cat("Total samples:", nrow(mdata), "\n")
  edata <- edata[, rownames(mdata)]
  column_total <- Matrix::colSums(edata)
  void <- sapply(names(feature_pattern), function(pname){
    grep(
      pattern = paste0(feature_pattern[[pname]], collapse = "|"),
      x = rownames(edata), ignore.case = TRUE, value = TRUE
    )
  })
  void <- void[sapply(void, length) > 0]
  if(length(void) > 0){
    if(verbose) cat("Retrieving sets:\n")
    if(verbose) str(void)
    pcts <- data.frame(sapply(void, function(feature_patterned){
      y <- round(Matrix::colSums(x = edata[feature_patterned, , drop = FALSE]) / column_total * 100, 2)
      y[is.na(y)] <- max(y, na.rm = TRUE)
      y
    }))
    mdata <- joindf(pcts, mdata)
  }else{ if(verbose) cat("No pattern found\n") }
  mdata
}

get_elbow <- function(
  x,
  y = NULL,
  threshold = .90,
  decide = FALSE,
  verbose = FALSE
) {
  if(is.null(y)){ y <- x; x <- 1:(length(y)) }
  mindexes <- sapply(threshold, function(z) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    thh <- abs(quantile(d2, probs = z))
    indices <- which(abs(d2) >= thh)
    if(verbose){
      cat('SDEV thh:', thh, '\n')
      cat("D'':", commas(round(d2[indices], 3)), '\n')
      cat('SDEVs:', commas(round(y[indices], 3)), '\n')
      cat('Chosen:', x[max(indices)], '->', y[max(indices)], '\n')
    }
    max(indices)
  })
  if(length(threshold) == 1)
    return(list(elbows = mindexes, ptable = data.frame(y = max(y), x = mindexes + 0.25, N = 1), elbow = mindexes))
  elfreqs <- table(mindexes)
  if(isTRUE(decide)){
    bestfit <- round(mean(as.numeric(names(elfreqs))))
  }else{ bestfit <- mean(mindexes) }
  qposition <- data.frame(y = max(y), N = elfreqs[as.character(mindexes)])
  colnames(qposition) <- c("y", "x", "N"); qposition$x <- as.numeric(as.character(qposition$x)) + 0.25
  list(elbows = mindexes, ptable = qposition, elbow = bestfit)
}

create_prefix <- function(configuration){ # based on norm and variable_features from the config
  new_prefix <- if(!grepl(configuration$norm, "sctransform")){
    paste0(
      "_mean", casefold(configuration$variable_features$mean.cutoff[1], upper = TRUE),
      "_pct", casefold(configuration$variable_features$percen, upper = TRUE)
    )
  }else{ "_sctransform" }
  tvar <- (!configuration$variable_features$method %in% c("vst", "mvp", "mean.var.plot") ||
    configuration$variable_features$mean.cutoff[2] != 8 ||
    any(configuration$variable_features$dispersion.cutoff != c(1, Inf))) && !file.exists(configuration$variable_features$file)
  new_prefix <- if(tvar){
    tvar <- sapply(configuration$variable_features, function(x) file.exists(as.character(x[[1]])) )
    tvar <- if(!any(tvar)) !names(configuration$variable_features) %in% "file" else rep(TRUE, length(configuration$variable_features))
    hash <- digest::digest(paste0(unlist(configuration$variable_features[tvar]), collapse = ""), "md5", serialize = FALSE)
    paste0(new_prefix, "_", hash)
  }else if(file.exists(configuration$variable_features$file)){
    paste0(new_prefix, "_", gsub("\\..*", "", basename(configuration$variable_features$file)))
  }else{
    new_prefix
  }
  new_prefix
}

complex_object_fetch <- function(
  object_or_dir,
  path = "./",
  id = ".rds", # it alswo tells you the type of data
  verbose = FALSE
){
  if(class(object_or_dir)[[1]] == "character"){
    tvar <- if(file.exists(object_or_dir) && !dir.exists(object_or_dir)){
      object_or_dir # id = gsub("([^.]+)\\.[[:alnum:]]+$", "\\1", object_or_dir);
    }else{
      path <- dirname(normalizePath(object_or_dir))
      suppressWarnings(system(paste0("ls ", path, "/.object_stem*"), intern = TRUE))
    }
    # tvar <- tvar[sapply(tvar, function(x) grepl(gsub(".*ct_stem", "", x), id) )]
    if(length(tvar) > 1) stop("More than one object templates match ", id)
    if(verbose) cat("Getting object:", tvar, "\n")
    object_or_dir <- readRDS(tvar)
  }
  if(verbose) cat("Building object:", class(object_or_dir)[[1]], "\n")
  if(verbose) cat("Path:", path, "\n")
  myfiles <- list.files(path = path, pattern = id, all.files = TRUE, full.names = TRUE)
  myfiles <- gsub("\\/{2,}", "/", myfiles)
  myfiles <- myfiles[grep("\\.object", myfiles)]
  names(myfiles) <- gsub(paste0(".object_stem|.object_|.object|", id, "|\\.rds"), "", basename(myfiles))
  myfiles <- myfiles[names(myfiles) != ""]
  for(component in names(myfiles)){
    if(!file.exists(myfiles[component])) next
    if(verbose) cat("Adding", myfiles[component], "\n")
    component_content <- readRDS(myfiles[component]) # this doesn't work sometimes (can't find why)
    eval(expr = parse(text = paste0("object_or_dir@", component, " <- component_content")))
  }
  return(object_or_dir)
}

# Get top markers
get_top_n <- function(
  x,
  cname = 'cluster',
  dpct = 'Dpct',
  orderby = 'p_val_adj',
  filter_pattern = "^rps|^rpl|^mt-",
  filter_neg = TRUE,
  n = Inf,
  verbose = FALSE
) {
  if(verbose) cat("Getting top", n, "features per", cname, "\n")
  groups <- names(table(x[, cname]))
  x$gene <- as.character(x$gene)
  if(!dpct %in% colnames(x)) x[, dpct] <- 100 * (x$pct.1 - x$pct.2)
  y <- as.data.frame(data.table::rbindlist(lapply(groups, function(z){
    if(verbose) cat('.')
    dat <- x[which(as.character(x[, cname]) == z), ]
    tvar <- isTRUE(filter_neg) && any(dat[, dpct] > 0)
    if(tvar) dat <- dat[dat[, dpct] > 0, ]
    tvar <- !grepl(filter_pattern, dat$gene, ignore.case = TRUE)
    dat <- dat[tvar, ] # Removing RPS and RPL features
    if(sum(!tvar)) dat$ribomito_genes <- sum(!tvar)
    data.table::setorderv(dat, orderby)
    # head(dat, n)
    dat[dat$gene %in% head(unique(dat$gene), n), ]
  }), fill = TRUE)); if(verbose) cat('\n')
  return(y)
}

ClassifyScoring <- function(
  object,
  features,
  set.ident = FALSE,
  name = 'Set',
  verbose = FALSE,
  ...
) {
  cc.columns <- grep(pattern = name, x = colnames(x = object[[]]), value = TRUE)
  if(is.null(names(features))) names(features) <- paste0("S", 1:length(features))
  if(name %in% cc.columns){ warning(name, " pre-computed"); return(object) }
  if(verbose) str(features)

  # Renaming columns colliding with previous signatures; first the 'name'
  cc.columns <- make.names(c(cc.columns, name), unique = TRUE); name <- cc.columns[length(cc.columns)]
  cc.columns <- make.names(c(colnames(x = object[[]]), names(features)), unique = TRUE); # now the 'classes'
  names(features) <- cc.columns[tail(1:length(cc.columns), length(names(features)))]

  classes <- names(features)
  if(verbose) cat("Calculating scores\n")
  ctrl_n <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  if(ctrl_n > 1000) ctrl_n <- 100
  if(verbose) cat("Background size:", ctrl_n, "\n")
  object <- AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl_n,
    ...
  )
  cc.columns <- grep(pattern = name, x = colnames(x = object[[]]), value = TRUE)
  cc.scores <- object[[cc.columns]]
  object@meta.data <- object@meta.data[, !colnames(object@meta.data) %in% cc.columns]
  # rm(object.cc); gc(verbose = FALSE) # 'object' was duplicated from AddModuleScore
  if(verbose) cat("Classification based on score.\nNone: all < 0; Undecided: max score > 1.\n")
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores) {
      if (all(scores < 0)) {
        return("None")
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(classes[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', classes, name)
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c(classes, name)]
  if(verbose) cat("Adding to meta data:", paste0(c(classes, name), collapse = ", "), "\n")
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    if(verbose) cat("Setting classes as identities\n")
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Class'
  }
  return(object)
}

# centered log-ratio (CLR) normalization - for each feature!
clr_function <- function(x) {
  apply(X = x, MARGIN = ifelse(nrow(x) < ncol(x), 1, 2), FUN = function(x){
    return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
  })
}

## Statistics and extra information
markers_summary <- function(
  marktab,
  annot,
  cluster_column = NULL,
  covar = NULL,
  datavis = NULL,
  datatype = "CTS",
  verbose = FALSE
){
  marktab <- remove.factors(marktab)
  if(casefold(class(annot)) == "seurat"){
    if(verbose) cat("Taking data from seurat object\n")
    datavis <- GetAssayData(annot)
    annot <- annot@meta.data
    datatype <- "SeuratNormalized"
  }
  annot <- remove.factors(annot)
  datavis <- as.matrix(datavis[rownames(datavis) %in% as.character(marktab$gene), ])
  if(is.null(cluster_column)) cluster_column <- colnames(annot)[grep("res\\.", colnames(annot))[1]]
  if(verbose) cat("Cluster column:", cluster_column, "\n")

  if(verbose) cat("Number of clusters per gene\n")
  ncluster <- t(sapply(unique(marktab$gene), function(x){
    y <- marktab[marktab$gene == x, "cluster"]
    if(length(y) < 1) return(NA)
    c(nCluster = as.character(length(y)), sCluster = paste0(y, collapse = "&"))
  }))
  marktab <- cbind(marktab, ncluster[marktab$gene, ])

  stattab <- stats_summary_table(
    mat = datavis,
    groups = make_list(x = annot, colname = cluster_column, grouping = TRUE),
    rnames = rownames(datavis),
    datatype = datatype,
    verbose = verbose
  )
  if(all(covar %in% colnames(annot)) && !is.null(covar)){
    if(verbose) cat("Using", show_commas(covar), "as covariate\n")
    annot$covar <- do.call('paste', c(annot[c(cluster_column, covar)], sep = "_"))
    stattab <- cbind(stattab, stats_summary_table(
      mat = datavis,
      groups = make_list(x = annot, colname = "covar", grouping = TRUE),
      rnames = rownames(datavis),
      datatype = datatype,
      verbose = verbose
    ))
    grps <- unique(annot$covar)
  }
  marktab$resolution <- cluster_column
  marktab$nCluster <- length(unique(marktab$cluster))
  marktab$cluster <- as.character(marktab$cluster)
  marktab$n_cells <- as.numeric(table(annot[, cluster_column])[as.character(marktab$cluster)])
  marktab$cluster_n <- paste0(
    marktab$cluster, "; N=", marktab$n_cells, "; ",
    paste0(round(marktab$n_cells / nrow(annot) * 100, 2), "%")
  )
  tvar <- gtools::mixedsort(unique(as.character(marktab$cluster_n)))
  marktab$cluster_n <- factor(marktab$cluster_n, levels = tvar)
  marktab$cPct <- round(marktab$n_cells / nrow(annot) * 100, 4)
  marktab$n_markers <- as.numeric(table(marktab$cluster)[marktab$cluster])
  if(verbose) cat("Calculating distances within clusters\n")
  kha <- sapply(unique(marktab$cluster), function(i){ # mean distance per cluster
    thesecells <- filters_subset_df(c(cluster_column, i), annot, verbose = FALSE)
    thesecells <- sample(thesecells, min(c(length(thesecells), 500)))
    mydata <- t(datavis[marktab[marktab$cluster == i, 'gene'], thesecells, drop = FALSE])
    tvar <- mean(dist(mydata), na.rm = TRUE)
  });
  marktab$dist <- kha[marktab$cluster]
  if(verbose) cat("Mean differences per gene\n")
  marktab$gene_mean <- sapply(1:nrow(marktab), function(i){ # mean per gene per cluster
    stattab[as.character(marktab[i, 'gene']), paste0(marktab[i, 'cluster'], '_mean', datatype)]
  }) # mean of the rest of clusters per gene per cluster
  marktab$bmean <- unlist(sapply(unique(marktab$cluster), function(i){
    thesecells <- filters_subset_df(c(cluster_column, paste0('-', i)), annot, verbose = F)
    rowMeans(datavis[marktab[marktab$cluster == i, 'gene'], thesecells, drop = FALSE], na.rm = TRUE)
  }))
  marktab$Dmean <- (marktab$gene_mean - marktab$bmean)
  marktab$Dpct <- 100 * (marktab$pct.1 - marktab$pct.2)
  tvar <- table(marktab$cluster, marktab$Dpct < 0)
  marktab$nDpct_neg <- if(ncol(tvar) == 2) tvar[as.character(marktab$cluster), "TRUE"] else 0
  marktab$Significance <- (-log10(marktab$p_val_adj))
  tvar <- max(marktab$Significance[is.finite(marktab$Significance)])
  marktab$Significance[is.infinite(marktab$Significance)] <- tvar
  marktab <- cbind(marktab, stattab[marktab$gene, ])
  #print(head(stattab[, head(1:ncol(stattab))]))

  stattab$gene <- rownames(stattab)
  marktab$gene_name <- paste0("'", marktab$gene)
  stattab$gene_name <- paste0("'", rownames(stattab))
  marktab <- data.table::rbindlist(list(marktab, stattab[!rownames(stattab) %in% marktab$gene, ]), fill = TRUE)
  marktab <- data.frame(marktab, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(marktab) <- make.names(marktab$gene, unique = TRUE)
  marktab
}

dgea_seurat = function(
  object,
  results_prefix = "markers",
  config_markers = list(test = "MAST", avg_logFC = 0.25, p_val_adj = 0.05),
  cluster_column = NULL,
  return_table = FALSE,
  verbose = TRUE
) {
  if(verbose) cat("---- Markers ----\n")
  Idents(object) <- cluster_column
  running_fag <- paste0(results_prefix, "_running")
  if(file.exists(running_fag)) next
  results_f = paste0(dirname(results_prefix), "/.", basename(results_prefix), ".rds")
  meansf = paste0(results_prefix, "_fc", config_markers$avg_logFC, "_padj",
    config_markers$p_val_adj, "_summary_stats.csv")
  if(file.exists(meansf)){
    if(verbose) cat(' - From final file\n')
    cmarkers = read.csv(meansf, row.names = 1, check.names = FALSE,
      stringsAsFactors = FALSE)
  }else{
    if(file.exists(results_f)){
      if(verbose) cat(' - From file\n')
      cmarkers <- readRDS(results_f)
      cmarkers$gene <- sub("'", "", cmarkers$gene_name)
    }else if(!file.exists(paste0(".", results_prefix, "_running"))){
      if(verbose) cat("- Calling FindAllMarkers\n")
      writeLines(text = results_prefix, con = running_fag)
      cmarkers <- FindAllMarkers(
        object = object,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0,
        min.diff.pct = 0.05,
        test.use = config_markers$test,
        return.thresh = 0.2,
        verbose = verbose
      ); timestamp()
      file.remove(running_fag); cmarkers$gene_name <- paste0("'", cmarkers$gene)
      saveRDS(object = cmarkers, file = results_f)
    }; tvar <- cmarkers$avg_logFC >= config_markers$avg_logFC
    tmp <- cmarkers$p_val_adj <= config_markers$p_val_adj
    cmarkers <- cmarkers[which(tvar & tmp), ]
    if(!file.exists(meansf) && !is.null(cluster_column)){
      cmarkers <- markers_summary(
        marktab = cmarkers,
        annot = object@meta.data,
        datavis = GetAssayData(object),
        cluster_column = cluster_column,
        datatype = "SeuratNormalized",
        verbose = verbose
      ); write.csv(cmarkers, file = meansf)
    }
  }; if(isTRUE(return_table)) return(cmarkers)
  return(invisible(x = NULL))
}

## utilities ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clust_similarity <- function(df){
  df[, 1] <- as.character(df[, 1])
  df[, 2] <- factor(df[, 2])
  df[, 3] <- factor(df[, 3])
  list(
    JI = clusteval::cluster_similarity(labels1 = as.numeric(df[, 2]), labels2 = as.numeric(df[, 3])),
    ARI = mclust::adjustedRandIndex(x = df[, 2], y = df[, 3]),
    NMI = NMI::NMI(X = df[, 1:2], Y = df[, c(1, 3)])$value
  )
}

make_cluster_names <- function(vec){
  vec <- as.character(vec)
  level_names <- names(sort(table(vec), decreasing = TRUE))
  levels <- 0:(length(level_names) - 1)
  names(levels) <- level_names
  factor(x = levels[vec], levels = unname(levels))
}

add_reduction <- function(
  object,
  reduction,
  reduction.key = "Dim_",
  reduction.name = 'reduction'
){
  object@reductions[[reduction.name]] <- object@reductions[[1]]
  scells <- rownames(object@reductions[[1]]@cell.embeddings)
  object@reductions[[reduction.name]]@cell.embeddings <- as.matrix(reduction[scells, ])
  colnames(object@reductions[[reduction.name]]@cell.embeddings) <- paste0(reduction.key, 1:ncol(reduction))
  object@reductions[[reduction.name]]@key <- reduction.key
  return(object)
}

show_commas <- function(x, hn = 3){
  tvar <- length(x)
  if(tvar > 1){
    tmp <- paste0(head(x[-tvar], hn), collapse = ', ')
    connect <- ifelse(tvar-1 > hn, ' ... and ', ' and ')
    tmp <- paste0(tmp, connect, x[tvar])
    if(tvar-1 > hn) tmp <- paste0(tmp, ' (', tvar, ')')
  }else{
    tmp <- x
  }
  return(tmp)
}
# make a list of rownames out a column's categories or
# return categories for each group in another column
make_list <- function(x, colname = colnames(x)[1], col_objects = NULL, grouping = FALSE){
  x <- remove.factors(x)
  tvar <- lapply(unique(x[, colname]), function(y){
    rnames <- filters_subset_df(list(c(colname, y)), x, v = F)
    if(is.null(col_objects)) return(rnames)
    x[rnames, col_objects]
  })
  names(tvar) <- unique(x[, colname])
  if(isTRUE(grouping)){
    tmp <- rep(names(tvar), sapply(tvar, length))
    names(tmp) <- unlist(tvar); tvar <- tmp
  }
  return(tvar)
}

ident_set_object <- function(x, verbose = TRUE){
  for (i in x) {
    if(exists("sc_tcells")){
      command <- paste0(i, " = ident_set_names(object = ", i, ", ident = ",
        i, "_ident, cluster_column = ", i, "_clust)")
      if(verbose) cat(command, "\n")
      base::eval(base::parse(text = command), envir = .GlobalEnv)
      command <- paste0(i, "_ident <- ident_colours(", i, "_ident, mdata = ",
        i, "@meta.data)")
      if(verbose) cat(command, "\n")
      base::eval(base::parse(text = command), envir = .GlobalEnv)
    }
  }
}
