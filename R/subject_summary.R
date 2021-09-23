#!/usr/bin/R

######################
# Donor-wise summary #
######################

# This code creates tables fo with metadata and proportions per donor

source('/home/ciro/scripts/handy_functions/devel/file_reading.R')
# readfile
source('/home/ciro/scripts/handy_functions/devel/filters.R')
# filters_subset_df, factormix, remove.factors
source('/home/ciro/scripts/handy_functions/devel/utilities.R')
# make_list, factormix, remove.factors

subject_summary <- function(
  metadata_list,
  subject_id,
  append_cols = NULL,# list(tables = length(metadata_list), names = "colnames")
  append_num = "RNA_snn_res", # sames as append_cols
  rename_cells = NULL, # c("col1", "col2", ..., "pattern", "replacement")
  filters = NULL,
  filters_seq = list(), # same as filters
  verbose = TRUE
){
  if(verbose) cat("---- Subject summary analysis ----\n")
  mdata_list = if(any(sapply(metadata_list, is.character))){
    if(verbose) str(metadata_list, max.level = 1)
    lapply(X = metadata_list, FUN = function(x){
      if(is.character(x)) readfile(x) else x
    })
  }else{ metadata_list }
  if(is.null(names(mdata_list))) names(mdata_list) <- paste0("Step", 1:length(mdata_list))

  if(!is.null(rename_cells)){ # THIS IS NOT REALLY NEEDED; ONLY SUBJECT NAMES IS RELEVANT
    if(verbose){ cat("** Redefining cell names\n"); str(rename_cells) }
    if(!is.list(rename_cells)) rename_cells <- list(rename_cells)
    if(length(rename_cells) == 1) rename_cells <- rep(rename_cells, length(mdata_list))
    if(is.null(names(rename_cells))) names(rename_cells) <- names(mdata_list)
    mdata_list <- lapply(
      X = setNames(nm = names(mdata_list)),
      FUN = function(x) rename_rows(mdata = mdata_list[[x]], rules = rename_cells[[x]])
    ); names(mdata_list) <- names(rename_cells)
    str(mdata_list[[1]][, (ncol(mdata_list[[1]])-3):ncol(mdata_list[[1]])])
    str(mdata_list[[2]][, (ncol(mdata_list[[2]])-3):ncol(mdata_list[[2]])])
  }

  if(verbose) cat("** Defining subject_id", paste0(subject_id, collapse = ","), "\n")
  subject_id_name <- paste0("subject_name", length(subject_id))
  mdata_list <- lapply(
    X = mdata_list,
    FUN = function(mdata){
      mdata$uniqname123 <- do.call(paste, c(mdata[, subject_id, drop = FALSE], sep = "-"))
      colnames(mdata) <- gsub("^uniqname123$", subject_id_name, colnames(mdata)); mdata
  })
  if(verbose && !is.null(filters)) cat("** Filtering\n")
  mdata_list <- filters_fun(mlist = mdata_list, rules = filters, verbose = verbose)

  if(verbose) cat("** Summarising", subject_id_name, "\n")
  summary_df <- summary_df_fun(mlist = mdata_list, subject_ids = subject_id_name, sname = ".N")

  if(length(filters_seq) > 0){
    if(verbose){ cat("** Sequencial filter\n"); str(filters_seq) }
    if(is.null(names(filters_seq))) names(filters_seq) <- paste0("filt", 1:length(filters_seq))
    for(i in names(filters_seq)){
      if(verbose) cat(i, "\n")
      mdata_list <- filters_fun(mlist = mdata_list, rules = filters_seq[[i]], verbose = verbose)
      summary_dfi <- summary_df_fun(mlist = mdata_list, subject_ids = subject_id_name, sname = i)
      summary_df <- cbind(summary_df, summary_dfi)
    }
  }

  if(length(append_cols) > 0){
    if(verbose) cat("** Including categories\n")
    summary_df <- add_columns(
      append_x = append_cols,
      mlist = mdata_list,
      type = "cols",
      subject_ids = subject_id_name,
      summary_df = summary_df,
      verbose = verbose
    )
  }

  if(length(append_num) > 0){
    if(verbose) cat("** Including numbers\n")
    summary_df <- add_columns(
      append_x = append_num,
      mlist = mdata_list,
      type = "num",
      subject_ids = subject_id_name,
      summary_df = summary_df,
      verbose = verbose
    )
  }

  tmp <- lapply(names(mdata_list), function(x) grep(paste0("^", x), colnames(summary_df)) )
  summary_df <- summary_df[, unique(unlist(tmp))]
  if(verbose > 1) str(summary_df)
  cats <- colnames(summary_df)[sapply(summary_df, class) %in% "character"]
  nums <- colnames(summary_df)[!sapply(summary_df, class) %in% "character"]
  summary_df <- remove.factors(cbind(Name = rownames(summary_df), summary_df))
  colnames(summary_df)[1] <- subject_id_name
  if(verbose) cat("** Summarising", paste0(cats, collapse = ","), "\n")
  summary_cats <- lapply(cats, function(caty){
    y <- plyr::ddply(
      .data = summary_df[summary_df[, caty] != "NA", ],
      .variable = caty, function(slice){
        sapply(nums, function(x) sum(slice[[x]], na.rm = TRUE) )
      })
    y
  })

  summary_df_final <- data.table::rbindlist(c(list(summary_df), summary_cats), fill = TRUE)
  summary_df_final[is.na(summary_df_final)] <- ""
  if(verbose) cat("---- ------- ------- -------- ----\n")

  return(summary_df_final)
}

# it also creates column "cellnames"
rename_rows <- function(
  mdata,
  rules = c("\\-[0-9]{1,}-", "-") # c("col1", "col2", ..., "pattern", "replacement")
) {
  cols <- rules[rules %in% colnames(mdata)]
  sufix = if(length(cols) > 0){
    if(any(grepl("^rev$", rules))) cols <- rev(cols)
    do.call(paste, c(mdata[, cols, drop = FALSE], sep = "-"))
  }else rep("", nrow(mdata))
  newnames = if(!any(grepl("^cols_only$", rules))){
    names = rownames(mdata)
    if(any(grepl("^rev$", rules))){ names <- sufix; sufix <- rownames(mdata) }
    paste0(names, "-", sufix)
  }else{ sufix }; rules_sub <- rules[!rules %in% colnames(mdata)]
  if(length(rules_sub) > 1){
    newnames = gsub(rules_sub[1], rules_sub[2], newnames)
  }
  mdata$cellnames_bk <- rownames(mdata)
  mdata$cellnames <- rownames(mdata) <- newnames
  mdata
}

summary_df_fun = function(
  mlist,
  subject_ids = "subject_name1",
  sname = ""
) {
  sumdf <- lapply(
    X = mlist,
    FUN = function(mdata){
      y <- as.data.frame.matrix(t(table(mdata[, subject_ids], useNA = "always")))
      return(y[, colSums(y) > 1])
  })
  sumdf <- t(data.table::rbindlist(l = sumdf, fill = TRUE))
  colnames(sumdf) <- paste0(names(mlist), sname)
  sumdf
}

filters_fun = function(
  mlist,
  rules = NULL,
  verbose = TRUE
) {
  lapply(
    X = mlist,
    FUN = function(mdata){
      filters_complex(mdata = mdata, filters = rules, verbose = verbose)[["annotation"]]
      # mdata[filters_subset_df(x = rules, df = mdata, verbose = verbose), ]
  })
}

add_columns = function(
  append_x,
  mlist,
  type = "cols",
  subject_ids = "subject_name1",
  summary_df = NULL,
  verbose = TRUE
) {
  if(!is.list(append_x)) append_x <- list(tables = length(mlist), names = append_x)
  if(is.null(append_x$tables)) append_x$tables <- length(mlist)
  if(verbose) str(append_x)
  summary_df_append <- lapply(
    X = mlist[append_x$tables],
    FUN = function(mdata){
      cols <- append_x$names[append_x$names %in% colnames(mdata)]
      if(length(cols) == 0){
        cols <- grep(pattern = append_x$names[1], x = colnames(mdata), value = TRUE)
      }
      if(length(cols) == 0) return(NULL)
      apply_fun = if(type == "cols") sapply else lapply
      y <- apply_fun(X = cols, FUN = function(x){
        if(type == "cols"){
          sapply(make_list(x = mdata, colname = subject_ids, col_objects = x),
            function(x) paste0(sort(unique(x)), collapse = "|") )
        }else{
          mdata$tmp <- if(is.factor(mdata[, x])) mdata[, x] else factormix(mdata[, x])
          as.data.frame.matrix(table(mdata[, c(subject_ids, 'tmp')]), check.names = FALSE)
        }
      })
      data.frame(y, stringsAsFactors = FALSE, check.names = FALSE)
  })
  summary_df_append <- summary_df_append[sapply(summary_df_append, length) > 0]
  summary_df_append <- if(!is.null(summary_df)){
    if(length(summary_df_append) > 0){
      for(i in 1:length(summary_df_append)){
        i_name <- names(summary_df_append[i])
        i_table <- summary_df_append[[i]]
        colnames(i_table) <- paste0(i_name, ".", colnames(i_table))
        summary_df <- cbind(summary_df, i_table[rownames(summary_df), , drop = FALSE])
      }
      summary_df
    }else{ summary_df }
  }else{ summary_df_append }
  summary_df_append
}
