#!/usr/bin/R

################
# Cluster tree #
################

# This script will create plots from clustering analyses using 'clustree'
# package

# source("/home/ciro/scripts/clustering/R/clustree.R")
# generate_tree(
#   dpath = "/home/ciro/large/simon/results/clustering/SiEs09_tfr_mouse_test2"
# )
source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
source("/home/ciro/scripts/handy_functions/devel/utilities.R") # column_collapse, joindf
source("/home/ciro/scripts/clustering/R/utilities.R") # ident_dimensions

suppressPackageStartupMessages({
  library(clustree)
})

clustree_plot = generate_tree = function(
  dpath,
  pattern_file = "meta.data|metadata_",
  pattern_dir = NULL,
  pattern_exclude = NULL,
  pattern_order = NULL, # vector of patterns
  column_cluster = "snn_res.", # pattern for cluster columns
  reduction = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2")),
  animation = NULL,
  output_dir = NULL,
  verbose = TRUE,
  ... # clustree_fetch_data parameters
){
  if(verbose) cat("---- Building cluster tree\n** Finding files\n")
  fnames <- if(any(dir.exists(dpath))){
    dnames <- gsub("\\/{2,}", "/", list.files(
      path = dpath[dir.exists(dpath)],
      pattern = pattern_file, all.files = TRUE, full.names = TRUE
    )); fnames <- dnames[!dir.exists(dnames)];
    dnames = dnames[dir.exists(dnames)]; dnames = dnames[!grepl("\\.$", dnames)]
    dnames <- list.files(path = dnames,
      pattern = pattern_file, full.names = TRUE, recursive = TRUE)
    c(dpath[!dir.exists(dpath)], dnames,
      grep(pattern = pattern_file, x = fnames, value = TRUE))
  }else{ dpath[file.exists(dpath)] }
  if(!is.null(pattern_exclude))
    fnames <- fnames[!grepl(pattern = pattern_exclude, x = fnames)]
  fnames <- normalizePath(fnames)

  if(is.null(output_dir))
    output_dir <- paste0(head(dirname(fnames), 1), '/clustree/')

  if(!is.null(pattern_order)){
    if(verbose) cat("Ordering\n")
    fnames <- if(pattern_order == "time"){
      fnames[order(file.info(fnames)$ctime)]
    }else{
      unlist(x = lapply(
        X = pattern_order,
        FUN = function(x) fnames[grepl(pattern = x, x = fnames)] ))
    }; fnames = unique(fnames)
  }

  if(verbose){
    cat("Writing output to:", output_dir, '\n')
    cat("Taking files:", fnames, sep = '\n')
  }; tvar <- !dir.exists(output_dir) && grepl("\\/$", output_dir)
  if(tvar) dir.create(output_dir)

  tmp = clustree_fetch_data(
    fnames, animation = animation, column_cluster = column_cluster,
    reduction = reduction, verbose = verbose, ...)
  summary_clust = tmp$mdata; summary_df = tmp$summary; sufix = tmp$type
  write.table(summary_df, file = paste0(output_dir, sufix, 'cluster_order.txt'),
    quote = FALSE, sep = "\t", row.names = FALSE)
  fname <- paste0(output_dir, sufix, 'cluster_annnotation.txt')
  write.table(summary_clust, file = fname, sep = "\t", row.names = FALSE)

  # Animation
  # https://www.rust-lang.org/tools/install
  # curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  # devtools::install_github("r-rust/gifski")
  # devtools::install_github('thomasp85/gganimate')
  # https://stackoverflow.com/questions/56389470/convert-multiple-png-to-gif-as-an-animation-in-r
  tvar <- requireNamespace("gganimate", quietly = TRUE) && !is.null(animation)
  if(tvar){
    if(verbose) cat("Creating animation\n"); library(gganimate)
    reductions <- ident_dimensions(colnames(summary_clust), animation[1])
    if(!is.na(animation[2])){
      reductions <- reductions[grepl(animation[2], names(reductions))]
    }
    set.seed(27); subsette = sample(rownames(summary_clust), 1000)
    summary_anim_list <- lapply(
      X = 1:2,
      FUN = function(x){
        tmp <- colnames(summary_clust)
        tvar <- tmp[!tmp %in% sapply(reductions, "[", x)]
        ddf = reshape2::melt(summary_clust[subsette, ], id.vars = tvar)
        tmp <- !colnames(ddf) %in% tvar
        colnames(ddf)[tmp] <- paste0(c("step_", "dim_"), x)
        return(ddf[, !grepl(animation[1], colnames(ddf))])
    })
    summary_anim <- joindf(
      summary_anim_list[[1]], summary_anim_list[[2]]
    )
    summary_anim$dim_1 <- as.numeric(summary_anim$dim_1)
    summary_anim$dim_2 <- as.numeric(summary_anim$dim_2)
    tvar <- grep(column_cluster[1], colnames(summary_anim),
      value = TRUE, invert = TRUE)
    df2plot <- reshape2::melt(summary_anim, id.vars = tvar)
    df2plot$Solution <- ident_combine(df2plot, c("step_1", "variable"))
    df2plot$Cluster <- df2plot$value
    aesy <- aes_string(x = "dim_1", y = "dim_2", color = "Cluster")
    p <- ggplot(data = df2plot, mapping = aesy) +
      geom_point() + theme_classic() + transition_states(states = Solution)
    if(requireNamespace("glue", quietly = TRUE)){
      library(glue)
      p <- p + ggtitle('{closest_state}', subtitle = '{frame} of {nframes}')
    }
    gganimate::anim_save(
      filename = paste0(sufix, ifelse(sufix == "", "", "_"), 'animation.gif'),
      animation = gganimate::animate(
        plot = p, duration = nlevels(df2plot$Solution),
        height = 7, width = 7, units = "in", res = 150),
      path = output_dir
    )
  }

  pp <- list()
  if(verbose){
    tvar <- sum(grepl(column_cluster[1], colnames(summary_clust)))
    cat("** Building trees...\nLayers:", tvar, "\n")
  }
  pp[['tree']] <- clustree(x = summary_clust, prefix = column_cluster[1])
  if(sufix == "markers_"){
    if(verbose) cat("Layouts for markers\n")
    pp[['sugi_tree']] <- clustree(
      summary_clust, prefix = column_cluster[1], layout = "sugiyama")
    pp[['sugi_nc_tree']] <- clustree(
      summary_clust, prefix = column_cluster[1], layout = "sugiyama",
      use_core_edges = FALSE)
  }
  if(!sufix == "markers_" && (length(reduction) > 0)){
    if(verbose) cat("On dimentional reductions\n")
    for(i in 1:length(reduction)){
      pname <- paste0(names(reduction[i]), '_overlay'); if(verbose) cat(pname)
      tvar <- !all(reduction[[i]] %in% colnames(summary_clust))
      if(tvar || pname %in% names(pp)){
        if(verbose) cat(" skip\n"); next
      }; if(verbose) cat('!\n')
      if(verbose) cat("- ", pname, "\n")
      print(ncol(summary_clust))
      pp[[pname]] <- clustree_overlay(
        summary_clust,
        prefix = column_cluster[1],
        x_value = reduction[[i]][1], y_value = reduction[[i]][2],
        label_nodes = ncol(summary_clust) <= (1 + 2 + 4) # barcode, x/y values
      )
    }
  }

  p_size <- lapply(gsub(".*_", "", names(pp)), function(x){
    if(grepl('^tree', x)) return(c(16, 15))
    if(grepl('^overlay', x)) return(c(14, 14))
  })

  if(verbose) cat("** Plotting\n")
  for(i in 1:length(pp)){
    fname <- paste0(output_dir, sufix, names(pp[i]), '.pdf')
    cat(names(pp[i]), "\n")
    # if(file.exists(fname) && i > 1) next
    pdf(fname, width = p_size[[i]][1], height = p_size[[i]][2]);
    print(pp[[i]])
    dev.off()
  }
}

clustree_fetch_data = function(
  fnames,
  pattern_pc = NULL,
  pattern_pct = ".*_([0-9]{1,})p_.*|.*_pct([0-9]{1,})_pc.*",
  column_cluster = "snn_res.", # pattern for cluster columns
  column_lib = "origlib", # library column so cell barcode names are consistent
  reduction = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2")),
  make_unique_columns = FALSE,
  sobject = NULL,
  animation = NULL,
  verbose = TRUE
) {
  if(verbose) cat("------ Fetching tree data ------\n")
  if(is.null(pattern_pc)){
    pattern_pc = c(".*kers.([0-9]{1,})PCs.*|.*_pc([0-9]{1,})_res.*",
      ".*ta_([0-9]{1,})PC.*|.*_pc([0-9]{1,})\\.rds")
    pattern_pc = pattern_pc[ifelse(grepl('\\/markers\\/|dgea', fnames[1]), 1, 2)]
  }
  if(verbose) cat("** Reading files (c = clusters, m = markers)\n")
  metas <- lapply(
    X = fnames,
    FUN = function(x){
      pcs <- gsub(pattern = pattern_pc, replacement = "\\1\\2", x = x)
      pct <- gsub(pattern = pattern_pct, replacement = "\\1\\2", x = x)
      y <- readfile(x, row.names = 1, v = FALSE)
      if(all(c('cluster', 'gene') %in% colnames(y))){ # markers
        if(verbose) cat("m")
        # choose column_collapse or summarise or any other you can find
        # check code from /home/ciro/scripts/handy_functions/R/clustering_clustree.R
        yy <- column_collapse(
          metab = y[, c('cluster', 'gene_name')],
          rname = 'gene_name',
          verbose = verbose
        ); rownames(yy) <- yy$gene_name
        yy$ngenes <- table(yy[, 'cluster'])[yy[, 'cluster']]
        splat <- strsplit(yy[, 'cluster'], "&")
        yy$nclusters <- sapply(splat, length)
        nlayers <- 1 #max(yy$nclusters)
        newy <- data.frame(sapply(1:nlayers, function(iclust){
          gsub(" {1,}", "", sapply(splat, "[", iclust))
        }), row.names = rownames(yy))
        colnames(newy) <- paste0(column_cluster[1], 1:ncol(newy))
        if(verbose) cat(paste0("c", ncol(y)))
        y <- cbind(yy, newy)
      }else{
        if(verbose) cat("c")
        if(column_lib %in% colnames(y)){
          tmp <- gsub('\\-[0-9]{1,}', '-', rownames(y))
          rownames(y) <- paste0(tmp, y[, column_lib])
        }; tvar <- paste0(c(column_cluster[1], animation[1]), collapse = "|")
        y <- y[, grepl(tvar, colnames(y)), drop = FALSE]
      };
      if(!is.na(column_cluster[2]))
        y <- y[, grepl(column_cluster[2], colnames(y)), drop = FALSE]
      tvar <- paste0(column_cluster[1], pct, pcs)
      colnames(y) <- gsub(column_cluster[1], tvar, colnames(y))
      if(!is.null(animation)){
        tvar <- grepl(animation[1], colnames(y))
        colnames(y)[tvar] <- paste0(pct, pcs, colnames(y)[tvar])
      }
      return(y)
  }); if(verbose) cat("\n")
  is_marker <- all(c("nclusters", "ngenes") %in% colnames(metas[[1]]))

  if(verbose) cat("** Creating resolution ID\n")
  metas <- lapply(1:length(metas), function(x){
    y <- metas[[x]]
    y$barcode <- rownames(y)
    if(isTRUE(make_unique_columns)){
      tvar <- paste0(column_cluster[1], x, "0")
      colnames(y) <- gsub(column_cluster[1], tvar, colnames(y))
    }
    return(y)
  }); # str(metas)

  mdata_clust <- Reduce(
    f = function(dtf1, dtf2) dplyr::full_join(dtf1, dtf2, by = 'barcode'),
    x = metas)
  tvar <- c('barcode', colnames(mdata_clust)[colnames(mdata_clust) != 'barcode'])
  mdata_clust <- mdata_clust[, tvar]
  colnames(mdata_clust) <- gsub(".x$|.y", "", colnames(mdata_clust))
  rownames(mdata_clust) <- mdata_clust$barcode
  tvar <- apply(X = mdata_clust, MARGIN = 2, FUN  = as.character)
  mdata_clust <- data.frame(
    tvar, row.names = rownames(mdata_clust),
    check.names = FALSE, stringsAsFactors = FALSE)
  mdata_clust[is.na(mdata_clust)] <- "Missed"

  if(!is.null(sobject) && !is_marker){
    if(is.character(sobject)) sobject <- readfile(sobject)
    tvar <- if(column_lib %in% colnames(sobject@meta.data)){
      tmp <- gsub('\\-.*', '-', rownames(sobject@meta.data))
      paste0(tmp, sobject@meta.data[, column_lib])
    }else{ rownames(sobject@meta.data) }
    if(any(rownames(mdata_clust) %in% tvar)){
      if(verbose) cat("Using object\n")
      ddf <- Seurat::FetchData(
        object = sobject, vars = c(unname(unlist(reduction)), column_lib))
      tvar <- sapply(reduction, function(x) all(x %in% colnames(ddf)) )
      reduction <- reduction[tvar]; if(verbose) str(reduction)
      if(length(reduction) > 0){
        tvar <- gsub('\\-.*', '-', rownames(ddf))
        rownames(ddf) <- paste0(tvar, ddf[, column_lib])
        mdata_clust <- joindf(mdata_clust, ddf)
      }
    }
  }else{ reduction = NULL }
  colnames(mdata_clust) = gsub(
    pattern = paste0(".*(", column_cluster[1], "[0-9]{,1}.*)"),
    replacement = "\\1",
    x = colnames(mdata_clust)
  )
  if(verbose) str(mdata_clust)
  if(verbose && is_marker) cat("%%%%%% These are markers %%%%%%\n")

  # Clustering order
  summary_df <- data.frame(
    N = 1:length(fnames),
    Name = fnames,
    PCT = gsub(pattern_pct, '\\1\\2', fnames),
    PC = gsub(pattern_pc, '\\1\\2', fnames)
  )
  if(verbose) cat("------ -------- ---- ---- ------\n")
  return(list(
    mdata = mdata_clust, summary = summary_df,
    type = if(is_marker) "markers_" else ""
  ))
}
