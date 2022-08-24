#!/usr/bin/R

######################
# Clustering: Seurat #
######################

optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
  optparse::make_option(
    opt_str = c("-p", "--percent"), type = "numeric",
    help = "Percentage of variance."
  ),
  optparse::make_option(
    opt_str = c("-n", "--n_comp"), type = "numeric", default = 50,
    help = "Total number of components to explore."
  ),
  optparse::make_option(
    opt_str = c("-c", "--chosen_comp"), type = "numeric",
    help = "Chosen components. Default: chooses with the 'get_elbow' function."
  ),
  optparse::make_option(
    opt_str = c("-r", "--resolution"), type = "character",
    help = paste0("Resolution(s). Usually when running just the markers.\n\t\t",
    "It doesn't replace resolution in the YAML.")
  ),
  optparse::make_option(
    opt_str = c("-m", "--do_markers"), type = "logical", default = TRUE,
    help = "Find markers."
  ),
  optparse::make_option(
    opt_str = c("-t", "--prefix"), type = "character",
    help = paste0("Prefix for this combination of parameters. Default isn\\t\t",
    "determined from the configuration file.")
  ),
  optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line
opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))
# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "/home/ciro/scripts/clustering/config.yaml"
}
config = yaml::read_yaml(opt$yaml)

if(suppressMessages(require(crayon))){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }

if(opt$verbose){
  cat(red_bold("-----------------------------------------\n"))
  cat(red_bold("----------- Seurat clustering -----------\n"))
  cat(red_bold("-----------------------------------------\n"))
}

#### Digesting parameters ###---------------------------------------------------
config$variable_features$mean.cutoff <- unlist(config$variable_features$mean.cutoff)
config$variable_features$dispersion.cutoff <- as.numeric(unlist(config$variable_features$dispersion.cutoff))
config$dim_reduction$tsne$perplexity <- unlist(config$dim_reduction$tsne$perplexity)
if(is.null(config$variable_features$file)) config$variable_features$file <- "none"
# substitutions
if(is.null(config$dim_reduction$base$n_comp)) config$dim_reduction$base$n_comp <- opt$n_comp
if(is.null(config$dim_reduction$base$chosen_comp) || !is.null(opt$chosen_comp)) config$dim_reduction$base$chosen_comp <- opt$chosen_comp
if(is.null(config$dim_reduction$base$chosen_comp)) config$dim_reduction$base$chosen_comp <- 0
if(is.null(config$variable_features$percent) || !is.null(opt$percent)) config$variable_features$percent <- opt$percent
if(is.null(config$markers$select) || !is.null(opt$resolution)) config$markers$select <- opt$resolution
if(!isTRUE(opt$do_markers)) config$markers$select <- "skipping_markers"
if(is.null(config$markers$select) && isTRUE(opt$do_markers)) config$markers$select <- "snn_res"
#### ######### ########## ###---------------------------------------------------
if(interactive()) str(config)
options(warn = -1)

if(opt$verbose) cat('\n----------------------- Loading functions -----------------------\n')
resources = c(
  "https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/filters.R", # filters_complex
  "/home/ciro/scripts/clustering/R/utilities.R", # get_source_data, get_elbow, complex_object_fetch, markers_summary
  "/home/ciro/scripts/clustering/R/plotting.R", # variable_features_report
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R" # stats_summary_table
)
for(i in resources){ source(i) }
#### ######### ####-------------------------------------------------------------

#### Directories structure ####-------------------------------------------------
output_dir <- paste0(config$output_dir, "/", config$project_name, "/")
output_dir <- gsub("\\/{2,}", "/", output_dir)
if(!grepl("scratch|beegfs", getwd())){
  cat("No scratch folder involved; careful about temp files...\n")
  dir.create(output_dir, recursive = TRUE); setwd(output_dir)
}
cat("Working in:", getwd(), "\n")
# The prefix is based on norm and variable_features from the config
global_prefix <- if(is.null(opt$prefix)) paste0("seurat", create_prefix(config)) else opt$prefix
#### ########### ######### ####-------------------------------------------------
### Log file ####---------------------------------------------------------------
register_log <- !interactive() && !grepl("scratch|beegfs", getwd())
log_history <- ""
if(register_log){
  ahash <- digest::digest(paste0(c(unlist(opt), unlist(config)), collapse = ""), "md5", serialize = FALSE)
  log_file <- paste0(
    output_dir, "/_log", ifelse(is.null(opt$prefix), "", paste0(opt$prefix, "_")), ahash
  );
  if(opt$verbose) cat("Check log file:", log_file, "\n")
  if(!file.exists(log_file)) file.create(log_file)
  out.file <- file(log_file, open = 'wt')
  sink(out.file) ; sink(out.file, type = 'message')
}
if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
### ### #### ####---------------------------------------------------------------

if(opt$verbose) cat('----------------------- Loading dependencies --------------------\n')
packages <- c("Seurat", "ggplot2", "cowplot")
loaded <- lapply(X = packages, FUN = function(x){
  if(opt$verbose) cat("*", x, "\n")
  suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
})
theme_set(theme_cowplot())

if(opt$verbose) cat('\n\n----------------------- Parameters ------------------------------\n')
str(opt)
str(config)

#### Getting data ####----------------------------------------------------------
filters_file <- paste0(".seurat_filters.rds")
these_filters <- if(file.exists(filters_file)) readRDS(filters_file)
init_object_name <- paste0(".object_stem_", global_prefix, ".rds")
# When you ask for new filters OR the object doesn't exists
if(!isTRUE(all.equal(these_filters, config$filtering)) || !file.exists(init_object_name)){
  if(opt$verbose) cat('----------------------- Loading data ----------------------------\n')
  Sys.time() # Find 10X directory (barcodes, gene_names, counts), Seurat object, CSV file, TXT, RDS
  expr_data <- get_source_data(
    xpath = config$input_expression,
    pj = config$project_name,
    metadata = config$metadata[1],
    merge_counts = grepl("meco|merge", config$project_name), # need to adjust function
    verbose = opt$verbose
  )
  Sys.time()
  config$input_expression <- expr_data$source
  meta_data <- remove.factors(expr_data$mdata)
  expr_data <- expr_data$edata

  addmetadataf <- if(length(config$metadata[-1]) > 0) config$metadata[-1] else "no_file"
  addmetadataf <- path.expand(addmetadataf)
  tvar <- file.exists(addmetadataf)
  if(any(tvar)){
    addmetadataf <- addmetadataf[tvar]
    if(opt$verbose) cat("Adding metadata from:", addmetadataf, sep = "\n")
    for(addmetadataf_i in addmetadataf){
      addannot <- remove.factors(readfile(addmetadataf_i))
      print(dim(meta_data))
      # partial matching problem addressed in /home/ciro/covid19/scripts/partial_matching.R
      meta_data <- joindf(x = meta_data, y = addannot)
    }
    print(dim(meta_data))
  }

  if(opt$verbose) cat('----------------------- Filtering data --------------------------\n')
  filtereddata <- filters_complex(
    mdata = meta_data,
    filters = lapply(names(config$filtering$subset), function(x) c(x, config$filtering$subset[[x]]) ),
    verbose = opt$verbose
  )
  cat("Preserving:", nrow(filtereddata[[1]]), "/", nrow(meta_data), "samples/cells\n")
  meta_data <- filtereddata[[1]]; rm(filtereddata)
  expr_data <- expr_data[, rownames(meta_data)]

  tvar <- grep(pattern = "orig.ident|seurat_clusters|^RNA|RNA$|^dim_", x = colnames(meta_data), value = TRUE, invert = TRUE)
  if(opt$verbose) cat("Keeping", length(tvar), "of", ncol(meta_data), "columns\n")
  meta_data <- meta_data[, tvar, drop = FALSE]
  filts_names = gsub("!|\\(|\\)", "", unlist(strsplit(unlist(config$filtering$subset), " ")))
  filts_names = unique(filts_names[filts_names %in% colnames(meta_data)])
  if(length(filts_names) > 0){
    void <- lapply(setNames(nm = filts_names), function(i){
      if(is.character(meta_data[, i])) table(meta_data[, i], useNA = 'always') else summary(meta_data[, i])
    }); print(void)
  }

  saveRDS(object = config$filtering, file = filters_file)
  log_history <- paste0(log_history, "Filt")
}

#### Seurat Pipeline ####-------------------------------------------------------
if(opt$verbose) cat('\n\n----------------------- Seurat pipeline steps -------------------\n')
if(!file.exists(init_object_name)){
  nSamples_min = ceiling(config$filtering$nSamples_expressed * ncol(expr_data))
  if(opt$verbose){
    cat("Samples/cells:", ncol(expr_data), "\n"); cat("Features:", nrow(expr_data), "\n")
    cat(
      '@ Initialize the Seurat object: keeping genes with >=',
      nSamples_min, 'samples/cells expressing them & cells with >=', 0,' expressed genes)...\n'
    )
  }; timestamp()
  mycells <- CreateSeuratObject(
    counts = expr_data,
    project = config$project_name, # gsub(".PBS.*", "", basename(getwd())),
    min.cells = nSamples_min,
    min.features = 0,
    meta.data = meta_data
  ); timestamp()
  if(opt$verbose){ cat("Samples/cells:", ncol(mycells), "\n");
  cat("Features:", nrow(mycells), "/", nrow(expr_data), "\n") }
  suppressWarnings(try(rm(expr_data, meta_data), silent = TRUE))

  ### Variables to regress ### ---
  regress_cc = if(any(config$regress_var %in% c("S.Score", "G2M.Score", "cellcycle"))){
    if("cellcycle" %in% config$regress_var){
      c("S.Score", "G2M.Score")
    }else{ config$regress_var[config$regress_var %in% c("S.Score", "G2M.Score")] }
  }
  config$regress_var <- config$regress_var[config$regress_var %in% colnames(mycells@meta.data)]
  config$regress_var <- unique(config$regress_var)
  tvar <- sapply(config$regress_var, function(x) length(unique(mycells@meta.data[, x])) )
  config$regress_var <- config$regress_var[tvar > 1]
  tvar <- sapply(config$regress_var, function(x) is.character(mycells@meta.data[, x]) )
  if(any(tvar)){
    if(opt$verbose) cat("Prepare (scaling) character variables for regression\n")
    if(opt$verbose) str(mycells@meta.data[, config$regress_var, drop = FALSE])
    for(i in config$regress_var[tvar]){
      mycells@meta.data$var2regress123 = if(is.numeric(mycells@meta.data[, i])){
        scale(mycells@meta.data[, i])
      }else{ scale(as.numeric(factor(mycells@meta.data[, i]))) }
      colnames(mycells@meta.data) <- gsub("var2regress123", paste0(i, ".reg"), colnames(mycells@meta.data))
    }; config$regress_var[tvar] <- paste0(config$regress_var[tvar], ".reg")
    if(opt$verbose) str(mycells@meta.data[, config$regress_var, drop = FALSE])
  }

  #### ------ -------- ####-----------------------------------------------------
  ### Filter highly variable features
  hvg_final <- rownames(mycells)
  tvar <- config$variable_features$file
  if(file.exists(tvar)){
    if(opt$verbose) cat("Subsetting features from file:\n", tvar, "\n")
    file_con <- file(description = tvar, open = "r")
    these_feats <- readLines(con = file_con); close(file_con)
    hvg_final <- hvg_final[hvg_final %in% these_feats]
  }

  if(!grepl(config$norm, "sctransform")){
    if(opt$verbose) cat('\n@@@@@@@@@ Normalising...\n'); timestamp()
    mycells <- NormalizeData(
      object = mycells,
      normalization.method = config$norm,
      scale.factor = 10000,
      verbose = opt$verbose
    )
    log_history <- paste0(log_history, "Norm")
  }else if(is.character(regress_cc)){
    if(opt$verbose) cat('\n@@@@@@@@@ SCTransform - no regression\n'); timestamp()
    mycells <- SCTransform(
      object = mycells,
      variable.features.n = config$variable_features$nfeatures,
      vars.to.regress = NULL,
      conserve.memory = FALSE,
      return.only.var.genes = TRUE,
      verbose = opt$verbose
    )
  }; if(opt$verbose) cat("Default Assay:", DefaultAssay(object = mycells), "\n")

  if(!grepl(config$norm, "sctransform")){
    if(opt$verbose) cat('\n@@@@@@@@@ Selecting features...\n'); timestamp()
    mycells <- FindVariableFeatures(
      object = mycells,
      selection.method = config$variable_features$method,
      nfeatures = config$variable_features$nfeatures,
      mean.cutoff = config$variable_features$mean.cutoff,
      dispersion.cutoff = config$variable_features$dispersion.cutoff,
      verbose = opt$verbose
    )
    hvg_df <- HVFInfo(mycells, selection.method = config$variable_features$method)
    disp_n <- colnames(hvg_df)[3]; hvg_df <- hvg_df[order(-hvg_df[, disp_n]), ] # order
    # turns mean and pct filters off in case you gave a file with HVGs
    mean_pct_filters <- !isTRUE(config$variable_features$file_only)
    tmp <- config$variable_features$exclude # Exclude, for example, "^MT" genes
    if(!is.null(tmp) && mean_pct_filters){ # if null, don't do it
      if(opt$verbose) cat("Excluding", tmp, "\n")
      tvar <- grep(pattern = tmp, x = rownames(hvg_df), value = TRUE, ignore.case = TRUE)
      if(opt$verbose)
        cat("N =", length(tvar), "-->", paste0(head(tvar, 10), collapse = ", "), "\n")
      hvg_df <- hvg_df[!rownames(hvg_df) %in% tvar, ]
    }
    if(mean_pct_filters){
      passed <- hvg_df[, 1] > config$variable_features$mean.cutoff[1]
      tvar <- paste("Mean >", config$variable_features$mean.cutoff[1])
      if(opt$verbose) cat(tvar, "filter:", sum(passed), "of", nrow(hvg_df), "\n")
      hvg_df <- hvg_df[passed, ]
    }
    if(!is.null(config$variable_features$percent) && mean_pct_filters){
      if(opt$verbose) cat("Taking", config$variable_features$percent, "%\n")
      hvg_df$cumulative = cumsum(hvg_df[, disp_n])
      hvg_df$cumulative_pct = round(hvg_df$cumulative / sum(hvg_df[, disp_n], na.rm = TRUE) * 100, 2)
      passed <- hvg_df$cumulative_pct <= config$variable_features$percent
      if(opt$verbose) cat("Number of features:", sum(passed), "of", nrow(hvg_df), "\n")
      hvg_df <- hvg_df[passed, ]; log_history <- paste0(log_history, "Pct")
    }
    hvg_final <- intersect(hvg_final, rownames(hvg_df))
    if(opt$verbose) cat("Taking", length(hvg_final), "/", length(VariableFeatures(mycells)), "features\n")
    tmp <- length(intersect(VariableFeatures(mycells), hvg_final))
    if(opt$verbose) cat("Overlap with initial selection:", tmp, "\n")
    VariableFeatures(object = mycells, assay = NULL) <- hvg_final
    if(opt$verbose) cat("Final number:", length(VariableFeatures(mycells)), "features\n")
  }

  if(is.character(regress_cc)){
    # Might not be really necessary: https://github.com/satijalab/seurat/issues/1679
    cat("Getting Cell Cycle scores\n")
    mycells <- CellCycleScoring(
      object = mycells,
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      set.ident = TRUE
    ); config$regress_var = c(config$regress_var, regress_cc)
  }

  if(!grepl(config$norm, "sctransform")){
    if(opt$verbose) cat('\n@@@@@@@@@ Scaling data...\n'); timestamp()
    mycells <- ScaleData(
      object = mycells,
      vars.to.regress = config$regress_var,
      block.size = 2000,
      verbose = opt$verbose
    ); log_history <- paste0(log_history, "Scale")
  }else{
    if(opt$verbose) cat('\n@@@@@@@@@ SCTransform...\n'); timestamp()
    library(sctransform)
    mycells <- SCTransform(
      object = mycells,
      variable.features.n = config$variable_features$nfeatures,
      vars.to.regress = config$regress_var[!config$regress_var %in% "nCount_RNA"],
      conserve.memory = FALSE,
      return.only.var.genes = TRUE,
      verbose = opt$verbose
    ); log_history <- paste0(log_history, "Sctrans")
  }
  if(opt$verbose) cat("Default Assay:", DefaultAssay(object = mycells), "\n")

  #### ------ -------- ####-----------------------------------------------------
  if(opt$verbose){
    cat('\n@@@@@@@@@ Linear dimensional reduction\n')
    cat('Computing:', casefold(config$dim_reduction$base$type, upper = TRUE), '\n')
    cat('Components:', config$dim_reduction$base$n_comp, '\n')
  }; timestamp()
  mycells <- RunPCA(
    object = mycells,
    features = VariableFeatures(mycells),
    npcs = config$dim_reduction$base$n_comp,
    nfeatures.print = 15
  ); log_history <- paste0(log_history, "LR")
  saveRDS(object = mycells, file = init_object_name)
  write(paste(Sys.time(), init_object_name), file = ".seurat_clustering_history", append = TRUE)
}else{
  mycells <- readRDS(init_object_name)
}; suppressWarnings(try(rm(expr_data, meta_data), silent = TRUE))
writeLines(text = init_object_name, con = paste0(".object_init_", global_prefix))

try({
  variable_features_report(
    object = mycells,
    cutoff = config$variable_features$mean.cutoff[1],
    prefix = global_prefix,
    verbose = opt$verbose
  )
})

#### ------ -------- ####-------------------------------------------------------
if(opt$verbose) cat("\n*** Selecting number of components\n")
elbows <- get_elbow(x = mycells@reductions$pca@stdev, threshold = seq(.7, .9, .02), decide = TRUE)
tvar <- table(elbows[[1]])
if(opt$verbose) cat('Elbow(s):', paste0(paste0(names(tvar), "(n=", tvar, ")"), collapse = ", "), '\n')
chosen_comp <- config$dim_reduction$base$chosen_comp[1] # possible iteration
if(chosen_comp < 1) chosen_comp <- elbows$elbow
cat("Selected:", chosen_comp, "\n")

p <- ElbowPlot(object = mycells, ndims = config$dim_reduction$base$n_comp) +
  ggtitle(label= casefold(gsub("_", " ", basename(getwd())), upper = TRUE)) +
  geom_vline(xintercept = elbows$elbows, linetype = "dotted", color = "gray") +
  geom_text(data = elbows$ptable, aes(x = x, y = y, label = N)) +
  geom_vline(xintercept = chosen_comp, linetype = "dashed", color = "red")
fname <- paste0(global_prefix, "_elbow_", chosen_comp, "of", config$dim_reduction$base$n_comp, ".pdf")
if(!file.exists(fname)){ pdf(fname, width = 9); print(p); graphics.off() }

initial_pc <- paste0(global_prefix, '_pc', config$dim_reduction$base$chosen_comp[1])
suffix_comp <- paste0(global_prefix, '_pc', chosen_comp)

#### ------ -------- ####-------------------------------------------------------
### Getting previous analyses info ###
# NOTE: complex_object_fetch is the heart of this step, where you retrieve the whole
# object without duplicating it. It's up to you which combination of parameters
# you use for your downstream analyses
# Files for speeding things up, aka, not re-calculate things: reductions, metadata, graphs, commands
reductionsf <- paste0(".object_reductions_", suffix_comp, ".rds")
graphsf <- paste0(".object_graphs_", suffix_comp, ".rds")
commandsf <- paste0(".object_commands_", suffix_comp, ".rds")
metadataf <- paste0(".object_meta.data_", suffix_comp, ".rds")
mycells <- complex_object_fetch(object_or_dir = mycells, id = paste0("_", suffix_comp), verbose = opt$verbose)
reductions <- mycells@reductions
compute_snn <- length(mycells@graphs) == 0

norm_type <- ifelse(grepl("sctransform", config$norm), "SCT_snn_res.", "RNA_snn_res.")
res_calc <- config$resolution[!paste0(norm_type, config$resolution) %in% colnames(mycells@meta.data)]

if(opt$verbose) cat('\n@@@@@@@@@ Clustering\n'); timestamp()
if(length(res_calc) > 0){
  if(opt$verbose) cat("Resolutions:", paste0(res_calc, collapse = ", "), "\n")
  mycells <- FindNeighbors(
    mycells, dims = 1:chosen_comp, compute.SNN = compute_snn, verbose = opt$verbose
  )
  mycells <- FindClusters(
    object = mycells, resolution = config$resolution, verbose = opt$verbose
  ); log_history <- paste0(log_history, "Clust")
  metadata <- mycells@meta.data; saveRDS(object = metadata, file = metadataf)
  graphs <- try(mycells@graphs); saveRDS(object = graphs, file = graphsf)
}

badpx <- Inf
updated <- FALSE
if(opt$verbose) cat('\n@@@@@@@@@ Non-linear dimensional reduction\n'); timestamp()
if(opt$verbose) cat("Reductions:", paste0(names(mycells@reductions), collapse = ", "), "\n")
pxities <- if(!is.null(config$dim_reduction$tsne)) unique(config$dim_reduction$tsne$perplexity)
n_neis <- if(!is.null(config$dim_reduction$umap)) unique(config$dim_reduction$umap$n.neighbors)
for(i in 1:max(length(pxities), length(n_neis))){
  if(!is.null(config$dim_reduction$tsne)){ ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    px <- ifelse(is.na(pxities[i]), tail(pxities, 1), pxities[i])
    if(grepl("auto|sqrt_n", px)) px <- round(sqrt(ncol(mycells)))
    dim_reduction_id <- paste0("tsne_pc", chosen_comp, "_px", px)
    cat(">>>>>>>>>>>>>>>>>>>>>>", dim_reduction_id, "\n")
    if(px < badpx && !dim_reduction_id %in% names(reductions)){
      mycells@reductions <- mycells@reductions[!names(mycells@reductions) %in% "tsne"]
      tmp <- try({ # perplexity: how well a probability distribution/model predicts a sample (lower better)
        mycells <- RunTSNE(
          object = mycells, reduction = config$dim_reduction$base$type,
          dims = 1:chosen_comp,
          perplexity = as.numeric(px),
          tsne.method = "FIt-SNE", fast_tsne_path = '/mnt/BioHome/ciro/bin/FIt-SNE2/bin/fast_tsne'
        ); "success"
      })#, silent = TRUE)
      if(class(tmp) == 'try-error'){
        str(tmp); badpx <- px <- round(sqrt(ncol(mycells)))
        if(opt$verbose){
          cat('\nTrying perplexity ~ N^(1/2)\n\n')
          cat("Reductions:", paste0(names(mycells@reductions), collapse = ", "), "\n")
        }
        mycells <- RunTSNE(
          object = mycells, reduction = config$dim_reduction$base$type,
          dims = 1:chosen_comp,
          perplexity = round(sqrt(ncol(mycells))),
          tsne.method = "FIt-SNE", fast_tsne_path = '/mnt/BioHome/ciro/bin/FIt-SNE2/bin/fast_tsne'
        ); rm(tmp)
      }else{ rm(tmp) }
      dim_reduction_id <- paste0("tsne_pc", chosen_comp, "_px", px); timestamp()
      cat(">>>>>>>>>>>>>>>>>>>>>>", dim_reduction_id, "calculated\n")
      eval(expr = parse(text = paste0("reductions$", dim_reduction_id, " <- mycells@reductions$tsne")))
      updated <- TRUE
    }
  }
  if(!is.null(config$dim_reduction$umap)){ ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nn <- ifelse(is.na(n_neis[i]), tail(n_neis, 1), n_neis[i])
    dim_reduction_id <- paste0("umap_pc", chosen_comp, "_nei", nn, "_dist", config$dim_reduction$umap$min.dist)
    cat(">>>>>>>>>>>>>>>>>>>>>>", dim_reduction_id, "\n")
    if(!dim_reduction_id %in% names(reductions)){
      mycells@reductions <- mycells@reductions[!names(mycells@reductions) %in% "umap"]
      mycells <- RunUMAP(
        object = mycells, reduction = config$dim_reduction$base$type,
        dims = 1:chosen_comp,
        n.neighbors = nn,
        min.dist = config$dim_reduction$umap$min.dist,
        verbose = opt$verbose
      ); timestamp()
      cat(">>>>>>>>>>>>>>>>>>>>>>", dim_reduction_id, "calculated\n")
      eval(expr = parse(text = paste0("reductions$", dim_reduction_id, " <- mycells@reductions$umap")))
      updated <- TRUE
    }
  }

  if(updated) saveRDS(object = reductions, file = reductionsf)
};
if(updated){
  commands <- mycells@commands; saveRDS(object = commands, file = commandsf)
  mycells@reductions <- reductions
  metadata <- mycells@meta.data
  for(dr in names(reductions)){
    mycells@reductions[[gsub("_.*", "", dr)]] <- mycells@reductions[[dr]]
    y <- FetchData(mycells, vars = paste0(mycells@reductions[[dr]]@key, 1:2))
    colnames(y) <- casefold(paste0("dim_", colnames(y), gsub(".*_pc[0-9]{1,}|pca", "", dr)))
    tvar <- !colnames(y) %in% colnames(metadata)
    if(any(tvar)) metadata <- cbind(metadata, y[, tvar, drop = FALSE])
  }
  str(metadata)
  saveRDS(object = metadata, file = metadataf)
}
writeLines(
  text = paste(config$dim_reduction$base$chosen_comp[1], ">", chosen_comp),
  con = paste0(".dr_", initial_pc)
)

#### ------ -------- ####-------------------------------------------------------
if(opt$verbose) cat('\n@@@@@@@@@ Markers\n'); timestamp()
nresolutions <- grep("RNA_snn_res.|SCT_snn_res.", colnames(mycells@meta.data), value = TRUE)
nresolutions <- grep(config$markers$select, nresolutions, value = TRUE)
tvar <- as.numeric(config$markers$select)
if(!is.na(tvar))
  nresolutions <- grep(paste0("s\\.", as.character(tvar)), nresolutions, value = TRUE)
if(opt$verbose) cat("Based on:", config$markers$select, "\n")
if(length(nresolutions) > 0){
  if(opt$verbose) cat("Resolutions:", paste0(nresolutions, collapse = ", "), "\n")
  log_history <- paste0(log_history, "Mark")
}
for(nres in nresolutions){
  Idents(mycells) <- nres
  if(opt$verbose){
    cat(">>>>>>>>>>>>>>>>>>>>>>", nres, "\n")
    print(table(Idents(mycells)))
  }
  suffix_res <- paste0(c(suffix_comp, gsub("RNA_snn_res.|SCT_snn_res.", "res", nres)), collapse = "_")
  if(nlevels(mycells@meta.data[, nres]) > 1){
    dir.create(suffix_res, showWarnings = FALSE)
    dgea_seurat(
      object = mycells,
      results_prefix = paste0(suffix_res, "/dgea_", config$markers$test),
      config_markers = config$markers,
      cluster_column = nres,
      return_table = FALSE,
      verbose = TRUE
    )
  }else{ cat("NO TEST PERFORMED\n") }
  fname <- paste0(".markers_", initial_pc, gsub("RNA_snn_res.|SCT_snn_res.", "_res", nres))
  writeLines(text = suffix_res, con = fname)
}; if(opt$verbose) cat(log_history, "\n")
#### ###### ######## ####-------------------------------------------------------

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
if(register_log){
  sink(type = 'message')
  sink()
  new_name <- gsub(paste0("(.*_log_)", ".*(_", ahash, ")"), paste0("\\1", log_history, "\\2"), log_file)
  try(file.rename(from = log_file, to = new_name))
}
