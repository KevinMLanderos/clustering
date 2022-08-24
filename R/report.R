#!/usr/bin/R

######################
# Clustering: report #
######################

optlist <- list(
  optparse::make_option(
    opt_str = c("-p", "--path"), type = "character", default = "./",
    help = "Path to clustering results."
  ),
  optparse::make_option(
    opt_str = c("-c", "--do_clusters"), type = "logical", default = TRUE,
    help = "Perform report for clustering solutions."
  ),
  optparse::make_option(
    opt_str = c("-m", "--do_markers"), type = "logical", default = TRUE,
    help = "Perform report for markers."
  ),
  optparse::make_option(
    opt_str = c("-i", "--include"), type = "character",
    help = "Pattern of files to include."
  ),
  optparse::make_option(
    opt_str = c("-e", "--exclude"), type = "character",
    help = "Pattern of files to exclude."
  ),
  optparse::make_option(
    opt_str = c("-f", "--config_file"), type = "character",, default = "config.yaml",
    help = "Configuration file used for clustering."
  ),
  optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

optparse <- optparse::OptionParser(option_list = optlist)
opt <- optparse::parse_args(optparse)
if(interactive()){ # Example/manually
  opt$path = "/home/ciro/large/pbtumor/results/clustering/deprecated/CD45pCD3p_xDoublets"
}

if(suppressMessages(require(crayon))){
  cyan = crayon::cyan; red_bold = crayon::red$bold
}else{ cyan = red_bold = c }

if(opt$verbose){
  cat(red_bold("-----------------------------------------\n"))
  cat(red_bold("----------- Clustering report -----------\n"))
  cat(red_bold("-----------------------------------------\n"))
}

if(!isTRUE(opt$do_clusters) && !isTRUE(opt$do_markers)){
  cat("You indicated that clusters and markers reports shouldn't be generated...\n");
  q(save = "no")
}

setwd(opt$path)
cat("Working in:", getwd(), "\n")
if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();

if(opt$verbose) cat(cyan("----------- Loading dependencies --------\n"))
packages <- c("ggplot2", "cowplot", "dplyr", "Seurat")
loaded <- lapply(X = packages, FUN = function(x){
  if(opt$verbose) cat("*", x, "\n")
  suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
})
theme_set(theme_cowplot())

if(opt$verbose) cat(cyan("----------- Loading functions -----------\n"))
resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # filters_columns is.file.finished show_commas (cluster_reports)
  "/home/ciro/scripts/clustering/R/plotting.R", # cluster_reports
  "/home/ciro/scripts/clustering/R/utilities.R", # get_top_n, marker_report
  "/home/ciro/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/handy_functions/devel/plots.R" # plot_pct getlegend mytheme
)
for(i in resources){ source(i) }

if(isTRUE(opt$do_clusters)){
  if(opt$verbose) cat(cyan("----------- Clustering reports ----------\n"))
  file_names = list.files(pattern = "meta.data", all.files = TRUE)
  if(!is.null(opt$include)) file_names = file_names[grepl(opt$include, file_names)]
  if(!is.null(opt$exclude)) file_names = file_names[!grepl(opt$exclude, file_names)]
  for(mdata_f in file_names){
    cat(mdata_f, "\n")
    mdata <- readfile(mdata_f)
    void <- cluster_reports(
      metadata = mdata,
      config_file = opt$config_file,
      prefix = gsub(".*meta.data_|.rds", "", mdata_f),
      verbose = opt$verbose
    )
    tool_name <- paste0(".", gsub(".*meta.data_([A-z]{1,})_.*", "\\1", mdata_f), "_finished_components")
    writeLines(text = ".finished", con = tool_name)
  }

  source("/home/ciro/scripts/clustering/R/clustree.R")
  void <- try(generate_tree(
    dpath = getwd(),
    # sobject = tail(list.files(pattern = "stem.*rds", all.files = TRUE), 1),
    # animation = "dim_umap_",
    output_dir = "clustree/"
  ))
}

if(isTRUE(opt$do_markers)){
  if(opt$verbose) cat(cyan("----------- Markers reports -------------\n"))
  opt$avg_logFC <- 0.25; opt$p_val_adj <- 0.05
  file_names = list.files(
    pattern = "_summary_stats|dgea_MAST.rds",
    full.names = TRUE, recursive = TRUE)
  if(!is.null(opt$include)) file_names = file_names[grepl(opt$include, file_names)]
  if(!is.null(opt$exclude)) file_names = file_names[!grepl(opt$exclude, file_names)]
  is_filtered <- grepl("_summary_stats", file_names)
  if(any(is_filtered)) file_names = file_names[is_filtered]
  for(cmarkers_f in file_names){
    cat(cmarkers_f, "\n")
    cmarkers <- readfile(cmarkers_f, row.names = 1)

    cmarkers_f1 <- if(!any(is_filtered)){
      tvar <- cmarkers$avg_log >= opt$avg_logFC & cmarkers$p_val_adj <= opt$p_val_adj
      cmarkers <- cmarkers[tvar, ]
      paste0(gsub(paste0("\\.", tools::file_ext(cmarkers_f)), "", cmarkers_f),
        "_fc", opt$avg_logFC, "_padj", opt$p_val_adj, ".csv")
    }else{ cmarkers_f }

    void <- marker_report(
      markers_df = cmarkers,
      file = cmarkers_f1,
      verbose = opt$verbose
    );
    tool_name <- paste0(".", gsub("..([A-z]{1,})_.*", "\\1", cmarkers_f), "_finished_markers")
    writeLines(text = ".finished", con = tool_name)
  }
}
if(opt$verbose) cat('----------------------- ------- ------- -------------------------\n')

if(opt$verbose) {
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
