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
    }; tvar <- cmarkers$avg_log >= config_markers$avg_logFC
    tmp <- cmarkers$p_val_adj <= config_markers$p_val_adj
    cmarkers <- cmarkers[which(tvar & tmp), ]
    if(!file.exists(meansf) && !is.null(cluster_column)){
      cmarkers <- markers_summary(
        marktab = cmarkers,
        annot = object@meta.data,
        datavis = expm1(object@assays$RNA@data) * 100,
        cluster_column = cluster_column,
        datatype = "SeuratNormalized",
        verbose = verbose
      ); write.csv(cmarkers, file = meansf)
    }
  }; if(isTRUE(return_table)) return(cmarkers)
  return(invisible(x = NULL))
}
