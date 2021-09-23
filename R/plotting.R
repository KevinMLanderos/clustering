#!/usr/bin/R

###################
# Seurat plotting #
###################

# This functions are designed to be help with visualisation of Seurat object

qc_violin <- function(
  dat,
  xax = "Data",
  yax = colnames(dat)[1],
  lb_filt = NULL,
  hb_filt = NULL,
  filtby = NULL
){
  if(!"Data" %in% colnames(dat)) dat$Data <- "SET"
  if(!is.factor(dat[, xax])){
    dat[is.na(dat[, xax]), xax] <- "NA"
    dat[, xax] <- factor(dat[, xax])
  }
  dat <- dat[-which.max(dat[, yax]), ] # when the violin is squashed

  if(!any(is.null(c(lb_filt, hb_filt)))){
    dat$filter <- dat[, yax] >= lb_filt[yax] & dat[, yax] <= hb_filt[yax]
    pcts <- paste0(table(dat[dat$filter, xax]), "/", table(dat[, xax]))
    names(pcts) <- levels(dat[, xax])
    dat$pct <- pcts[as.character(dat[, xax])]
  }
  p <- simple_violin(dat = dat, x = xax, y = yax) +
    theme(
      axis.text.x = element_text(hjust = 1, angle = 45),
      legend.position = "none",
      axis.title.x = element_blank()
    ) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  if('pct' %in% colnames(dat)){
    ymax <- max(dat[, yax], na.rm = TRUE)
    ll <- c(lb_filt[yax], min(dat[, yax], na.rm = TRUE))
    hl <- c(hb_filt[yax], ymax); coly <- ifelse(isTRUE(filtby[yax]), 'red', '#808080')
    p <- p + geom_hline(yintercept = max(ll), linetype = "dashed", color = coly) +
      geom_hline(yintercept = min(hl), linetype = "dashed", color = coly) +
      geom_text(
        aes_string(label = 'pct', y = ymax), # position = position_dodge(0.9),
        nudge_x = ifelse(nlevels(dat[, xax]) < 3, 0.3, 0),
        vjust = 1.6, angle = 45
      )
  }
}

simple_violin <- function(
  dat,
  xax = 1,
  yax = 2,
  sample_it = TRUE
){ # from filters sample_even
  if(is.numeric(xax)) xax <- colnames(dat)[xax]
  if(is.numeric(yax)) xax <- colnames(dat)[yax]

  # tvar <- sample(1:nrow(dat), min(ifelse(nrow(dat) < 20000, nrow(dat)/3, nrow(dat)/5), 15000))
  tvar <- if(isTRUE(sample_it)){ set.seed(27); sample_even(dat, xax, maxln = -1000, v = FALSE) }else rownames(dat)
  dat <- dat[order(dat[, xax]), ]
  d2show <- dat[tvar, ]

  aesy <- aes_string(x = xax, y = yax, fill = xax)
  dp <- ggplot(dat, aesy) +
    geom_jitter(
      data = d2show, mapping =  aes_string(x = xax, y = yax), inherit.aes = FALSE,
      shape = 16, position = position_jitter(0.2), color = 'black'
    ) +
    geom_violin(aes_string(color = xax), trim = TRUE, scale = 'width', alpha = 0.8)+
    geom_boxplot(width=0.1, fill = "white", alpha = 0.25, outlier.shape = NA)
  if(length(unique(dat[, xax])) < 10){
    dp <- dp + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1")
  }
  dp
}

#' Scatter plot of features and their thresholds
#'
#' This function creates a scatter plot of two (third optional) features or
#' variables with their respective thresholds
#'
#' @param data Annotation (data.frame) or object containing the variables
#' @param variables Vector of variables to plot.
#' @param thresh Vector/table/data.frame with 'variables' columns and two values of
#' @param verbose Show progress.
#' @keywords ROC
#'
#' @return Returns a data.frame with the signature (generates a lot of reports)
#'
#' @importFrom
#'
#' @export
#'
#' @examples
#' feature_scatter(data = meta_data, thresholds = summ[, qc_filters$apply], prefix = "0_qc")
#'

feature_scatter <- function(
  data,
  variables,
  thresh = NULL, # matrix/data.frame: thresholds rows and qc variables as as columns
  verbose = FALSE
){
  if(casefold(class(data)) == "seurat"){
    tvar <- unique(c(variables, colnames(data@meta.data)))
    data <- FetchData(data, vars = tvar); gc()
  }

  # p <- plot_corr( # from handy_functions/devel/plotting.R
  #   df = data, var1 = variables[1],  var2 = variables[2], addvar = variables[3],
  #   add_line = FALSE, log2t = 'auto', return_plot = TRUE, v = verbose
  # )
  my_mapping = if(is.na(variables[3])){
    aes_string(x = variables[1], y = variables[2])
  }else{ aes_string(x = variables[1], y = variables[2], color = variables[3]) }
  p <- ggplot(data = data, mapping = my_mapping) +
    geom_point(size = 0.7, shape = 20) +
    theme_minimal()
  scale_x_fun = if(diff(range(data[, variables[1]])) > 100)
    scale_x_log10(breaks = scales::pretty_breaks(n = 10))
  else scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  scale_y_fun = if(diff(range(data[, variables[2]])) > 100)
    scale_y_log10(breaks = scales::pretty_breaks(n = 10))
  else scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  p <- p + scale_x_fun + scale_y_fun +
    theme(axis.text.x = element_text(hjust = 1, angle = 90))

  if(!is.null(thresh)){ # when no threshold is given
    if(verbose) cat("Adding boundaries\n")
    for(i in seq(1, nrow(thresh), 2)){
      p <- p + annotate(geom = "rect",
        xmin = thresh[i, variables[1]], xmax = thresh[i + 1, variables[1]],
        ymin = thresh[i, variables[2]], ymax = thresh[i + 1, variables[2]],
        alpha = 0.1, color = ifelse(i + 1 == nrow(thresh), 'red', 'black'), fill = 'white'
      )
    }
  }
  p; #graphics.off()
}

qc_scatters <- function(
  data,
  thresholds,
  theme_extra = NULL,
  prefix = "0_qc",
  pdf_width = 10, pdf_height = 10,
  verbose = FALSE
){
  qc_variables = if(is.null(dim(thresholds))) thresholds else colnames(thresholds)
  apfilters <- unique(c(qc_variables, "final_STAR_counts"))
  apfilters <- gtools::combinations(
    length(apfilters), r = 2,
    v = apfilters,
    set = TRUE, repeats.allowed = FALSE
  )
  apfilters_names = c()
  for(i in 1:nrow(apfilters)){
    axes <- c(apfilters[i, ], "uniquely_mapped_reads_perc") # var. of interes
    axes <- unique(c(axes, "percent.mt")) # var. of interes
    axes <- axes[axes %in% colnames(data)]
    axes <- axes[axes %in% qc_variables]
    tvar <- paste0(axes, collapse = "_") %in% apfilters_names
    if(length(axes) < 2 || tvar) next
    apfilters_names <- c(apfilters_names, paste0(axes, collapse = "_"))
    if(verbose) cat(show_commas(axes))
    p <- try(feature_scatter(
      data = data,
      variables = axes,
      thresh = if(!is.null(dim(thresholds))) thresholds[, axes],
      verbose = !TRUE
    ), silent = TRUE)
    if(class(p)[1] != "gg"){
      if(verbose) cat(" - failed\n"); next
    }; if(verbose) cat(" - plotting\n")
    fname <- paste0(paste0(c(prefix, axes), collapse = "_"), ".pdf")
    if(is.file.finished(fname)) next
    if(!is.null(theme_extra)) p <- p + theme_extra
    pdf(fname, width = pdf_width, height = pdf_height)
    print(try(p, silent = TRUE))
    graphics.off()
  }
}

ncounts_outsies <- function(
  object,
  redu = c("PC_1", "PC_2", "nFeature_RNA", "nCount_RNA" , "percent.mt")
){
  tvar <- if(casefold(class(object)) == "seurat")
    FetchData(object = object, vars = redu) else object[, redu]
  tvar$cell_name <- rownames(tvar); x <- order(abs(tvar[, 1]))
  y <- order(abs(tvar[, 2])); outsies <- tvar[unique(c(tail(x), tail(y))), ];
  rm(x, y)
  void <- ggplot(
    data = tvar,
    mapping = aes_string(redu[1], redu[2])) +
    labs(x = redu[1], y = redu[2]) + theme_classic()
  nplots <- lapply(redu[-c(1:2)], function(x){
    void + geom_point(aes_string(colour = x), size = 0.7) +
      scale_colour_gradientn(colours = rainbow(7)) +
      ggrepel::geom_text_repel(
        data = outsies, aes(label = cell_name), size = 0.75)
  })
  print(cowplot::plot_grid(plotlist = nplots))
  return(outsies)
}

## Visualise fit
mean_var <- function(object){
  hvf.info <- HVFInfo(object, selection.method = "vst")
  not.const <- hvf.info$variance > 0
  fit <- loess(
     formula = log10(x = variance) ~ log10(x = mean),
     data = hvf.info[not.const, ],
     span = 0.3
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  not.const <- sample(which(not.const), min(c(100, sum(not.const) * 0.01)))
  p <- ggplot(hvf.info, aes(log10(mean), log10(variance))) + geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE, color = 'red', method = "loess") +
    geom_point(
      data = hvf.info[not.const, ],
      mapping = aes(log10(mean), log10(variance.expected)),
      alpha = 0.6, size = 0.1, color = "green") +
    labs(title = 'Global mean-variance LOESS')
  return(p )
}

cluster_reports <- function(
  metadata,
  confounders = "orig",
  resolutions = "_snn_res.",
  norm_type = "RNA_snn_res.|SCT_snn_res.",
  qc_names = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  prefix = NULL,
  gcols = NULL,
  prop.normalise = TRUE,
  config_file = "config.yaml",
  verbose = FALSE
){ # make_grid from devel/utilities
  if(verbose) cat("  ------- Plotting clusters -------\n")
  # Getting variables
  if(any(!confounders %in% colnames(metadata))){
    confounders_tag <- confounders
    confounders <- filters_columns(metadata,
      include = confounders[1], duplicate_rm = FALSE, verbose = verbose)
  }else{ confounders_tag <- "rando1234" }
  if(!any(resolutions %in% colnames(metadata))){
    resolutions <- filters_columns(metadata,
      include = resolutions[1], duplicate_rm = FALSE, verbose = verbose)
  }
  confounders <- unique(confounders[confounders %in% colnames(metadata)])
  tvar <- sapply(confounders, function(x) length(table(metadata[, x])) )
  confounders <- confounders[tvar > 1]
  resolutions <- unique(resolutions[resolutions %in% colnames(metadata)])
  qc_names <- qc_names[qc_names %in% colnames(metadata)]
  reductions <- ident_dimensions(colnames(metadata), "dim_")
  reductions <- reductions[grep("umap", names(reductions))]

  if(verbose) cat("  Quality Control metrics\n")
  if(length(qc_names)){
    dir.create(paste0(prefix, "_qc/"), showWarnings = FALSE)
    for(rdim in names(reductions)){
      fname <- paste0(prefix, "_qc/", rdim, ".pdf")
      if(is.file.finished(fname)) next
      pdf(fname, width = 12, height = 10)
      x <- try(ncounts_outsies(metadata, c(reductions[[rdim]], qc_names)))
      graphics.off()
    }
  }

  if(verbose) cat("  Per resolution report ----\n")
  for(nres in resolutions){
    if(verbose) cat("  @", nres, '\n')
    metadata$Identity <- factor(
      metadata[, nres],
      levels = gtools::mixedsort(unique(as.character(metadata[, nres]))))
    if(nlevels(metadata$Identity) == 1) if(verbose) cat("nothing to report :0\n")
    if(verbose) print(table(metadata$Identity))
    resdir <- paste0(prefix, "_", gsub("RNA_snn_res.|SCT_snn_res.", "res", nres), '/')
    dir.create(resdir, showWarnings = FALSE)

    centers <- list()
    for(rdim in names(reductions)){
      metadata$x1234 <- metadata[, reductions[[rdim]][1]];
      metadata$y1234 <- metadata[, reductions[[rdim]][2]]
      centers[[rdim]] <- metadata %>% group_by(Identity) %>%
        summarize(x = median(x = x1234), y = median(x = y1234))
    }; metadata <- metadata[, !colnames(metadata) %in% c("x1234", "y1234")]

    if(file.exists(config_file)){
      source("/home/ciro/scripts/clustering/R/subject_summary.R")
      config_str = yaml::read_yaml(config_file)
      if(!is.null(config_str$subject_summary)){
        tmp <- gsub("orig|orig\\.", "", config_str$subject_summary$subject_id)
        fname <- paste0(resdir, paste0(tmp, collapse = "_"), "_summary.csv")
        if(!file.exists(fname)){
          config_str$subject_summary$append_cols$names <- confounders
          ddf <- subject_summary(
            metadata_list = c(config_str$subject_summary$metadata_list,
                list(filt = metadata)),
            subject_id = config_str$subject_summary$subject_id,
            append_cols = config_str$subject_summary$append_cols,
            append_num = nres,
            verbose = verbose > 1
          )
          colnames(ddf) <- gsub("orig|orig\\.", "", colnames(ddf))
          if(verbose){ str(ddf); cat(fname, "\n") }
          write.csv(ddf, file = fname, quote = FALSE, row.names = FALSE)
        }
      }
    }

    if(verbose) cat("  Dimentional reductions\n")
    for(rdim in names(reductions)){
      if(verbose) cat("  -", rdim, "\n")
      fname = paste0(resdir, "dimentional_reduction_", rdim, '.pdf')
      if(is.file.finished(fname)) next
      p <- plot_grids(metadata, x = reductions[[rdim]][1],
        y = reductions[[rdim]][2], color = "Identity",
        centers = centers[[rdim]]) + NoLegend() + labs(x = "Dim 1", y = "Dim 2")
      pdf(fname); print(p); graphics.off()
      # fname <- paste0(resdir, "dimentional_reduction_", rdim, '_grid.pdf')
      # pdf(fname, height = 12, width = 12)
      # print(p + facet_wrap(~ Identity) + theme_bw())
      # graphics.off()
    }

    if(verbose) cat("  QC plots\n")
    fname <- paste0(resdir, '_qc_violins.pdf')
    if(!is.file.finished(fname) && length(qc_names) > 0){
      if(verbose) cat("   - Violin\n")
      pp <- lapply(qc_names, function(qcvar){
        y <- suppressMessages(
          simple_violin(dat = metadata, x = "Identity", y = qcvar) +
          scale_color_manual(values = v2cols(levels(metadata$Identity), gcols)) +
          scale_fill_manual(values = v2cols(levels(metadata$Identity), gcols)) +
            RotatedAxis() + NoLegend()
          )
        y
      })
      pdf(fname, width = 10, height = 7)
      void <- lapply(pp, function(x) print(x) )
      graphics.off()
    }
    fname <- list.files(path = resdir, pattern = "_qc_scatter")
    if((!all(is.file.finished(fname)) || length(fname) == 0) && length(qc_names) > 0){
      if(verbose) cat("   - Scatter\n")
      void <- try(qc_scatters(
        dat = metadata, thresholds = qc_names,
        theme_extra = facet_wrap(~ Identity),
        prefix = paste0(resdir, "_qc_scatter"),
        pdf_width = 15, pdf_height = 15, verbose = verbose
      ))
    }

    if(verbose) cat("  Proportions\n")
    for(orig in confounders){
      if(verbose) cat("   #", orig)
      df_plots <- metadata[!is.na(metadata[, orig]), ]
      df_plots <- df_plots[df_plots[, orig] != "", ]
      identity_levels <- gtools::mixedsort(unique(as.character(df_plots[, orig])))
      df_plots$var_interest <- factor(df_plots[, orig], levels = identity_levels)
      identities <- table(df_plots$var_interest) # needs the size per identity
      if(length(identities) == 1){ if(verbose) cat(" (N = 1)\n"); next }
      fname <- paste0(resdir, 'proportions_', sub("orig\\.|orig", "", orig))
      fnames <- paste0(
        fname, c(".csv", paste0("_", names(reductions), ".pdf"))
      )
      if(all(is.file.finished(fnames))){ if(verbose) cat(" done\n"); next }
      if(verbose) cat("\n")

      do_i_norm <- !grepl("ht_id", orig)
      tvar <- (min(identities) < 500) && do_i_norm
      thesecells <- if(((min(identities) / sum(identities)) < 0.1) && tvar){
        if(verbose)
          cat("     Minimum size cluster is < 10 %, sampling to a median of",
            median(identities), "per group\n")
        sample_even(annot = df_plots, cname = orig,
          maxln = -median(identities), v = FALSE)
      }else if(do_i_norm){
        if(verbose) cat("     Random sampling\n")
        sample_even(annot = df_plots, cname = orig, v = FALSE)
      }else{ rownames(df_plots) }
      thesecells <- unname(thesecells)

      p_bar <- plot_pct(
        x = df_plots, groups = c("var_interest", "Identity"),
        normalise = prop.normalise && do_i_norm, return_table = TRUE
      )
      write.csv(p_bar$table, file = fnames[1], quote = FALSE)
      p_bar <- p_bar$plot +
        guides(fill = list(ncol = make_grid(nlevels(df_plots$var_interest))[2])) +
        scale_fill_manual(values = v2cols(names(identities), gcols)) + coord_flip();
      p_pie <- plot_pct(
        x = df_plots, groups = c("Identity", "var_interest"),
        normalise = FALSE, type = "pie") +
        scale_fill_manual(values = v2cols(levels(df_plots$Identity), gcols)) +
        labs(fill = NULL) + guides(colour = guide_legend(override.aes = list(size = 6)))

      if(verbose) cat("     dimensional reductions")
      for(rdim in names(reductions)){ # Visualising with sampled data
        fname_i = paste0(fname, "_", rdim, ".pdf"); if(verbose) cat(" -", rdim)
        if(is.file.finished(fname_i)) next
        pp <- plot_grids(df_plots[thesecells, ], x = reductions[[rdim]][1],
          y = reductions[[rdim]][2], color = "Identity",
          centers = centers[[rdim]],
          facet = list("var_interest", make_grid(identities)[1]))
        pp2 <- plot_grids(df_plots[thesecells, ], x = reductions[[rdim]][1],
          y = reductions[[rdim]][2], color = "var_interest",
          centers = centers[[rdim]])
        if(nlevels(df_plots$var_interest) == 2)
          pp <- pp + theme(plot.margin = margin(0, 3, 0, 3, "cm"))
        pdf(fname_i, width = 15, height = 10)
        print(cowplot::plot_grid(pp, p_bar + NoLegend(), rel_widths = c(2, 0.8), ncol = 2))
        print(cowplot::plot_grid(pp2, p_pie + NoLegend(), rel_widths = c(1.5, 0.8), ncol = 2))
        graphics.off()
      }; cat("\n")
    }; cat("\n\n")
  }# end of resolutions
  if(verbose) cat("  ------- -------- -------- -------\n")
  return(NULL)
}
plot_grids = function(
  mdata, x, y, color, centers = NULL, facet = NULL, colours = NULL
) {
  color_i = names(table(mdata[, color]))
  p <- ggplot(mdata, aes_string(x = x, y = y, color = color)) +
    geom_point(size = 0.7) + labs(x = NULL, y = NULL, color = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    scale_color_manual(
      values = v2cols(color_i, colours), labels = color_i
    ) + theme_classic() + theme(
      panel.border = element_rect(fill = NA),
      strip.text.x = element_text(face = "bold")
    )
  if(!is.null(facet)){
    if(is.na(facet[2])) facet[2] = make_grid(table(mdata[, facet[[1]]]))[1]
    p <- p + facet_wrap(facets = paste0("~", facet[[1]]), ncol = as.numeric(facet[[2]]))
  }
  if(!is.null(centers)){
    p <- p +
      geom_point(data = centers, mapping = aes_string(x = "x", y = "y"), size = 0, alpha = 0, inherit.aes = FALSE) +
      geom_text(data = centers, mapping = aes_string(x = "x", y = "y", label = "Identity"), inherit.aes = FALSE)
  }
  return(p)
}

marker_report <- function(
  markers_df,
  top_n = 20,
  file = "_summary_stats.csv",
  return_plot = FALSE,
  verbose = FALSE
){
  if(verbose) cat("  ------- Plotting markers -------\n")
  return_list = list()
  if(is.null(markers_df$cluster_n)) markers_df$cluster_n <- markers_df$cluster
  tvar <- gtools::mixedsort(unique(as.character(markers_df$cluster_n)))
  markers_df$cluster_n <- factor(markers_df$cluster_n, levels = tvar)
  if(is.null(markers_df$Dpct))
    markers_df$Dpct <- 100 * (markers_df$pct.1 - markers_df$pct.2)
  if(is.null(markers_df$Significance)){
    markers_df$Significance <- (-log10(markers_df$p_val_adj))
    tvar <- max(markers_df$Significance[is.finite(markers_df$Significance)])
    markers_df$Significance[is.infinite(markers_df$Significance)] <- tvar
  }
  deltaSlope_skip <- is.null(markers_df$Dmean)
  trail <- paste0("_summary_stats|\\.", tools::file_ext(file))

  if(verbose && !is.null(markers_df$Dpct)){
    cat("Negative deltas:"); print(table(markers_df$cluster, markers_df$Dpct < 0))
  }

  if(!deltaSlope_skip){
    top_features <- get_top_n(
      x = markers_df, filter_pattern = "XY123",
      filter_neg = FALSE, n = top_n, verbose = verbose
    )
    deltaSlope <- ggplot(
      data = top_features,
      mapping = aes(x = Dpct, y = log2(Dmean), size = Significance, color = avg_logFC)
    ) + geom_point() + facet_wrap(~ cluster_n, scale = 'free') +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(x = "%(CLUSTER) - %(REST)", y = "mean(CL) - mean(REST)") +
      xlim(min(top_features$Dpct, na.rm = TRUE), 100) +
      scale_size(breaks = scales::pretty_breaks(n = 5)) + theme_minimal()
    if(requireNamespace("ggrepel", quietly = TRUE)){
      deltaSlope = deltaSlope + ggrepel::geom_text_repel(
        aes(label = gene), size = 2, color = 'red',
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )
    }
    return_list$deltaSlope = deltaSlope
  }

  top_features <- get_top_n(x = markers_df, n = top_n, verbose = verbose)
  top_features$gene <- as.character(gsub(".ENSG.*", "", top_features$gene))
  if(requireNamespace("dplyr", quietly = TRUE)){
    top_features = dplyr::mutate(
      .data = top_features,
      gene = tidytext::reorder_within(x = gene, by = -Dpct, cluster_n)
    )
  }
  deltas = ggplot(
    data = top_features,
    mapping = aes(x = gene, y = Dpct, size = Significance, color = avg_logFC)
  ) + geom_point() + facet_wrap(~ cluster_n, scale = 'free_x') +
    scale_size(breaks = scales::pretty_breaks(n = 5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = NULL, y = "%(CLUSTER) - %(REST)") +
    ylim(min(top_features$Dpct, na.rm = TRUE), 100) + theme_classic() +
    theme(
      axis.text.x = element_text(face = "bold", hjust = 1, angle = 45),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      panel.grid.major.x = element_line(size = 0.1, colour = "#f4f4f4"),
      panel.border = element_rect(fill = NA),
      strip.text.x = element_text(face = "bold")
    )
  if(requireNamespace("tidytext", quietly = TRUE)){
    deltas = deltas + tidytext::scale_x_reordered()
  }
  return_list$delta_percentage = deltas
  if(requireNamespace("UpSetR", quietly = TRUE)){
    updf <- as.data.frame.matrix(table(markers_df[, c("gene", "cluster_n")]))
    updf <- updf[, gtools::mixedsort(colnames(updf))]
    return_list$markers_overlap = UpSetR::upset(
      data = updf, sets.x.label = "Number of markers")
  }

  if(!return_plot){
    for(i in names(return_list)){
      fname_i <- paste0(gsub(trail, "", file), "_", i, ".pdf")
      if(verbose) cat(fname_i, "\n")
      pdf(fname_i, height = 12, width = 14, onefile = FALSE);
      print(return_list[[i]]); graphics.off()
    }
  }

  if(verbose) cat("  ------- -------- ------- -------\n")
  if(return_plot) return(return_list) else invisible(x = NULL)
}

variable_features_report <- function(
  object,
  prefix = NULL,
  cutoff = NULL,
  top_n = 30,
  verbose = TRUE
){
  if(verbose) cat("  ------- Plotting variable features -------\n")
  plot_and_df <- feature_variance_pct(
    object = object,
    cutoff = cutoff,
    verbose = verbose
  )
  pdf(paste0(prefix, "_variable_features_percentage.pdf")); print(plot_and_df$plot); graphics.off()
  write.csv(plot_and_df$table, file = paste0(prefix, "_variable_features_percentage.csv"))

  # Identify the 10 most highly variable features
  top10 <- head(rownames(plot_and_df$table)[plot_and_df$table$variable], top_n)#head(VariableFeatures(object), top_n)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlotp(
    object = object,
    variable_features = rownames(plot_and_df$table)[plot_and_df$table$variable]
  )
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  pdf(paste0(prefix, "_variable_features_top", top_n, ".pdf"), width = 10, height = 7);
  print(plot2); graphics.off()

  if(verbose) cat("  ------- -------- -------- -------- -------\n")
}

#' Generate cumulative percentage of highly variable genes
#'
#' This function creates a plot with the cumulative percentage of  the variance/
#' dispersion and the number of genes in order to visualise a plateau
#'
#' @param object Seurat object or data.frame with HVG information (need a
#' logical column 'variable' specifying the if a feature is a highly variable)
#' @param smethod Method used, 'mvp', 'vst', etc.
#' @param var_column Column with variance values.
#' @param var_order Column to order genes with, it takes 'var_column.'
#' @param cutoff Mean threshold used.
#' @param v Verbose
#'
#' @return Returns a plot with the variance explained by the HVGs and generates
#'
#' @importFrom
#'
#' @export
#'
#' @examples
#' plot_and_df <- feature_variance_pct(object = pbmc)
#'

feature_variance_pct <- function(
  object,
  smethod = NULL,
  var_column = NULL,
  var_order = NULL,
  cutoff = NULL,
  frac_cutoff = NULL,
  verbose = TRUE
){
  if(class(object) == "Seurat"){
    if(is.null(smethod)){
      smethod <- unique(sub("\\..*", "", names(object@assays$RNA@meta.features)))
    }
    hdf <- HVFInfo(object, selection.method = smethod)
    hdf$variable <- rownames(hdf) %in% VariableFeatures(object)
  }else if(class(object) == "seurat"){
    hdf <- object@hvg.info; hdf$variable <- rownames(hdf) %in% object@var.genes
    smethod <- "mvp"; # var_column <- "gene.dispersion"; var_order <- "gene.dispersion.scaled"
  }else{
    hdf <- object
  }; smethod <- paste0(smethod, collapse = ", ")
  if(verbose) cat("Total number genes:", nrow(hdf), "\n")
  if(verbose) cat("Methods:", ifelse(is.null(smethod), "none", smethod), "\n")
  colnames(hdf) <- sub(".*mean.*", "mean", colnames(hdf)) # re-naming mean column
  if(is.null(cutoff)) cutoff <- min(hdf$mean, na.rm = TRUE)
  # selecting the dispersion/variance column
  if(is.null(var_column)) var_column <- tail(grep("dispersion|variance", colnames(hdf)), 1)
  if(is.numeric(var_column)) var_column <- colnames(hdf)[var_column]
  if(isTRUE(grepl("scaled", var_column))){
    hdf[, paste0("shifted.", var_column)] <- hdf[, var_column] + abs(min(hdf[, var_column]))
    var_column <- paste0("shifted.", var_column)
  }
  if(is.null(var_order)) var_order <- var_column # column for ordering
  if(verbose) cat("Dispersion/variance column:", var_column, "\n")
  if(verbose) cat("Ordering by:", var_order, "\n") # ' is added to avoid Excel problems
  hdf <- cbind(gene_name = paste0("'", rownames(hdf)), hdf)

  hdf <- hdft <- hdf[order(-hdf[, var_order]), ] # ordering by dispersion/variance
  if(verbose) print(head(hdf[hdf$variable, -1, drop = FALSE], 20)); #print(summary(hdf[, ])) }
  # % of variance out of the TOTAL variance
  hvg_pct_total <- round(sum(hdf[hdf$variable, var_column]) / sum(hdf[, var_column], na.rm = TRUE) * 100, 2)
  hdf$cumulative_total = cumsum(hdf[, var_column]) # cumulative (and %) across ALL genes
  hdf$cumulative_total_pct = round(hdf$cumulative_total / sum(hdf[, var_column], na.rm = TRUE) * 100, 2)
  hdf$N_total = 1:nrow(hdf); hdf$lmean <- log2(hdf$mean + 1); hdf$lmean[!hdf$variable] <- NA
  if(!is.null(frac_cutoff) && class(object) == "Seurat"){
    hdf$exprFrac <- stats_summary_table(
      mat = object@assays$RNA@counts, rnames = rownames(hdf), moments = "p", v = verbose
    )
    hdf <- hdf[hdf$exprFrac > (frac_cutoff), ] # FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }else if(!is.null(frac_cutoff)){
    hdf <- hdf[hdf[, names(frac_cutoff)] > (frac_cutoff), ]
  }
  hdf <- hdf[hdf$mean > cutoff, ] # FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  hdf$cumulative = cumsum(hdf[, var_column]) # cumulative (and %) with genes avobe a threshold of mean expression
  hdf$cumulative_pct = round(hdf$cumulative / sum(hdf[, var_column], na.rm = TRUE) * 100, 2)
  hdf$N = 1:nrow(hdf); hdft <- cbind(hdft, hdf[rownames(hdft), !colnames(hdf) %in% colnames(hdft)])
  # write.table(hdf, file = 'hvGenes_pcts.csv', sep = ',', row.names = FALSE)
  hvg_pct <- round(sum(hdf[hdf$variable, var_column]) / sum(hdf[, var_column], na.rm = TRUE) * 100, 2)
  hvg_label <- paste0(
    "From all ", nrow(hdft), " genes: ",  hvg_pct_total, "% (", sum(hdft$variable),
    " HVGs)\nFrom ", nrow(hdf), " with mean > ", cutoff, ": ", hvg_pct, "% (",
    sum(hdf$variable), " HVGs)\nMethod: ", smethod
  )
  if(any(head(!hdf$variable, sum(hdf$variable)))){
    capt <- paste0("Warning: the 'top' genes may not be selected based on'", var_column,
    "'\nSeurat v2 dispersions tends to have this effect; or\nYou could've applied another filter before selecting features.\n")
  }else{ capt <- NULL }
  if(verbose) cat(hvg_label, "\n")
  coulsn = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
  hvg_pct <- max(hdf$cumulative_pct[hdf$variable]) # last HGV genes
  ngenes <- nrow(hdf[hdf$cumulative_pct <= hvg_pct, ]) # max($N)
  p <- ggplot(hdf, aes(x = N, y = cumulative_pct, color = lmean)) + geom_point() +
    labs(
      title = "Cumulative variance of Highly Variable Genes",
      x = paste("Number of genes, order:", var_order),
      y = paste("Cumulative percentage:", gsub("\\.", " ", var_column)),
      subtitle = hvg_label, color = "Mean", caption = capt
    ) + ylim(c(0, 100)) +
    scale_colour_gradientn(colours = coulsn, limits = c(0, max(hdf$lmean))) + theme_classic() +
    geom_hline(aes(yintercept = hvg_pct), linetype = "dashed", color = "red", size = 1.3) +
    geom_text(
      x = max(hdf$N) / 2, y = hvg_pct - 2.5,
      label = paste(ngenes, "genes = ", hvg_pct, "%"),
      size = 7, color = "black"
    )
  return(list(plot = p, table = hdft))
}


# mean-variance from HVFInfo
VariableFeaturePlotp <- function(
  object,
  cols = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000'),
  pt.size = 1.5,
  log = NULL,
  axes = NULL,
  selection.method = NULL,
  variable_features = NULL,
  assay = NULL
) {
  hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
  variable_features <- if(is.null(variable_features)){
    unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1
  }else if(is.character(variable_features)){
    unlist(rownames(hvf.info) %in% variable_features) + 1
  }else{ variable_features }
  var.status <- c("no", "yes")[variable_features]
  hvf.info <- hvf.info[, c(1, 3)]
  hvf.info$pct <- round(Matrix::rowMeans(GetAssayData(object, assay = assay)[rownames(hvf.info), ] > 0) * 100, 1)
  if(is.null(axes)){
    log <- !is.null(log) || (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    axes <- c('mean', 'variance.standardized', 'pct')
  }else{ log <- FALSE }
  axis.labels <- sapply(axes, switch,
    variance.standardized = "Standardized Variance",
    mean = "Average Expression",
    dispersion.scaled = "Dispersion",
    gmean = "Geometric Mean of Expression",
    pct = "Percentage expressing",
    "log10(mean)" = "Log10(Mean)",
    residual_variance = "Residual Variance")
  if(is.null(axis.labels[1])) axis.labels[1] <- axes[1]
  axis.labels[3] <- gsub(" ", "\n", axis.labels[3])
  hvf.info$tmp <- hvf.info[, axes[3]]
  hvf.info[var.status == "no", "tmp"] <- NA
  if(grepl("pct", axes[3])) lms <- c(0, 100) else lms <- range(hvf.info[, axes[3]])
  plot <- ggplot(data = hvf.info, aes_string(x = axes[1], y = axes[2], color = "tmp")) +
    geom_point(size = pt.size) +
    scale_color_gradientn(name = axis.labels[3], colours = cols, na.value = 'gray70', limits = lms, guide = "colorbar") +
    theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(subtitle = paste(c("Non-variable", "Variable"), "count:", table(var.status), collapse = "\n"),
      x = axis.labels[1], y = axis.labels[2])
  if (log) {
      plot <- plot + scale_x_log10()
  }
  return(plot)
}
