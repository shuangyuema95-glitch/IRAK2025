#' @title plotConditionStack
#' @description Compute cell-type proportions within each condition and draw a single stacked barplot
#'              (x-axis: conditions, y-axis: within-condition proportion, fill: cell type).
#' @param seurat_obj input Seurat object with metadata columns for condition and cell type
#' @param condition_col input Name of condition column in metadata (default "condition")
#' @param celltype_col input Name of cell-type column in metadata (default "celltype")
#' @param conditions Character vector of conditions to include (default: all found in metadata)
#' @param colors Optional color vector or named vector for cell types. If NULL a default palette is used.
#' @param exclude_celltypes Optional vector of cell-type names to exclude
#' @param order_celltypes Whether to order cell types by overall abundance (default TRUE)
#' @param show_percent Whether to draw percent labels inside stacked bars (default TRUE)
#' @param label_threshold Minimal proportion to show label (default 0.02 = 2%)
#' @param title Optional plot title
#' @return A list with elements: $plot (ggplot object) and $data (data.frame used for plotting)


library(Seurat)
library(dplyr)
library(ggplot2)
library(scales) 
library(rlang)

plotConditionStack <- function(seurat_obj,
                               condition_col = "condition",
                               celltype_col = "celltype",
                               conditions = NULL,
                               colors = NULL,
                               exclude_celltypes = NULL,
                               order_celltypes = TRUE,
                               show_percent = TRUE,
                               label_threshold = 0.02,
                               title = NULL) {
  # metadata
  meta <- as.data.frame(seurat_obj@meta.data)
  
  # choose conditions if not specified
  if (is.null(conditions)) {
    conditions <- unique(as.character(meta[[condition_col]]))
  }
  
  # filter metadata to chosen conditions
  meta <- meta[meta[[condition_col]] %in% conditions, , drop = FALSE]
  
  # optionally exclude certain cell types
  if (!is.null(exclude_celltypes)) {
    meta <- meta[!meta[[celltype_col]] %in% exclude_celltypes, , drop = FALSE]
  }
  
  # compute counts per (condition, cellType)
  counts <- meta %>%
    group_by(condition = .data[[condition_col]], cellType = .data[[celltype_col]]) %>%
    summarise(n = n(), .groups = "drop")
  
  # compute within-condition proportions
  counts <- counts %>%
    group_by(condition) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  # order cell types by overall abundance if requested
  if (order_celltypes) {
    tot <- counts %>%
      group_by(cellType) %>%
      summarise(total = sum(n), .groups = "drop") %>%
      arrange(desc(total))
    cell_order <- tot$cellType
    counts$cellType <- factor(counts$cellType, levels = cell_order)
  } else {
    counts$cellType <- factor(counts$cellType)
  }
  
  # prepare fill scale depending on provided colors
  fill_scale <- NULL
  ct_levels <- levels(counts$cellType)
  if (!is.null(colors)) {
    # if named vector and matches cell types, use direct mapping
    if (!is.null(names(colors)) && all(ct_levels %in% names(colors))) {
      mapped <- colors[ct_levels]
      fill_scale <- scale_fill_manual(values = mapped)
    } else if (length(colors) >= length(ct_levels)) {
      mapped <- colors[seq_len(length(ct_levels))]
      names(mapped) <- ct_levels
      fill_scale <- scale_fill_manual(values = mapped)
    } else {
      # fallback
      fill_scale <- scale_fill_brewer(palette = "Paired")
      warning("Provided 'colors' too short or not named; using default palette.")
    }
  } else {
    fill_scale <- scale_fill_brewer(palette = "Paired")
  }
  
  # base plot
  p <- ggplot(counts, aes(x = condition, y = prop, fill = cellType)) +
    geom_col(position = "stack", width = 0.7) +
    fill_scale +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Proportion", title = title) +
    theme_classic() +
    theme(
      legend.key.height = unit(0.3, 'cm'),
      legend.key.width  = unit(0.3, 'cm'),
      legend.title = element_blank(),
      axis.text.x = element_text(colour = "black", size = 10),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5)
    )
  
  # add percent labels (only for segments >= label_threshold to avoid clutter)
  if (isTRUE(show_percent)) {
    label_df <- counts %>%
      mutate(label = paste0(round(prop * 100, 1), "%")) %>%
      filter(prop >= label_threshold)
    
    if (nrow(label_df) > 0) {
      p <- p + geom_text(data = label_df,
                         aes(label = label),
                         position = position_stack(vjust = 0.5),
                         size = 3,
                         colour = "black")
    }
  }
  
  return(list(plot = p, data = counts))
}
