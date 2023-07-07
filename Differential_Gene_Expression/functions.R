library(DESeq2)
library(ggplot2)

removeComma <- function(s) {as.numeric(gsub(",", "", s, fixed = TRUE))} 

pal_gm <- colorRampPalette(c("magenta", "black", "green"))(100)

pal_extremes <- c(rgb(173, 41, 182, maxColorValue = 255), # magenta
                  rgb(86, 194, 57, maxColorValue = 255))  # green

pal_volcano <- c(pal_extremes, "gray80")
names(pal_volcano) <- c("Significantly down-regulated",
                        "Significantly up-regulated",
                        "Not significant")

# ggplot theme used throughout  
theme_pub <- function(base_size = 11, base_family = "") {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.title = element_text(colour = "black", size = rel(1.2)),
      legend.title = element_text(colour = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = rel(0.7)),
      axis.ticks = element_line(colour = "black", size = rel(0.7))
    )
}
## PCA plot of a DEseq2 object with samples labeled  
plotPCAWithSampleNames <- function(x, intgroup, ntop=500){
  
  
  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  names = colnames(x)
  
  # factor of groups
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  plot_data <- as.data.frame(pca$x) %>% select(PC1, PC2) %>% mutate(fac = fac) %>% rownames_to_column("samples")
  
  # plot
  p <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = fac), alpha = 0.5, size=4) +
    geom_text_repel(aes(label = samples), 
                    max.overlaps = 20,
                    box.padding   = 0.5,
                    point.padding = 0.5,
                    segment.alpha = 0.25,
                    #nudge_y       = 20,
                    segment.size  = 0.2,
                    #segment.color = "grey50",
                    #size          = 4 #, direction     = "x"
    ) + 
    scale_color_discrete(name = paste(intgroup, collapse =" : ")) +
    xlab(paste("PC1 (", variance[1], "%)", sep="")) +
    ylab(paste("PC2 (", variance[2], "%)", sep=""))
  return(p)
  
}

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


# Volcano Plot from the multi-comparison Deseq2 results   
results_volcano <- function(full_lrt_df, cond1, cond2, gois){
  lfc_column <- paste("log2FoldChange", cond1, cond2, sep="_")
  padj_column <- paste("padj", cond1, cond2, sep="_")
  res_df <-  full_lrt_df %>% select(gene, all_of(!!lfc_column), all_of(!!padj_column)) %>% rename(log2FoldChange = lfc_column, padj = !!padj_column)
  
  res_df <- res_df %>% 
    dplyr::mutate(log_padj = -log10(padj)) %>% 
    dplyr::mutate(group = case_when(
      (padj < 0.05 & log2FoldChange <= 0) ~ "Significantly down-regulated",
      (padj < 0.05 & log2FoldChange >= 0) ~ "Significantly up-regulated",
      TRUE ~ "Not significant")
      ) %>% 
    dplyr::mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>% 
    dplyr::arrange(log_padj)
  p1 <- ggplot(data = res_df, aes(x = log2FoldChange, y = log_padj)) +
    geom_point(aes(colour = group), alpha = 0.8) +
    scale_colour_manual(values = pal_volcano) +
    theme_pub() +
    xlab("log2(fold change)") +
    ylab("-log10(adjusted p value)") +
    geom_text_repel(data = res_df %>% dplyr::filter(gene %in% gois), aes(label = gene),
                       box.padding   = 1.0,
                       #nudge_x       = 0.3,
                       force         = 100,
                       point.padding = 0.5,
                       #nudge_y       = 20,
                       segment.size  = 0.2,
                       segment.color = "grey50",
                       segment.alpha = 0.5,
                       size          = 4 #, direction     = "x"
  )
  return(p1)
}

# Fetch Genes Associated with a GO term.  
getGenes <- function(go){
  suppressMessages( genes <- AnnotationDbi::get(go, org.Mm.egGO2ALLEGS))
  suppressMessages( symbols <- AnnotationDbi::select( org.Mm.eg.db, keys=genes, keytype='ENTREZID', columns=c('SYMBOL')))
  geneList <- list(unique(as.character(symbols$SYMBOL))) 
  return(geneList)
}

# modification of DESeq2 plotCounts  
plotCounts_gg <- function(dds, gene, normalized = FALSE){
  plotCounts(dds, gene=gene,
             normalized = FALSE,
             returnData = TRUE,
             intgroup=c("culture", "cell_type"),
             transform = FALSE) %>% 
    ggplot(aes(x = culture, y = count)) +
    geom_point(aes(colour = cell_type)) + 
    ggtitle(gene)
  
}


# Heatmap of top genes of a PCA   
pcaHeatmap <- function(vsd, num_genes, prin_comp, ntop=500){
  
  rv = rowVars(assay(vsd))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(vsd)[select,]))
  
  pca_genes <- as.data.frame(pca[["rotation"]]) %>%
    select(.dot=prin_comp) %>%
    dplyr::rename("PC"=1) %>%
    dplyr::arrange(-abs(PC)) %>%
    dplyr::slice(1:num_genes) 
  
  mat <- assay(vsd)[row.names(pca_genes),]
  mat <- mat - rowMeans(mat)
  
  
  pheatmap(mat,
           col = pal_gm,
           scale = "row",
           #cluster_rows=FALSE,
           #clustering_distance_rows = "euclidean",
           #clustering_distance_cols = "euclidean",
           #annotation_row = selected %>% select(cluster),
           annotation_names_col = FALSE,
           #annotation_colors = list(cell_type = pal_cell_type),
           #annotation_legend = TRUE,
           fontsize = 8,
           border_color = NA,
           #cutree_rows = 4,
           cluster_cols = TRUE) 
}

resultsHeatmap <- function(results,vsd, meta_data, num_degs){
  if(num_degs > nrow(results)){
    num_degs = nrow(results)
  }
  selected <- results[1:num_degs,]
  mat <- assay(vsd)[ rownames(selected), ]
  mat <- mat - rowMeans(mat)
  annotation_col <- meta_data %>% select(culture)
  rownames(annotation_col) <- row.names(meta_data)
  pheatmap(mat,
           col = PurpleAndYellow(),
           trace = "none",
           fontsize_row = 7,
           annotation_col = annotation_col, 
           scale = "row",
           border_color = NA,
           width = 4,
           height = 4
  ) 
}

valHeatmap <- function(gene_list, vsd, meta_data){
  colfunc <- colorRampPalette(c("magenta", "black", "yellow"))
  mat <- vsd[ gene_list, ]
  mat <- mat - rowMeans(mat)
  annotation_col <- meta_data %>% dplyr::select(culture)
  rownames(annotation_col) <- row.names(meta_data)
  p <- pheatmap(mat,
                col = colfunc(100),
                trace = "none",
                fontsize_row = 7,
                annotation_col = annotation_col, 
                scale = "row",
                border_color = NA,
                width = 4,
                height = 4
  ) 
  return(p)
}

valCounts <- function(gene_list){
  val_norm <- norm_counts[gene_list,] %>% rownames_to_column("gene")
  val_long <- melt(val_norm, variable.name="sample", value.name="norm_counts") %>% left_join(meta_data %>%
                                                                                               rownames_to_column("sample") %>%
                                                                                               select(sample, culture), by = "sample")
  ggplot(val_long, aes(x = culture, y = norm_counts)) +  geom_point() +
    facet_wrap(~gene, scales = "free")
}

# function to fetch Entrez ID from gene symbol
entrz <- function(x){
  entz_list <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[["ENTREZID"]]
  return(entz_list)
}