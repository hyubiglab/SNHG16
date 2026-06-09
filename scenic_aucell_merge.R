library(ggplot2)
library(plyr)

dirs <- list.files('./', pattern = 'sample')

aucell_avg <- data.frame(matrix(ncol = 0, nrow = 11))
rownames(aucell_avg) <- c("Epithelial cells" , "Malignant cells", "Macrophages", "cDC", "Mast cells", 
                          "B cells", "Plasma cells", "T cells", "Fibroblasts", "Stromal cells", "Endothelial cells")
aucell_max <- aucell_avg
regulon_counts <- aucell_avg

for (dir in dirs) {
  aucell <- read.delim(paste(c(dir, '/4.auc_10kb.txt'), collapse = ''), row.names = 1, check.names = F)
  colnames(aucell) <- unlist( lapply(colnames(aucell), function(x) unlist(strsplit(as.character(x), split = '\\('))[1] ) )
  label <- read.delim(paste(c('tables/labels/label.merged.', dir, '.txt'), collapse = ''), row.names = 1)

  if ( identical(rownames(label), rownames(aucell)) ) {
    
    for (ct in levels(label$Cell_type_v2)) {
      barcodes <- rownames(subset(label, Cell_type_v2 == ct))
      aucell_ct <- aucell[barcodes, ]
      
      for (grn in colnames(aucell_ct)) {
        if ( grn %in% colnames(aucell_avg) ) {
          if (is.na(aucell_avg[ct, grn])) {
            aucell_avg[ct, grn] <- mean(aucell_ct[, grn])
            aucell_max[ct, grn] <- mean(aucell_ct[, grn])
            regulon_counts[ct, grn] <- 1
          } else {
            aucell_avg[ct, grn] <- mean(c(aucell_avg[ct, grn], mean(aucell_ct[, grn])))
            aucell_max[ct, grn] <- max(c(aucell_avg[ct, grn], mean(aucell_ct[, grn])))
            regulon_counts[ct, grn] <- regulon_counts[ct, grn] +1
          }
        } else {
          aucell_avg[ct, grn] <- mean(aucell_ct[, grn])
          aucell_max[ct, grn] <- mean(aucell_ct[, grn])
          regulon_counts[ct, grn] <- 1
        }
      } # grn
    } # ct
    
  }
}
  
head(aucell_avg, n=3); dim(aucell_avg)
head(aucell_max, n=3); dim(aucell_max)
head(regulon_counts)

freq <- data.frame(t(regulon_counts[1, ]))
colnames(freq) <- 'Count'
head(freq)

summary(freq$Count)

ggplot(freq, aes(Count)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = c(summary(freq$Count)[3:4]), col = 'red2') +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  labs(x = 'TF frequency', y = 'Count') +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('scenic_aucell_merge.TFfreq.pdf', units = 'cm', width = 4, height = 3)

label <- read.delim('../tmp/label.merged.txt', row.names = 1)
head(label)
trunc(summary(label$Cell_subtype_v2)*0.2)


#atleast2 <- colnames(regulon_counts[1, regulon_counts[1, ] != 1]); atleast2
atleast3 <- colnames(regulon_counts[1, regulon_counts[1, ] > 2]); atleast3
atleast4 <- colnames(regulon_counts[1, regulon_counts[1, ] > 3]); atleast4
atleast5 <- colnames(regulon_counts[1, regulon_counts[1, ] > 4]); atleast5

setdiff(atleast3, atleast4)
setdiff(atleast4, atleast5)


library(pheatmap)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

use_tfs <- atleast5 ###

pheatmap(as.matrix(aucell_avg[, use_tfs]), 
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         treeheight_row = 15, treeheight_col = 15, cellwidth = 7, cellheight = 7, fontsize_row = 7, fontsize_col = 7, border_color = NA)#, filename = 'scenic_aucell_merge.avg.atleast5.pdf')
dev.off()
pheatmap(as.matrix(aucell_avg[, use_tfs]), scale = 'column',
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         treeheight_row = 15, treeheight_col = 15, cellwidth = 7, cellheight = 7, fontsize_row = 7, fontsize_col = 7, border_color = NA)#, filename = 'scenic_aucell_merge.avg.atleast5.scaled.pdf')
dev.off()

pheatmap(as.matrix(aucell_max[, use_tfs]), 
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         treeheight_row = 15, treeheight_col = 15, cellwidth = 7, cellheight = 7, fontsize_row = 7, fontsize_col = 7, border_color = NA)#, filename = 'scenic_aucell_merge.max.atleast5.pdf')
dev.off()
pheatmap(as.matrix(aucell_max[, use_tfs]), scale = 'column',
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         treeheight_row = 15, treeheight_col = 15, cellwidth = 7, cellheight = 7, fontsize_row = 7, fontsize_col = 7, border_color = NA)#, filename = 'scenic_aucell_merge.max.atleast5.scaled.pdf')
dev.off()


rep_four <- setdiff(atleast4, atleast5)
pheatmap(as.matrix(aucell_avg[, rep_four]), scale = 'column',
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
         treeheight_row = 15, treeheight_col = 15, cellwidth = 7, cellheight = 7, fontsize_row = 7, fontsize_col = 7, border_color = NA)#, filename = 'scenic_aucell_merge.avg.atleast5.scaled.pdf')
dev.off()


