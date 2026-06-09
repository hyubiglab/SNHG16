library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

myObj_flt <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
head(myObj_flt@meta.data); nrow(myObj_flt@meta.data) # 98428 cells
summary(myObj_flt@meta.data)
getPalette <- colorRampPalette(brewer.pal(9, "Spectral"))

levels(myObj_flt@meta.data$Cell_type_v2)
getPalette(11)
ggplot(myObj_flt@meta.data, aes(reorder(Cell_type_v2, -nFeature_RNA), percent.mt, fill = reorder(Cell_type_v2, -nFeature_RNA))) +
  geom_boxplot(outlier.size = .1) +
  scale_fill_manual(values = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                                  "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                                  "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                                  "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) +
  labs(x = '', y = 'Mito (%)') +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid = element_blank(),
        legend.position = 'none')
#ggsave('stats/2-5.stats.mitoContents.byCell_type_v2.pdf', units = 'cm', width = 5, height = 5)

ggplot(myObj_flt@meta.data, aes(reorder(Cell_type_v2, -nFeature_RNA), nFeature_RNA, fill = reorder(Cell_type_v2, -nFeature_RNA))) +
  geom_boxplot(outlier.size = .1) +
  scale_fill_manual(values = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                                  "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                                  "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                                  "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) +
  labs(x = '', y = 'Gene count') +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid = element_blank(),
        legend.position = 'none')
#ggsave('stats/2-5.stats.nFeature_RNA.byCell_type_v2.pdf', units = 'cm', width = 5, height = 5)

ggplot(myObj_flt@meta.data, aes(reorder(Cell_type_v2, -nFeature_RNA), nCount_RNA, fill = reorder(Cell_type_v2, -nFeature_RNA))) +
  geom_boxplot(outlier.size = .1) +
  scale_fill_manual(values = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                                  "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                                  "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                                  "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) +
  labs(x = '', y = 'UMI count') +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid = element_blank(),
        legend.position = 'none')
#ggsave('stats/2-5.stats.nCount_RNA.byCell_type_v2.pdf', units = 'cm', width = 5, height = 5)


### cluster statistics ###
Celltype_count <- data.frame(Cell_type_v2 = levels(myObj_flt@meta.data$Cell_type_v2), count = summary(myObj_flt@meta.data$Cell_type_v2))
Celltype_count$Cell_type_v2 <- factor(Celltype_count$Cell_type_v2, levels = levels(myObj_flt@meta.data$Cell_type_v2))
ggplot(Celltype_count, aes(reorder(Cell_type_v2, -count), count, fill = Cell_type_v2)) +
  geom_bar(stat = 'identity', col = 'black') +
  scale_fill_manual(values = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                                  "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                                  "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                                  "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) +
  labs(x = '', y = 'Count') +
  theme_bw(base_size = 7) +
  theme(axis.text = element_text(colour = 'black', size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid = element_blank(),
        legend.position = 'none')
#ggsave('stats/2-5.stats.cell_type_v2.pdf', units = 'cm', width = 4, height = 4)



### dimension reduction ###
DimPlot(myObj_flt, reduction = 'umap', group.by = 'Cell_type_v2', cols = getPalette(11)) + labs(x = 'UMAP 1', y = 'UMAP 2') +
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype_v2.seed211123159.pdf', units = 'cm', height = 10, width = 14)
DimPlot(myObj_flt, reduction = 'umap', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'UMAP 1', y = 'UMAP 2') + theme_void() + theme(legend.position = 'none') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype_v2.seed211123159.void.png', units = 'cm', height = 8, width = 8)

myObj_flt <- subset(myObj_flt, Class != 'Normal')
DimPlot(myObj_flt, reduction = 'umap', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'UMAP 1', y = 'UMAP 2') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype_v2.seed211123159.tumor.pdf', units = 'cm', height = 10, width = 14)
DimPlot(myObj_flt, reduction = 'umap', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'UMAP 1', y = 'UMAP 2') + theme_void() + theme(legend.position = 'none') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype_v2.seed211123159.tumor.void.png', units = 'cm', height = 8, width = 8)


myObj_flt <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
myObj_flt <- subset(myObj_flt, Class == 'Normal'); myObj_flt@meta.data <- droplevels(myObj_flt@meta.data)
DimPlot(subset(myObj_flt, Class == 'Normal'), reduction = 'umap', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'UMAP 1', y = 'UMAP 2') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype.seed211123159.normal.pdf', units = 'cm', height = 10, width = 14)
DimPlot(subset(myObj_flt, Class == 'Normal'), reduction = 'umap', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'UMAP 1', y = 'UMAP 2') + theme_void() + theme(legend.position = 'none') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'umap')[, 1]), max(Embeddings(myObj_flt, reduction = 'umap')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'umap')[, 2]), max(Embeddings(myObj_flt, reduction = 'umap')[, 2]))
#ggsave('umap/umap.mitoCut.1_2.celltype.seed211123159.normal.void.png', units = 'cm', height = 8, width = 8)


myObj_flt <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')

DimPlot(myObj_flt, reduction = 'tsne', group.by = 'Cell_type_v2', cols = getPalette(11)) + labs(x = 'tSNE 1', y = 'tSNE 2') +
  xlim(min(Embeddings(myObj_flt, reduction = 'tsne')[, 1]), max(Embeddings(myObj_flt, reduction = 'tsne')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'tsne')[, 2]), max(Embeddings(myObj_flt, reduction = 'tsne')[, 2]))
#ggsave('tsne/tsne.mitoCut.1_2.celltype_v2.seed211123105.pdf', units = 'cm', height = 10, width = 14)
DimPlot(myObj_flt, reduction = 'tsne', group.by = 'Cell_type_v2', pt.size = .5,
        cols = list("Epithelial cells" = "#D53E4F", "Malignant cells" = "#ED6345", 
                    "Macrophages" = "#F99455", "cDC" = "#FDC271", "Mast cells" = "#FEE695", 
                    "B cells" = "#FFFFBF", "Plasma cells" = "#EBF79F", "T cells" = "#C2E69F", 
                    "Fibroblasts" = "#8FD2A4", "Stromal cells" = "#5BB6A9", "Endothelial cells" = "#3288BD")) + 
  labs(x = 'tSNE 1', y = 'tSNE 2') + theme_void() + theme(legend.position = 'none') + 
  xlim(min(Embeddings(myObj_flt, reduction = 'tsne')[, 1]), max(Embeddings(myObj_flt, reduction = 'tsne')[, 1])) + 
  ylim(min(Embeddings(myObj_flt, reduction = 'tsne')[, 2]), max(Embeddings(myObj_flt, reduction = 'tsne')[, 2]))
#ggsave('tsne/tsne.mitoCut.1_2.celltype_v2.seed211123105.void.png', units = 'cm', height = 8, width = 8)



### lncRNA ###
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

myObj_flt <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
head(myObj_flt@meta.data); nrow(myObj_flt@meta.data) # 98428

novel_lncRNAs <- read.delim('tmp/novel_lncRNAs.txt', header = F)
novel_lncRNAs$V2 <- as.character(novel_lncRNAs$V2)
head(novel_lncRNAs); nrow(novel_lncRNAs) # 20971
known_lncRNAs <- read.delim('tmp/known_lncRNAs.txt', header = F)
known_lncRNAs$V2 <- as.character(known_lncRNAs$V2)
known_lncRNAs$V4 <- paste(known_lncRNAs$V2, gsub(known_lncRNAs$V1, pattern = '_', replacement = '-'), sep = '-')
head(known_lncRNAs); nrow(known_lncRNAs) # 16944
lncRNAs <- c(known_lncRNAs$V4, novel_lncRNAs$V2); length(lncRNAs) # 37915

known_pcgs <- read.delim('tmp/known_pcgs.txt', header = F)
known_pcgs$V2 <- as.character(known_pcgs$V2)
known_pcgs$V4 <- paste(known_pcgs$V2, gsub(known_pcgs$V1, pattern = '_', replacement = '-'), sep = '-')
head(known_pcgs); nrow(known_pcgs) # 20109
20109 + 37915 # 58024

#
variable_lncRNAs <- intersect(VariableFeatures(object = myObj_flt), lncRNAs)
length(variable_lncRNAs) # 266

myObj_flt_hvf <- myObj_flt@assays$RNA@meta.features
head(myObj_flt_hvf)

length(intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), known_pcgs$V4)) # 1567
length(intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), known_lncRNAs$V4)) # 119
length(intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), novel_lncRNAs$V2)) # 147

myObj_flt_hvf$class <- 'N.S.'
myObj_flt_hvf[ intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), known_pcgs$V4) , 'class'] <- 'PCGs'
myObj_flt_hvf[ intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), known_lncRNAs$V4) , 'class'] <- 'Known lncRNAs'
myObj_flt_hvf[ intersect(rownames(subset(myObj_flt_hvf, vst.variable == TRUE)), novel_lncRNAs$V2) , 'class'] <- 'Novel lncRNAs'
myObj_flt_hvf[ rownames(subset(myObj_flt_hvf, class == 'N.S.' & vst.variable == TRUE)) , 'class'] <- 'Others'
myObj_flt_hvf$class <- factor(myObj_flt_hvf$class, levels = c('Known lncRNAs', 'Novel lncRNAs', 'PCGs', 'Others', 'N.S.'))

head(myObj_flt_hvf)

ggplot(myObj_flt_hvf, aes(log10(vst.mean), vst.variance.standardized, col = class)) +
  geom_point(position = 'identity') +
  scale_color_manual(values = list('Others' = 'grey60',
                                   'PCGs' = 'steelblue2',
                                   'Known lncRNAs' = 'red2',
                                   'Novel lncRNAs' = 'red4',
                                   'N.S.' = 'grey90')) +
  labs(x = 'log10(average expression)', y = 'Standardized variance') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('stats/2-3.FindVariableFeatures.mitoCut.GeneClass.pdf', units = 'cm', width = 12, height = 8)


ggplot(subset(myObj_flt_hvf, class != 'N.S.'), aes(class, vst.variance.standardized, fill = class)) +
  geom_boxplot(outlier.size = .5) +
  scale_fill_manual(values = list('Others' = 'grey60',
                                  'PCGs' = 'steelblue2',
                                  'Known lncRNAs' = 'red2',
                                  'Novel lncRNAs' = 'red4',
                                  'N.S.' = 'grey90')) +
  labs(x = '', y = 'Standardized variance') +
  scale_y_continuous(trans='log10', breaks = c(1, 2, 5, 10, 100)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#ggsave('stats/2-3.standardized_variance_boxplot.GeneClass.pdf', units = 'cm', width = 8, height = 6)

t.test(subset(myObj_flt_hvf, class == 'PCGs')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'Known lncRNAs')$vst.variance.standardized) # p-value = 3.298e-16

t.test(subset(myObj_flt_hvf, class == 'PCGs')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'Novel lncRNAs')$vst.variance.standardized) # p-value = 6.579e-07

t.test(subset(myObj_flt_hvf, class == 'Known lncRNAs')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'Novel lncRNAs')$vst.variance.standardized) # 0.387

t.test(subset(myObj_flt_hvf, class == 'Others')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'PCGs')$vst.variance.standardized) # 0.0005047

t.test(subset(myObj_flt_hvf, class == 'Others')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'Known lncRNAs')$vst.variance.standardized) # 1.179e-07

t.test(subset(myObj_flt_hvf, class == 'Others')$vst.variance.standardized,
       subset(myObj_flt_hvf, class == 'Novel lncRNAs')$vst.variance.standardized) # 8.276e-07


###
library(reshape2)
head(myObj_flt_hvf)
summary(myObj_flt_hvf$vst.variance.standardized)
nrow(myObj_flt_hvf)/10


myObj_flt_hvf_sorted <- myObj_flt_hvf[order(myObj_flt_hvf$vst.variance.standardized, decreasing = F), ]
head(myObj_flt_hvf_sorted)

myObj_flt_hvf_sorted$class_all <- 'Others'
myObj_flt_hvf_sorted[ intersect(rownames(myObj_flt_hvf_sorted), known_pcgs$V4) , 'class_all'] <- 'PCGs'
myObj_flt_hvf_sorted[ intersect(rownames(myObj_flt_hvf_sorted), known_lncRNAs$V4) , 'class_all'] <- 'Known lncRNAs'
myObj_flt_hvf_sorted[ intersect(rownames(myObj_flt_hvf_sorted), novel_lncRNAs$V2) , 'class_all'] <- 'Novel lncRNAs'
myObj_flt_hvf_sorted$class_all <- factor(myObj_flt_hvf_sorted$class_all, levels = c('Known lncRNAs', 'Novel lncRNAs', 'PCGs', 'Others'))
summary(myObj_flt_hvf_sorted$class_all)


count_by_bin <- data.frame(matrix(nrow = 10, ncol = 4))
colnames(count_by_bin) <- c('Known lncRNAs', 'Novel lncRNAs', 'PCGs', 'Others')
head(count_by_bin)

for (i in c(1:10)) {
  start <- 6066*(i-1) +1
  end <- 6066*i
  
  tmp_bin <- myObj_flt_hvf_sorted[c(start:end), ]
  #summary(tmp_bin$class_all)
  count_by_bin[i, ] <- summary(tmp_bin$class_all)
}
count_by_bin

count_by_bin <- data.frame(Bin = rownames(count_by_bin), count_by_bin, check.names = F)
count_by_bin_m <- melt(count_by_bin)
count_by_bin_m$Bin <- factor(count_by_bin_m$Bin, levels = c(1:10))
head(count_by_bin_m)

ggplot(count_by_bin_m, aes(Bin, value, col = variable, group = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('Others' = 'grey60', 'PCGs' = 'steelblue2',
                                'Known lncRNAs' = 'red2', 'Novel lncRNAs' = 'red4'), 'Class') +
  ylim(0, max(count_by_bin_m$value)) +
  labs(x = 'Variance bin', y = 'Gene count') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('stats/2-6.standardized_variance_byBin.GeneClass.pdf', units = 'cm', width = 10, height = 4)

count_by_bin_m_prop <- count_by_bin_m
count_by_bin_m_prop[count_by_bin_m_prop$variable == 'Known lncRNAs', 'value'] <- subset(count_by_bin_m_prop, variable == 'Known lncRNAs')$value/sum(subset(count_by_bin_m_prop, variable == 'Known lncRNAs')$value)*100
count_by_bin_m_prop[count_by_bin_m_prop$variable == 'Novel lncRNAs', 'value'] <- subset(count_by_bin_m_prop, variable == 'Novel lncRNAs')$value/sum(subset(count_by_bin_m_prop, variable == 'Novel lncRNAs')$value)*100
count_by_bin_m_prop[count_by_bin_m_prop$variable == 'PCGs', 'value'] <- subset(count_by_bin_m_prop, variable == 'PCGs')$value/sum(subset(count_by_bin_m_prop, variable == 'PCGs')$value)*100
count_by_bin_m_prop[count_by_bin_m_prop$variable == 'Others', 'value'] <- subset(count_by_bin_m_prop, variable == 'Others')$value/sum(subset(count_by_bin_m_prop, variable == 'Others')$value)*100
count_by_bin_m_prop

ggplot(count_by_bin_m_prop, aes(Bin, value, col = variable, group = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('Others' = 'grey60', 'PCGs' = 'steelblue2',
                                'Known lncRNAs' = 'red2', 'Novel lncRNAs' = 'red4'), 'Class') +
  ylim(0, max(count_by_bin_m_prop$value)) +
  labs(x = 'Variance bin', y = 'Gene count (%)') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('stats/2-6.standardized_variance_byBin_prop.GeneClass.pdf', units = 'cm', width = 10, height = 4)



### Q3 above
summary(myObj_flt_hvf_sorted$vst.variance.standardized)[5]
myObj_flt_hvf_sorted_q3 <- subset(myObj_flt_hvf_sorted, vst.variance.standardized >= summary(myObj_flt_hvf_sorted$vst.variance.standardized)[5]) # Q3 = 0.9999576
nrow(myObj_flt_hvf_sorted_q3) # 15702

count_by_bin_q3 <- data.frame(matrix(nrow = 10, ncol = 4))
colnames(count_by_bin_q3) <- c('Known lncRNAs', 'Novel lncRNAs', 'PCGs', 'Others')
head(count_by_bin_q3)

for (i in c(1:10)) {
  start <- 1570*(i-1) +1
  if (i == 10) {
    end <- nrow(myObj_flt_hvf_sorted_q3)
  } else {
    end <- 1570*i
  }
  
  tmp_bin <- myObj_flt_hvf_sorted[c(start:end), ]
  #summary(tmp_bin$class_all)
  count_by_bin_q3[i, ] <- summary(tmp_bin$class_all)
}
count_by_bin_q3

count_by_bin_q3 <- data.frame(Bin = rownames(count_by_bin_q3), count_by_bin_q3, check.names = F)
count_by_bin_q3_m <- melt(count_by_bin_q3)
count_by_bin_q3_m$Bin <- factor(count_by_bin_q3_m$Bin, levels = c(1:10))
head(count_by_bin_q3_m)

ggplot(count_by_bin_q3_m, aes(Bin, value, col = variable, group = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('Others' = 'grey60', 'PCGs' = 'steelblue2',
                                'Known lncRNAs' = 'red2', 'Novel lncRNAs' = 'red4'), 'Class') +
  ylim(0, max(count_by_bin_q3_m$value)) +
  labs(x = 'Variance (>=1) bin', y = 'Gene count') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('stats/2-6.standardized_variance_byBin.q3above.GeneClass.pdf', units = 'cm', width = 10, height = 4)


count_by_bin_q3_m_prop <- count_by_bin_q3_m
count_by_bin_q3_m_prop[count_by_bin_q3_m_prop$variable == 'Known lncRNAs', 'value'] <- subset(count_by_bin_q3_m_prop, variable == 'Known lncRNAs')$value/sum(subset(count_by_bin_q3_m_prop, variable == 'Known lncRNAs')$value)*100
count_by_bin_q3_m_prop[count_by_bin_q3_m_prop$variable == 'Novel lncRNAs', 'value'] <- subset(count_by_bin_q3_m_prop, variable == 'Novel lncRNAs')$value/sum(subset(count_by_bin_q3_m_prop, variable == 'Novel lncRNAs')$value)*100
count_by_bin_q3_m_prop[count_by_bin_q3_m_prop$variable == 'PCGs', 'value'] <- subset(count_by_bin_q3_m_prop, variable == 'PCGs')$value/sum(subset(count_by_bin_q3_m_prop, variable == 'PCGs')$value)*100
count_by_bin_q3_m_prop[count_by_bin_q3_m_prop$variable == 'Others', 'value'] <- subset(count_by_bin_q3_m_prop, variable == 'Others')$value/sum(subset(count_by_bin_q3_m_prop, variable == 'Others')$value)*100
count_by_bin_q3_m_prop
ggplot(count_by_bin_q3_m_prop, aes(Bin, value, col = variable, group = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('Others' = 'grey60', 'PCGs' = 'steelblue2',
                                'Known lncRNAs' = 'red2', 'Novel lncRNAs' = 'red4'), 'Class') +
  ylim(0, max(count_by_bin_q3_m_prop$value)) +
  labs(x = 'Variance (>=1) bin', y = 'Gene count (%)') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'))
#ggsave('stats/2-6.standardized_variance_byBin.q3above_prop.GeneClass.pdf', units = 'cm', width = 10, height = 4)


