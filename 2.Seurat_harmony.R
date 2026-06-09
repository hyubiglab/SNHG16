library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(harmony)

crc <- readRDS('tmp/mergedObj.Rds')
nrow(crc@meta.data) # 98428 cells
Idents(crc) <- 'Cell_type_v2'

dir.create('stats')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')
dir.create('degs')

plt <- ggplot(crc@meta.data, aes(Library, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.byLibrary.pdf', units = 'cm', width = 50, height = 8)
plt <- ggplot(crc@meta.data, aes(Library, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nFeature_RNA.byLibrary.pdf', units = 'cm', width = 50, height = 8)

plt <- ggplot(crc@meta.data, aes(Library, percent.mt)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.byLibrary.pdf', units = 'cm', width = 50, height = 8)

plt <- ggplot(crc@meta.data, aes(percent.mt, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.density_byLibrary.pdf', units = 'cm', width = 24, height = 6)
plt <- ggplot(crc@meta.data, aes(nCount_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.density_byLibrary.pdf', units = 'cm', width = 24, height = 6)
plt <- ggplot(crc@meta.data, aes(nFeature_RNA, col = Library)) + geom_density() +
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nFeature_RNA.density_byLibrary.pdf', units = 'cm', width = 24, height = 6)


plot1 <- FeatureScatter(object = crc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library')
plot2 <- FeatureScatter(object = crc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-3.stats.factorCombinations.mitoCut.byLibrary.pdf', units = 'cm', width = 50, height = 10)

plot1 <- FeatureScatter(object = crc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Cell_type_v2')
plot2 <- FeatureScatter(object = crc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Cell_type_v2')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-3.stats.factorCombinations.mitoCut.byCell_type_v2.pdf', units = 'cm', width = 30, height = 10)


######
crc.mtcut.norm <- NormalizeData(object = crc, normalization.method = "LogNormalize", scale.factor = 10000)
crc.mtcut.norm <- FindVariableFeatures(object = crc.mtcut.norm, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = crc.mtcut.norm), 10)
plot1 <- VariableFeaturePlot(object = crc.mtcut.norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2
ggsave('stats/2-3.FindVariableFeatures.mitoCut.pdf', units = 'cm', width = 17, height = 10)

crc.mtcut.norm <- ScaleData(object = crc.mtcut.norm, vars.to.regress = c('nCount_RNA', 'percent.mt'))
crc.mtcut.norm <- RunPCA(object = crc.mtcut.norm, features = VariableFeatures(crc.mtcut.norm), npcs = 120)
crc.mtcut.norm <- JackStraw(object = crc.mtcut.norm, num.replicate = 100, dims = 120)
crc.mtcut.norm <- ScoreJackStraw(object = crc.mtcut.norm, dims = 1:120)
saveRDS(crc.mtcut.norm, 'tmp/crc.mtcut.norm.Rds')

JackStrawPlot(crc.mtcut.norm, dims = 1:120) # 97 PC
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 45, height = 12)

### harmony ###
crc.mtcut.norm <- RunHarmony(crc.mtcut.norm, group.by.vars = "Library")
saveRDS(crc.mtcut.norm, 'tmp/crc.mtcut.norm.harmony.Rds')


### START HERE 11/23 ###
### UMAP ###
for (seed in c(211123150:211123169)){
  crc.mtcut.norm <- RunUMAP(crc.mtcut.norm, dims=1:97, reduction = "harmony", reduction.key='UMAP', n.components=2, min.dist=0.3, seed.use = seed)
  DimPlot(crc.mtcut.norm, reduction = 'umap', pt.size = .5, label = T, label.size = 2)
  ggsave(paste0(c('umap/compare/umap.1_2.celltype_v2.seed', seed, '.png'), collapse = ''), units = 'cm', width = 15, height = 10)
} # 211123157, 211123159

### tSNE ###
for (seed in c(211123100:211123109)){
  crc.mtcut.norm <- RunTSNE(crc.mtcut.norm, dims=1:97, reduction = "harmony", reduction.key='tSNE', dim.embed=2, seed.use = seed)
  DimPlot(crc.mtcut.norm, reduction = 'tsne', pt.size = .5, label = T, label.size = 2)
  ggsave(paste0(c('tsne/compare/tsne.1_2.celltype_v2.seed', seed, '.png'), collapse = ''), units = 'cm', width = 15, height = 10)
} # 211123105


crc.mtcut.norm <- RunUMAP(crc.mtcut.norm, dims=1:97, reduction = "harmony", reduction.key='UMAP', n.components=2, min.dist=0.3, seed.use = 211123159)
DimPlot(crc.mtcut.norm, reduction = 'umap', pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.mitoCut.1_2.celltype_v2.seed211123159.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=crc.mtcut.norm, group.by = "Library", reduction = 'umap')
ggsave('umap/umap.mitoCut.1_2.libraries.seed211123159.pdf', units = 'cm', width = 25, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'nCount_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.nCount_RNA.seed211123159.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.nFeature_RNA.seed211123159.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'percent.mt', cols = c("grey","red"), reduction = "umap")
ggsave('umap/umap.mitoCut.1_2.percent.mt.seed211123159.pdf', units = 'cm', width = 11, height = 10)


crc.mtcut.norm <- RunTSNE(crc.mtcut.norm, dims=1:97, reduction = "harmony", reduction.key='tSNE', dim.embed=2, seed.use = 211123105)
DimPlot(crc.mtcut.norm, reduction = 'tsne', pt.size = .5)
ggsave('tsne/tsne.mitoCut.1_2.celltype_v2.seed211123105.pdf', units = 'cm', width = 15, height = 10)
DimPlot(object=crc.mtcut.norm, group.by = "Library", reduction = 'tsne')
ggsave('tsne/tsne.mitoCut.1_2.libraries.seed211123105.pdf', units = 'cm', width = 25, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'nCount_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.nCount_RNA.seed211123105.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'nFeature_RNA', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.nFeature_RNA.seed211123105.pdf', units = 'cm', width = 12, height = 10)
FeaturePlot(crc.mtcut.norm, features = 'percent.mt', cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.percent.mt.seed211123105.pdf', units = 'cm', width = 11, height = 10)


head(crc.mtcut.norm@meta.data, n = 3)
writeLabel <- as.matrix(crc.mtcut.norm@meta.data)
writeLabel <- data.frame(Barcode = rownames(writeLabel), writeLabel, check.names = F)
write.table(writeLabel, 'crc.harmony.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)
saveRDS(crc.mtcut.norm, 'tmp/crc.mtcut.norm.harmony.Rds')


### DEGs ###
# cell type
Idents(crc.mtcut.norm) <- 'Cell_type_v2'
crc.mtcut.norm.markers <- FindAllMarkers(object = crc.mtcut.norm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
crc.mtcut.norm.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- crc.mtcut.norm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = crc.mtcut.norm, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/markers.Cell_type_v2.MAST.pdf', units = 'cm', width = 50, height = 40)
crc.mtcut.norm.markers$Symbol <- unlist(lapply(crc.mtcut.norm.markers$gene, function (x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] ))
write.table(crc.mtcut.norm.markers, 'degs/markers.Cell_type_v2.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 


# cell subtype
Idents(crc.mtcut.norm) <- 'Cell_subtype_v2'
crc.mtcut.norm.markers <- FindAllMarkers(object = crc.mtcut.norm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
crc.mtcut.norm.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- crc.mtcut.norm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = crc.mtcut.norm, features = top10$gene, angle = 90, size = 2, raster = T, draw.lines = F)
ggsave('degs/markers.Cell_subtype_v2.MAST.pdf', units = 'cm', width = 50, height = 40)
crc.mtcut.norm.markers$Symbol <- unlist(lapply(crc.mtcut.norm.markers$gene, function (x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] ))
write.table(crc.mtcut.norm.markers, 'degs/markers.Cell_subtype_v2.MAST.txt', sep = '\t', quote = F, col.names = T, row.names = F) 



Idents(crc.mtcut.norm) <- 'Cell_type_v2'
crc.mtcut.norm <- ScaleData(object = crc.mtcut.norm, vars.to.regress = c('nCount_RNA', 'percent.mt'), features = rownames(crc.mtcut.norm))
saveRDS(crc.mtcut.norm, 'tmp/crc.mtcut.norm.harmony_allscaled.Rds')

