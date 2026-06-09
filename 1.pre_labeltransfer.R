library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


anno_prov <- readRDS('/home/sangho/2020_newbuild/projects/BRL_cancers/crc_gtfv4.2/scRNA-seq/2.seurat/anno_provided/tmp/crc.anno_provided.Rds') 
head(anno_prov@meta.data); nrow(anno_prov@meta.data) # 82674
anno_miss <- readRDS('/home/sangho/2020_newbuild/projects/BRL_cancers/crc_gtfv4.2/scRNA-seq/2.seurat/anno_missed/tmp/crc.anno_missing.Rds')
head(anno_miss@meta.data, n=3); nrow(anno_miss@meta.data) # 18010

anno_prov@meta.data$Sample <- NULL
anno_prov@meta.data$orig.ident <- 'anno_provided'
anno_prov@meta.data <- anno_prov@meta.data[, c(1:6, 9, 7, 8)]
anno_prov@meta.data <- droplevels(anno_prov@meta.data)
head(anno_prov@meta.data, n=3) # 82674

anno_miss@meta.data$orig.ident <- 'anno_missing'
anno_miss@meta.data$Patient <- unlist(lapply(anno_miss@meta.data$Library, function (x) unlist(strsplit(as.character(x), split = '-'))[1] ))
anno_miss@meta.data$Patient <- factor(anno_miss@meta.data$Patient, levels = c('KUL24', 'KUL27', 'KUL29'))
anno_miss@meta.data$Class <- unlist(lapply(anno_miss@meta.data$Library, function (x) paste0(unlist(strsplit(as.character(x), split = '-'))[-1], collapse = '_') ))
anno_miss@meta.data$Class <- mapvalues(anno_miss@meta.data$Class, from = unique(anno_miss@meta.data$Class), to = c('Normal', 'Border', 'Tumor'))
anno_miss@meta.data$Class <- factor(anno_miss@meta.data$Class, levels = c('Normal', 'Border', 'Tumor'))
anno_miss@meta.data <- anno_miss@meta.data[, c(1:4,6,7,5)]
head(anno_miss@meta.data) # 18026

nrow(anno_prov@meta.data) # 82674
nrow(anno_miss@meta.data) # 18010

head(anno_prov@meta.data, n=3)
head(anno_miss@meta.data, n=3) 

### Transfer ###
dir.create('tmp')
transfer.anchors <- FindTransferAnchors(reference = anno_prov, query = anno_miss, dims = 1:30)
saveRDS(transfer.anchors, 'tmp/transfer.anchors.Rds')
predictions_subtype <- TransferData(anchorset = transfer.anchors, refdata = anno_prov@meta.data$Cell_subtype, dims = 1:30)
saveRDS(predictions_subtype, 'tmp/predictions_subtype.Rds')
head(predictions_subtype)
nrow(subset(predictions_subtype, prediction.score.max >= 0.5)) # 15913

anno_miss@meta.data$Cell_subtype <- as.character(predictions_subtype$predicted.id)
anno_miss@meta.data$predictionscore <- predictions_subtype$prediction.score.max


### clean-up and matcehd colnames
anno_miss <- subset(anno_miss, predictionscore >= 0.5)
anno_miss@meta.data$predictionscore <- NULL
nrow(anno_miss@meta.data) # 15913


c("Stromal cells", "T cells", "T cells", "Myeloids", "B cells", "B cells", "T cells", 
  "Stromal cells", "Malignant cells", "Stromal cells", "Malignant cells", "Mast cells", 
  "Myeloids", "B cells", "Myeloids", "Stromal cells", "Myeloids", "Malignant cells", 
  "Epithelial cells", "Stromal cells", "Stromal cells", "Stromal cells", "Stromal cells", "Myeloids", 
  "T cells", "T cells", "T cells", "Stromal cells", "Myeloids", "Stromal cells", "Stromal cells", "Epithelial cells", 
  "Epithelial cells",  "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells")

anno_miss@meta.data$Cell_type <- mapvalues(anno_miss@meta.data$Cell_subtype, from = unique(anno_miss@meta.data$Cell_subtype),
                                           to = c("Stromal cells", "T cells", "T cells", "Myeloids", "B cells", "B cells", "T cells", 
                                                  "Stromal cells", "Malignant cells", "Stromal cells", "Malignant cells", "Mast cells", 
                                                  "Myeloids", "B cells", "Myeloids", "Stromal cells", "Myeloids", "Malignant cells", 
                                                  "Epithelial cells", "Stromal cells", "Stromal cells", "Stromal cells", "Stromal cells", "Myeloids", 
                                                  "T cells", "T cells", "T cells", "Stromal cells", "Myeloids", "Stromal cells", "Stromal cells", "Epithelial cells", 
                                                  "Epithelial cells",  "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells"))

anno_miss@meta.data$Cell_type <- factor(anno_miss@meta.data$Cell_type, levels = c("Epithelial cells", "Malignant cells", 'Myeloids', "Mast cells", "B cells", 'T cells', 'Stromal cells'))
anno_miss@meta.data <- anno_miss@meta.data[, c(1:7,9,8)]

removecells <- rownames(subset(anno_miss@meta.data, Class == "Normal" & Cell_type == "Malignant cells"))
length(removecells) # malignant cells in normal tissue; 159
anno_miss <- subset(anno_miss, cells = setdiff(rownames(anno_miss@meta.data), removecells))
anno_miss@meta.data <- droplevels(anno_miss@meta.data)
nrow(anno_miss@meta.data) # 15754

saveRDS(anno_miss, '/home/sangho/2020_newbuild/projects/BRL_cancers/crc_gtfv4.2/scRNA-seq/2.seurat/anno_missed/tmp/crc.anno_missing.pred.Rds')

anno_prov@meta.data <- droplevels(anno_prov@meta.data)
anno_prov@meta.data$Cell_type <- as.character(anno_prov@meta.data$Cell_type)
anno_prov@meta.data[rownames(subset(anno_prov@meta.data, Cell_subtype %in% c('CMS1', 'CMS2', 'CMS3', 'CMS4'))), 'Cell_type'] <- 'Malignant cells'
anno_prov@meta.data$Cell_type <- factor(anno_prov@meta.data$Cell_type, levels = c("Epithelial cells", "Malignant cells", 'Myeloids', "Mast cells", "B cells", 'T cells', 'Stromal cells'))
head(anno_prov@meta.data)
nrow(anno_prov@meta.data) # 82674

saveRDS(anno_prov, '/home/sangho/2020_newbuild/projects/BRL_cancers/crc_gtfv4.2/scRNA-seq/2.seurat/anno_provided/tmp/crc.anno_provided.Rds') 



### merge datasets ###
nrow(anno_prov@meta.data) # 82674
nrow(anno_miss@meta.data) # 15754

identical(colnames(anno_prov@meta.data), colnames(anno_miss@meta.data))
merged_label <- rbind(anno_prov@meta.data, anno_miss@meta.data)# 98428

saveRDS(merged_label, '/home/sangho/2020_newbuild/projects/BRL_cancers/crc_gtfv4.2/scRNA-seq/2.seurat/merged/tmp/merged.label.Rds')


mergedObj <- merge(x = anno_prov, y = anno_miss)
head(mergedObj@meta.data, n=3); nrow(mergedObj@meta.data) # 98428

mergedObj@meta.data$orig.ident <- factor(mergedObj@meta.data$orig.ident, levels = c("anno_provided", "anno_missing"))
mergedObj@meta.data$Cell_type <- factor(mergedObj@meta.data$Cell_type, levels = levels(anno_prov@meta.data$Cell_type))
mergedObj@meta.data$Cell_subtype <- factor(mergedObj@meta.data$Cell_subtype, levels = levels(anno_prov@meta.data$Cell_subtype))
mergedObj@meta.data$Library <- factor(mergedObj@meta.data$Library, levels = sort(unique(mergedObj@meta.data$Library)))
mergedObj@meta.data$Patient <- factor(mergedObj@meta.data$Patient, levels = sort(unique(mergedObj@meta.data$Patient)))
mergedObj@meta.data$Class <- factor(mergedObj@meta.data$Class, levels = c('Normal', 'Border', 'Tumor'))
mergedObj@meta.data$Cell_subtype_v2 <- mapvalues(mergedObj@meta.data$Cell_subtype,
                                                 from = sort(as.character(unique(mergedObj@meta.data$Cell_subtype))),
                                                 to = c(sort(as.character(unique(mergedObj@meta.data$Cell_subtype)))[1:31], 'SPP1+', 'SPP1+', 
                                                        sort(as.character(unique(mergedObj@meta.data$Cell_subtype)))[34:43]) )

subtypes <- c("Epithelial cells", "Goblet cells", "Intermediate", "Tuft cells", "Stem-like/TA", # Epithelial
              "BEST4+ Enterocytes", "Mature Enterocytes", "Mature Enterocytes type 1", "Mature Enterocytes type 2", # Epithelial
              "CMS1", "CMS2", "CMS3", "CMS4", # Malignant cells
              "Proliferating", "Anti-inflammatory", "Pro-inflammatory", "SPP1+", "cDC", "Mast cells", # Myeloids
              "CD19+CD20+ B", "IgA+ Plasma", "IgG+ Plasma", "Unspecified Plasma", # B cells
              "CD4+ T cells", "T helper 17 cells", "gamma delta T cells", "T follicular helper cells", "Regulatory T cells", "CD8+ T cells", "NK cells", # T cells
              "Stromal 1", "Stromal 2", "Stromal 3", "Myofibroblasts",  # Stromal
              "Enteric glial cells", "Pericytes", "Proliferative ECs", "Smooth muscle cells", # Stromal
              "Stalk-like ECs", "Tip-like ECs", "Lymphatic ECs") # Stromal

mergedObj@meta.data$Cell_subtype_v2 <- factor(mergedObj@meta.data$Cell_subtype_v2, levels = subtypes)


mergedObj@meta.data$Cell_type_v2 <- mapvalues(mergedObj@meta.data$Cell_subtype_v2,
                                              from = levels(mergedObj@meta.data$Cell_subtype_v2),
                                              to = c("Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", # Epithelial
                                                     "Epithelial cells", "Epithelial cells", "Epithelial cells", "Epithelial cells", # Epithelial
                                                     "Malignant cells", "Malignant cells", "Malignant cells", "Malignant cells", # Malignant
                                                     "Macrophages", "Macrophages", "Macrophages", "Macrophages", "cDC", "Mast cells", # Myeloids
                                                     "B cells", "Plasma cells", "Plasma cells", "Plasma cells", # B cells
                                                     "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", # T cells
                                                     "Fibroblasts", "Fibroblasts", "Fibroblasts", "Fibroblasts",  # Stromal
                                                     "Stromal cells", "Stromal cells", "Endothelial cells", "Stromal cells", # Stromal
                                                     "Endothelial cells", "Endothelial cells", "Endothelial cells")) # Stromal

mergedObj@meta.data$Cell_type_v2 <- factor(mergedObj@meta.data$Cell_type_v2, 
                                           levels = c("Epithelial cells", "Malignant cells", "Macrophages", "cDC", "Mast cells", 
                                                      "B cells", "Plasma cells", "T cells", "Fibroblasts", "Stromal cells", "Endothelial cells"))
mergedObj@meta.data <- mergedObj@meta.data[, c(1:9, 11, 10)]
summary(mergedObj@meta.data)

mergedObj@meta.data$Condition <- mapvalues(mergedObj@meta.data$Class, from = levels(mergedObj@meta.data$Class), to = c('Normal', 'Tumor', 'Tumor'))
mergedObj@meta.data <- mergedObj@meta.data[, c(1:6, 12, 7:11)]
summary(subset(mergedObj@meta.data, Cell_type_v2 == 'Malignant cells')$Condition)
# Normal  Tumor 
# 0       17508 
head(mergedObj@meta.data); nrow(mergedObj@meta.data) # 98428

saveRDS(mergedObj, 'tmp/mergedObj.Rds')


