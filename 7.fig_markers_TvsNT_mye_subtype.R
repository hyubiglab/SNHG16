library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

myObj <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
head(myObj@meta.data); summary(myObj@meta.data)
sigs <- list.files('./degs_myesub_tvsnt/', pattern = 'markers.Cell_celltype.*txt'); sigs

sigs <- c("markers.Cell_celltype.Proliferating.T_89_vs_NT_462.txt", "markers.Cell_celltype.Pro-inflammatory.T_642_vs_NT_2997.txt", "markers.Cell_celltype.Anti-inflammatory.T_457_vs_NT_389.txt",
          "markers.Cell_celltype.SPP1+.T_255_vs_NT_3968.txt", "markers.Cell_celltype.cDC.T_262_vs_NT_597.txt")

signature_merged <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(signature_merged) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", "Symbol", "Cell_type")
signature_ncrnas_merged <- signature_merged
top_n <- 3

knownlncrnas <- read.delim('tmp/known_lncRNAs.txt', header = F)
head(knownlncrnas)
novellncrnas <- read.delim('tmp/novel_lncRNAs.txt', header = F)
head(novellncrnas)
mye_subtype_markers <- read.delim('degs_myesub_tvsnt/Mye_subtypes_markers.txt')
summary(mye_subtype_markers$cluster)
head(mye_subtype_markers)

for (sigfile in sigs) {
  
  subtype <- unlist(strsplit(sigfile, split = '\\.'))[3]
  sig <- read.delim(paste0(c('./degs_myesub_tvsnt/', sigfile), collapse = ''))
  sig$Cell_type <- subtype
  sig$gene_symbol <- unlist( lapply(sig$gene, function (x) unlist(strsplit(as.character(x), split = '-ENSG'))[1] ) )
  
  sig_ncrnas <- subset(sig, sig$gene_symbol %in% c(as.character(knownlncrnas$V2), as.character(novellncrnas$V2)) )
  sig_ncrnas <- subset(sig_ncrnas, gene %in% as.character(subset(mye_subtype_markers, cluster == subtype)$gene))
  sig_ncrnas <- sig_ncrnas[order(sig_ncrnas$cluster, decreasing = T), ]
  
  if (nrow(sig_ncrnas) != 0) {
    signature_ncrnas_merged <- rbind(signature_ncrnas_merged, sig_ncrnas)
  }
  
  sig <- subset(sig, !gene_symbol %in% c(as.character(knownlncrnas$V2), as.character(novellncrnas$V2)))
  sig <- subset(sig, gene %in% as.character(subset(mye_subtype_markers, cluster == subtype)$gene))
  sig_top_n <- sig[c(rownames(na.omit(subset(sig, cluster == 'Tumor')[1:top_n,])), rownames(na.omit(subset(sig, cluster == 'Normal')[1:top_n,]))), ]
  signature_merged <- rbind(signature_merged, sig_top_n)
}

rownames(signature_merged) <- c(1:nrow(signature_merged))
signature_merged[rownames(subset(signature_merged, cluster == 'Normal')), 'avg_logFC'] <- -signature_merged[rownames(subset(signature_merged, cluster == 'Normal')), 'avg_logFC']
head(signature_merged); nrow(signature_merged)

rownames(signature_ncrnas_merged) <- c(1:nrow(signature_ncrnas_merged))
signature_ncrnas_merged[rownames(subset(signature_ncrnas_merged, cluster == 'Normal')), 'avg_logFC'] <- -signature_ncrnas_merged[rownames(subset(signature_ncrnas_merged, cluster == 'Normal')), 'avg_logFC']
head(signature_ncrnas_merged); nrow(signature_ncrnas_merged)


mye_subtype_markers_flt <- subset(mye_subtype_markers, cluster %in% c('Proliferating', 'Pro-inflammatory', 'Anti-inflammatory', 'SPP1+', 'cDC'))
mye_subtype_markers_flt <- subset(mye_subtype_markers_flt, Symbol %in% c(as.character(knownlncrnas$V2), as.character(novellncrnas$V2)))
remove_genes <- subset(subset(mye_subtype_markers, Symbol %in% as.character(mye_subtype_markers_flt$Symbol)), 
                       !cluster %in% c("Proliferating", "Pro-inflammatory", "Anti-inflammatory", "SPP1+", "cDC"))$Symbol
remove_genes <- as.character(remove_genes); length(remove_genes)
removed_markers <- subset(mye_subtype_markers, Symbol %in% remove_genes)
signature_ncrnas_merged <- subset(signature_ncrnas_merged, !Symbol %in% remove_genes)


### remove genes called for multiple times
duplicated <- data.frame(table(signature_ncrnas_merged$gene))
duplicated <- duplicated[duplicated$Freq > 1, ]; duplicated
duplicated <- duplicated[duplicated$Freq > 100, ]
head(duplicated); sum(duplicated$Freq)

signature_ncrnas_merged <- subset(signature_ncrnas_merged, !gene %in% duplicated$Var1)
head(signature_ncrnas_merged); nrow(signature_ncrnas_merged)

unique(signature_ncrnas_merged$Cell_type)
signature_ncrnas_merged
#saveRDS(signature_ncrnas_merged, 'tmp/fig_markers_TvsNT_mye_subtype.signature_ncrnas_merged.Rds')

############
head(myObj@meta.data)
myObj@meta.data$Cell_type_site <- paste(myObj@meta.data$Cell_type_v2, myObj@meta.data$Condition)
#sort(unique(myObj@meta.data$Cell_type_site))
myObj@meta.data$Cell_type_site <- factor(myObj@meta.data$Cell_type_site,
                                            levels = c("Epithelial cells Tumor", "Malignant cells Tumor", "Epithelial cells Normal", 
                                                       "Macrophages Tumor", "Macrophages Normal", "cDC Tumor", "cDC Normal", "Mast cells Tumor", "Mast cells Normal", 
                                                       "B cells Tumor", "B cells Normal", "Plasma cells Tumor", "Plasma cells Normal", "T cells Tumor", "T cells Normal", 
                                                       "Fibroblasts Tumor", "Fibroblasts Normal", "Stromal cells Tumor", "Stromal cells Normal", "Endothelial cells Tumor", "Endothelial cells Normal"))
levels(myObj@meta.data$Cell_type_site)

myObj@meta.data$Mye_subtypes <- as.character(myObj@meta.data$Cell_type_v2)
myObj@meta.data[rownames(subset(myObj@meta.data, Cell_type == 'Myeloids')), 'Mye_subtypes'] <- as.character(paste(subset(myObj@meta.data, Cell_type == 'Myeloids')$Cell_subtype_v2,
                                                                                                                  subset(myObj@meta.data, Cell_type == 'Myeloids')$Condition))
sort(unique(myObj@meta.data$Mye_subtypes))

myObj@meta.data$Mye_subtypes <- factor(myObj@meta.data$Mye_subtypes,
                                       levels = c(levels(myObj@meta.data$Cell_type_v2)[1:2], 
                                                  "Proliferating Tumor", "Proliferating Normal", "Pro-inflammatory Tumor", "Pro-inflammatory Normal", 
                                                  "Anti-inflammatory Tumor", "Anti-inflammatory Normal", "SPP1+ Tumor", "SPP1+ Normal", "cDC Tumor", "cDC Normal",
                                                  levels(myObj@meta.data$Cell_type_v2)[5:11]))
levels(myObj@meta.data$Mye_subtypes)
Idents(myObj) <-'Mye_subtypes'

source('dotplot_expression.R')

### lncRNA ###
lnc_count <- data.frame(table(signature_ncrnas_merged$gene))
subset(lnc_count, Freq > 1)
lnc_redun <- subset(signature_ncrnas_merged, gene %in% subset(lnc_count, Freq > 1)$Var1)
lnc_redun <- lnc_redun[order(lnc_redun$gene), ]
lnc_redun

dotplot_custom(object = myObj, scale_size = 5, clstlv = rev(levels(myObj@meta.data$Mye_subtypes)), 
               features = rev( c(setdiff(unique(as.character(signature_ncrnas_merged$gene)), 
                                         as.character(subset(lnc_count, Freq > 1)$Var1)), 
                                 intersect(unique(as.character(signature_ncrnas_merged$gene)), 
                                           as.character(subset(lnc_count, Freq > 1)$Var1))) ) )
#ggsave('markers/lncRNA.Myeloids.cellsubtype.MAST.TvsNT.pdf', units = 'cm', width = 12, height = 10)


### PCG ###
pcg_count <- data.frame(table(signature_merged$gene))
subset(pcg_count, Freq > 1)
pcg_redun <- subset(signature_merged, gene %in% subset(pcg_count, Freq > 1)$Var1)
pcg_redun <- pcg_redun[order(pcg_redun$gene), ]
pcg_redun

dotplot_custom(object = myObj, scale_size = 5, clstlv = rev(levels(myObj@meta.data$Mye_subtypes)), 
               features = rev( c(setdiff(unique(as.character(signature_merged$gene)), 
                                         as.character(subset(pcg_count, Freq > 1)$Var1)), 
                                 intersect(unique(as.character(signature_merged$gene)), 
                                           as.character(subset(pcg_count, Freq > 1)$Var1))) ) )
#ggsave('markers/pcg.Myeloids.cellsubtype.MAST.TvsNT.pdf', units = 'cm', width = 18, height = 10)


