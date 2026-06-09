library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
source('dotplot_expression.R')


myObj_flt <- readRDS('tmp/crc.mtcut.norm.harmony.Rds')
head(myObj_flt@meta.data); nrow(myObj_flt@meta.data) # 98428


#dir.create('markers')
#######
# PCG #
#######
dotplot_custom(object = myObj_flt, scale_size = 5, clstlv = rev(levels(myObj_flt@meta.data$Cell_type_v2)),
               features = c("PLVAP-ENSG00000130300.9-5", "CLDN5-ENSG00000184113.9-7", "RAMP2-ENSG00000131477.11-5",
                            "TAGLN-ENSG00000149591.17-5", "ACTA2-ENSG00000107796.13-8", "MYL9-ENSG00000101335.10-4",
                            "COL1A1-ENSG00000108821.14-5", "COL1A2-ENSG00000164692.18-4", "COL3A1-ENSG00000168542.16-5",
                            "CD8A-ENSG00000153563.15-5", "CD3E-ENSG00000198851.10-4", "CD3D-ENSG00000167286.9-3",
                            "IGHG1-ENSG00000211896.7-7", "JCHAIN-ENSG00000132465.11-6", "MZB1-ENSG00000170476.16-5", 
                            "CD79B-ENSG00000007312.13-6", "CD79A-ENSG00000105369.10-5", "MS4A1-ENSG00000156738.18-5",
                            "TPSB2-ENSG00000197253.13-5", "TPSAB1-ENSG00000172236.17-4", "CPA3-ENSG00000163751.4-4",
                            "CLEC10A-ENSG00000132514.13-5", "CD1E-ENSG00000158488.16-7", "CD1C-ENSG00000158481.13-4",
                            "SPP1-ENSG00000118785.14-6", "C1QA-ENSG00000173372.17-5", "S100A9-ENSG00000163220.11-4", 
                            "KRT8-ENSG00000170421.12-8", "KRT19-ENSG00000171345.13-6", "EPCAM-ENSG00000119888.11-4"))
#ggsave('markers/pcg.celltype_v2.pdf', units = 'cm', width = 18, height = 8)


###
novel_lncRNAs <- read.delim('tmp/novel_lncRNAs.txt', header = F)
novel_lncRNAs$V2 <- as.character(novel_lncRNAs$V2)
head(novel_lncRNAs); nrow(novel_lncRNAs) # 20971
known_lncRNAs <- read.delim('tmp/known_lncRNAs.txt', header = F)
known_lncRNAs$V2 <- as.character(known_lncRNAs$V2)
head(known_lncRNAs); nrow(known_lncRNAs) # 16944

celltype_degs <- read.delim('degs/markers.Cell_type_v2.MAST.txt'); head(celltype_degs)
celltype_degs_nc <- subset(celltype_degs, Symbol %in% c(known_lncRNAs$V2, novel_lncRNAs$V2))
head(celltype_degs_nc); nrow(celltype_degs_nc) # 135 genes


#################
# known lncRNAs #
#################
knownlnc_count <- data.frame(table(subset(celltype_degs_nc, Symbol %in% known_lncRNAs$V2)$gene))
subset(knownlnc_count, Freq > 1)
subset(celltype_degs_nc, gene %in% subset(knownlnc_count, Freq == 1)$Var1)
knownlnc_redun <- subset(celltype_degs_nc, gene %in% subset(knownlnc_count, Freq > 1)$Var1)
knownlnc_redun <- knownlnc_redun[order(knownlnc_redun$gene), ]
knownlnc_redun

dotplot_custom(object = myObj_flt, scale_size = 5, clstlv = rev(levels(myObj_flt@meta.data$Cell_type_v2)), 
               features = rev( c(setdiff(unique(subset(celltype_degs_nc, Symbol %in% known_lncRNAs$V2)$gene), 
                                         as.character(subset(knownlnc_count, Freq > 1)$Var1)), 
                                 intersect(unique(subset(celltype_degs_nc, Symbol %in% known_lncRNAs$V2)$gene), 
                                           as.character(subset(knownlnc_count, Freq > 1)$Var1))) ) )
#ggsave('markers/lncRNA.known.celltype_v2.MAST.pdf', units = 'cm', width = 24, height = 9)


#################
# novel lncRNAs #
#################
novellnc_count <- data.frame(table(subset(celltype_degs_nc, Symbol %in% novel_lncRNAs$V2)$gene))
subset(novellnc_count, Freq > 1)
subset(celltype_degs_nc, gene %in% subset(novellnc_count, Freq == 1)$Var1)
novellnc_redun <- subset(celltype_degs_nc, gene %in% subset(novellnc_count, Freq > 1)$Var1)
novellnc_redun <- novellnc_redun[order(novellnc_redun$gene), ]
novellnc_redun

dotplot_custom(object = myObj_flt, scale_size = 5, clstlv = rev(levels(myObj_flt@meta.data$Cell_type_v2)),
               features = rev( c(setdiff(unique(subset(celltype_degs_nc, Symbol %in% novel_lncRNAs$V2)$gene), 
                                         as.character(subset(novellnc_count, Freq > 1)$Var1)),
                                 intersect(unique(subset(celltype_degs_nc, Symbol %in% novel_lncRNAs$V2)$gene), 
                                           as.character(subset(novellnc_count, Freq > 1)$Var1))) ) )
#ggsave('markers/lncRNA.novel.celltype_v2.MAST.pdf', units = 'cm', width = 18, height = 9)


