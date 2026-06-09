save.image('~/Programs/R/workspace/BRL_FigS_UPF1_target_analysis.RData')
load('~/Programs/R/workspace/BRL_FigS_UPF1_target_analysis.RData')
library(DESeq2)
library(ggplot2)
library(pheatmap)


DLD1_siCont.cm <- read.table('~/Projects/BRL/SNHG16/UPF1_target_analysis/250617_siUPF1_RNA_seq/A13/3.FeatureCount/A13_output.txt', 
                             sep = '\t', header = T)
DLD1_siUPF1.cm <- read.table('~/Projects/BRL/SNHG16/UPF1_target_analysis/250617_siUPF1_RNA_seq/A14/3.FeatureCounts/A14_output.txt', 
                             sep = '\t', header = T)
SW480_siCont.cm <- read.table('~/Projects/BRL/SNHG16/UPF1_target_analysis/250617_siUPF1_RNA_seq/A15/3.FeatureCounts/A15_output.txt', 
                              sep = '\t', header = T)
SW480_siUPF1.cm <- read.table('~/Projects/BRL/SNHG16/UPF1_target_analysis/250617_siUPF1_RNA_seq/A16/3.FeatureCounts/A16_output.txt', 
                              sep = '\t', header = T)

rownames(DLD1_siCont.cm) <- DLD1_siCont.cm$Geneid
colnames(DLD1_siCont.cm) <- c('Gene_name', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count')
rownames(DLD1_siUPF1.cm) <- DLD1_siUPF1.cm$Geneid
colnames(DLD1_siUPF1.cm) <- c('Gene_name', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count')

rownames(SW480_siCont.cm) <- SW480_siCont.cm$Geneid
colnames(SW480_siCont.cm) <- c('Gene_name', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count')
rownames(SW480_siUPF1.cm) <- SW480_siUPF1.cm$Geneid
colnames(SW480_siUPF1.cm) <- c('Gene_name', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count')

Rseq.cm <- data.frame(DLD1_siCont.cm$Gene_name, DLD1_siCont.cm$Count, DLD1_siUPF1.cm$Count,
                      SW480_siCont.cm$Count, SW480_siUPF1.cm$Count)
rownames(Rseq.cm) <- Rseq.cm$Gene_name
colnames(Rseq.cm) <- c('Gene_name', 'DLD1_siCont', 'DLD1_siUPF1', 'SW480_siCont', 'SW480_siUPF1')
Rseq.cm <- Rseq.cm[ , -1]

Rseq2.cm <- data.frame(DLD1_siCont.cm$Gene_name, DLD1_siUPF1.cm$Count, DLD1_siCont.cm$Count,
                       SW480_siUPF1.cm$Count, SW480_siCont.cm$Count)
rownames(Rseq2.cm) <- Rseq2.cm$Gene_name
colnames(Rseq2.cm) <- c('Gene_name', 'DLD1_siUPF1', 'DLD1_siCont', 'SW480_siUPF1', 'SW480_siCont')
Rseq2.cm <- Rseq2.cm[ , -1]

coldata <- data.frame(sample = colnames(Rseq.cm),
                      condition = c('siCont', 'siUPF1', 'siCont', 'siUPF1'))
rownames(coldata) <- colnames(Rseq.cm)
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = Rseq.cm,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "siCont", "siUPF1"))
summary(res)
summary(dds)
result_05 <- results(dds, alpha = 0.05, contrast=c("condition", 'siCont', 'siUPF1'))
summary(result_05)
test.deseq <- data.frame(result_05)
test.deseq$Gene_name <- rownames(test.deseq)

SNHGs.deseq <- test.deseq
SNHGs.deseq <- SNHGs.deseq[grepl('^SNHG', SNHGs.deseq$gene_name), ]
SNHGs.deseq$gene_name <- rownames(SNHGs.deseq)

bar.df <- data.frame(SNHGs.deseq$log2FoldChange, SNHGs.deseq$pvalue, SNHGs.deseq$gene_name)
colnames(bar.df) <- c('log2FoldChange', 'pvalue', 'gene_name')

genes <- rownames(bar.df)
coldata <- colnames(bar.df)
for(g in genes){
  for(c in coldata){
    if(is.na(bar.df[g, c])){
      bar.df[g, c] <- 0
    }
    if(bar.df[g, 'pvalue'] == 0){
      bar.df[g, 'pvalue'] <- 1
    }
  }
}

gene_order <- c(
  "SNHG25", "SNHG6", "SNHG3", "SNHG29", "SNHG8", "SNHG5", "SNHG32", "SNHG9", "SNHG1",
  "SNHG19", "SNHG17", "SNHG20", "SNHG21", "SNHG15", "SNHG4", "SNHG12", "SNHG16", "SNHG10",
  "SNHG7", "SNHG14", "SNHG26", "SNHG18", "SNHG30", "SNHG11", "SNHG31", "SNHG22", "SNHG27", "SNHG28"
)

bar.df <- bar.df[order(match(bar.df$gene_name, gene_order)), ]
bar.df$gene_name <- factor(bar.df$gene_name, levels = rev(gene_order))
rownames(bar.df) <- bar.df$gene_name

coldata <- data.frame(sample = colnames(Rseq2.cm),
                      condition = c('siUPF1', 'siCont', 'siUPF1', 'siCont'))
rownames(coldata) <- colnames(Rseq2.cm)
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = Rseq2.cm,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "siUPF1", "siCont"))
summary(res)
summary(dds)
result_05 <- results(dds, alpha = 0.05, contrast=c("condition", 'siUPF1', 'siCont'))
summary(result_05)
test.deseq <- data.frame(result_05)
test.deseq$Gene_name <- rownames(test.deseq)

SNHGs.deseq <- test.deseq
SNHGs.deseq <- SNHGs.deseq[grepl('^SNHG', SNHGs.deseq$Gene_name), ]

bar.df <- data.frame(SNHGs.deseq$log2FoldChange, SNHGs.deseq$pvalue, SNHGs.deseq$Gene_name)
colnames(bar.df) <- c('log2FoldChange', 'pvalue', 'gene_name')

genes <- rownames(bar.df)
coldata <- colnames(bar.df)
for(g in genes){
  for(c in coldata){
    if(is.na(bar.df[g, c])){
      bar.df[g, c] <- 0
    }
    if(bar.df[g, 'pvalue'] == 0){
      bar.df[g, 'pvalue'] <- 1
    }
  }
}



# Figure visualization
upf1_2.bar <- ggplot(bar.df, aes(x = log2FoldChange, y = gene_name)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
upf1_2.bar


# For 단측 검정
# 양측 검정 결과 얻기 (siCont vs siUPF1)
res <- results(dds, contrast = c("condition", "siUPF1", "siCont"))

# Upregulated gene만 대상으로 단측 p-value 계산
# p-value를 절반으로 나누고, log2FoldChange가 양수일 때만 유의하게 간주
res$one_sided_pvalue <- ifelse(res$log2FoldChange > 0, res$pvalue / 2, 1)

res <- data.frame(res)
test.deseq <- data.frame(res)
test.deseq$Gene_name <- rownames(test.deseq)
SNHGs.deseq <- test.deseq
SNHGs.deseq <- SNHGs.deseq[grepl('^SNHG', SNHGs.deseq$Gene_name), ]


# Visualization
bar.df <- data.frame()
bar.df <- data.frame(SNHGs.deseq$log2FoldChange, SNHGs.deseq$one_sided_pvalue, SNHGs.deseq$Gene_name)
colnames(bar.df) <- c('log2FoldChange', 'pvalue', 'gene_name')
rownames(bar.df) <- bar.df$gene_name

genes <- rownames(bar.df)
coldata <- colnames(bar.df)
for (g in genes) {
  for (c in coldata) {
    # NA를 0으로 대체
    if (is.na(bar.df[g, c])) {
      bar.df[g, c] <- 0
    }
  }
  
  # pvalue가 0인 경우 1로 대체 (NA가 아님을 먼저 확인)
  if (!is.na(bar.df[g, "pvalue"]) && bar.df[g, "pvalue"] == 0) {
    bar.df[g, "pvalue"] <- 1
  }
}

gene_order <- c(
  "SNHG25", "SNHG6", "SNHG3", "SNHG29", "SNHG8", "SNHG5", "SNHG32", "SNHG9", "SNHG1",
  "SNHG19", "SNHG17", "SNHG20", "SNHG21", "SNHG15", "SNHG4", "SNHG12", "SNHG16", "SNHG10",
  "SNHG7", "SNHG14", "SNHG26", "SNHG18", "SNHG30", "SNHG11", "SNHG31", "SNHG22", "SNHG27", "SNHG28"
)

bar.df <- bar.df[order(match(bar.df$gene_name, gene_order)), ]
bar.df$gene_name <- factor(bar.df$gene_name, levels = rev(gene_order))
rownames(bar.df) <- bar.df$gene_name

bar.df$label <- ""
bar.df$label[bar.df$pvalue <= 0.05] <- "*"
bar.df$label[bar.df$pvalue <= 0.01] <- "**"
bar.df$label[bar.df$pvalue <= 0.001] <- "***"
bar.df$label[bar.df$pvalue <= 0.0001] <- "****"

# bar plot
new_upf1.bar <- ggplot(bar.df, aes(x = log2FoldChange, y = gene_name)) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = label, 
                x = ifelse(log2FoldChange > 0, log2FoldChange + 0.1, log2FoldChange - 0.1)),
            size = 4, hjust = ifelse(bar.df$log2FoldChange > 0, 0, 1)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
new_upf1.bar

pdf(file = '~/Projects/BRL/UPF1_target_analysis/Final_figures/new_UPF1_target_barplot.pdf', 
    width = 3, height = 5)
new_upf1.bar
dev.off()


# Changing to heatmaps
# Cell lines SNHGs expression
# Count matrixes recall
head(Rseq.cm)

# DLD1 first
# Calculate TPM, log2 normalize and draw heatmap
# Bringing gene lengths for TPM
DLD1_siCont.cm
DLD1_siUPF1.cm
dld1.cont <- data.frame(DLD1_siCont.cm$Gene_name, DLD1_siCont.cm$Count, DLD1_siCont.cm$Length)
dld1.si <- data.frame(DLD1_siUPF1.cm$Gene_name, DLD1_siUPF1.cm$Count, DLD1_siUPF1.cm$Length)

colnames(dld1.cont) <- c('Gene_name', 'Count', 'Length')
rownames(dld1.cont) <- dld1.cont$Gene_name
colnames(dld1.si) <- c('Gene_name', 'Count', 'Length')
rownames(dld1.si) <- dld1.si$Gene_name

# TPM calculation
dld1.cont$len_kb <- dld1.cont$Length / 1000
dld1.cont$rpk <- dld1.cont$Count / dld1.cont$len_kb
dld1.cont$tpm <- (dld1.cont$rpk / sum(dld1.cont$rpk, na.rm = TRUE)) * 1e6
dld1.cont$log <- log2(dld1.cont$tpm + 1)

dld1.si$len_kb <- dld1.si$Length / 1000
dld1.si$rpk <- dld1.si$Count / dld1.si$len_kb
dld1.si$tpm <- (dld1.si$rpk / sum(dld1.si$rpk, na.rm = TRUE)) * 1e6
dld1.si$log <- log2(dld1.si$tpm + 1)

dld1.fc <- data.frame(dld1.cont$Gene_name, dld1.si$log - dld1.cont$log)
colnames(dld1.fc) <- c('Gene_name', 'log2FC')
rownames(dld1.fc) <- dld1.fc$Gene_name

# SNHGs subsetting
dld1_snhg.fc <- dld1.fc[grepl('^SNHG', dld1.fc$Gene_name), ]
gene_order <- c(
  "SNHG25", "SNHG6", "SNHG3", "SNHG29", "SNHG8", "SNHG5", "SNHG32", "SNHG9", "SNHG1",
  "SNHG19", "SNHG17", "SNHG20", "SNHG21", "SNHG15", "SNHG4", "SNHG12", "SNHG16", "SNHG10",
  "SNHG7", "SNHG14", "SNHG26", "SNHG18", "SNHG30", "SNHG11", "SNHG31", "SNHG22", "SNHG27", "SNHG28"
)

# Visualization
dld1_snhg.fc$Gene_name <- factor(dld1_snhg.fc$Gene_name, levels = gene_order)
dld1_snhg.fc <- dld1_snhg.fc[order(dld1_snhg.fc$Gene_name), ]
plot.df <- dld1_snhg.fc[, "log2FC", drop = FALSE]
rownames(plot.df) <- dld1_snhg.fc$Gene_name
max_val <- max(abs(plot.df), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 100)

dld1_snhgs.heat <- pheatmap(
  plot.df,
  breaks = breaks,
  color = colorRampPalette(c("blue", "white", "red"))(length(breaks)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = 'none',
  cellwidth = 15, 
  cellheight = 15, 
  legend_breaks = seq(round(-max_val), round(max_val), by = 1)
)


# Next SW480
# Calculate TPM, log2 normalize and draw heatmap
# Bringing gene lengths for TPM
SW480_siCont.cm
SW480_siUPF1.cm
sw480.cont <- data.frame(SW480_siCont.cm$Gene_name, SW480_siCont.cm$Count, SW480_siCont.cm$Length)
sw480.si <- data.frame(SW480_siUPF1.cm$Gene_name, SW480_siUPF1.cm$Count, SW480_siUPF1.cm$Length)

colnames(sw480.cont) <- c('Gene_name', 'Count', 'Length')
rownames(sw480.cont) <- sw480.cont$Gene_name
colnames(sw480.si) <- c('Gene_name', 'Count', 'Length')
rownames(sw480.si) <- sw480.si$Gene_name

# TPM calculation
sw480.cont$len_kb <- sw480.cont$Length / 1000
sw480.cont$rpk <- sw480.cont$Count / sw480.cont$len_kb
sw480.cont$tpm <- (sw480.cont$rpk / sum(sw480.cont$rpk, na.rm = TRUE)) * 1e6
sw480.cont$log <- log2(sw480.cont$tpm + 1)

sw480.si$len_kb <- sw480.si$Length / 1000
sw480.si$rpk <- sw480.si$Count / sw480.si$len_kb
sw480.si$tpm <- (sw480.si$rpk / sum(sw480.si$rpk, na.rm = TRUE)) * 1e6
sw480.si$log <- log2(sw480.si$tpm + 1)

sw480.fc <- data.frame(sw480.cont$Gene_name, sw480.si$log - sw480.cont$log)
colnames(sw480.fc) <- c('Gene_name', 'log2FC')
rownames(sw480.fc) <- sw480.fc$Gene_name

# SNHGs subsetting
sw480_snhg.fc <- sw480.fc[grepl('^SNHG', sw480.fc$Gene_name), ]
gene_order <- c(
  "SNHG25", "SNHG6", "SNHG3", "SNHG29", "SNHG8", "SNHG5", "SNHG32", "SNHG9", "SNHG1",
  "SNHG19", "SNHG17", "SNHG20", "SNHG21", "SNHG15", "SNHG4", "SNHG12", "SNHG16", "SNHG10",
  "SNHG7", "SNHG14", "SNHG26", "SNHG18", "SNHG30", "SNHG11", "SNHG31", "SNHG22", "SNHG27", "SNHG28"
)

# Visualization
sw480_snhg.fc$Gene_name <- factor(sw480_snhg.fc$Gene_name, levels = gene_order)
sw480_snhg.fc <- sw480_snhg.fc[order(sw480_snhg.fc$Gene_name), ]
plot.df <- sw480_snhg.fc[, "log2FC", drop = FALSE]
rownames(plot.df) <- sw480_snhg.fc$Gene_name
max_val <- max(abs(plot.df), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 100)

sw480_snhgs.heat <- pheatmap(
  plot.df,
  breaks = breaks,
  color = colorRampPalette(c("blue", "white", "red"))(length(breaks)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = 'none',
  cellwidth = 15, 
  cellheight = 15, 
  legend_breaks = seq(round(-max_val), round(max_val), by = 1)
)

print(dld1_snhgs.heat)
print(sw480_snhgs.heat)

# Combining the plots
combined_fc <- merge(dld1_snhg.fc, sw480_snhg.fc, by = "Gene_name")
colnames(combined_fc) <- c('Gene_name', 'DLD1', 'SW480')
combined_plot.df <- combined_fc[, -1]
rownames(combined_plot.df) <- combined_fc$Gene_name
combined_plot.df <- combined_plot.df[gene_order, ]
max_val <- max(abs(combined_plot.df), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 100)

pdf('~/Projects/BRL/SNHG16/UPF1_target_analysis/Final_figures/new_UPF1_target_heatmap.pdf', 
    width = 3, height = 10)
pheatmap(
  combined_plot.df,
  breaks = breaks,
  color = colorRampPalette(c("blue", "white", "red"))(length(breaks)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = 'none',
  cellwidth = 15, 
  cellheight = 15,
  legend_breaks = seq(round(-max_val), round(max_val), by = 1)
)
dev.off()


###### Updating GAS5 and ZFAS1 ######
target_genes <- c("GAS5", "ZFAS1")

dld1_snhg.fc <- dld1.fc[
  grepl('^SNHG', dld1.fc$Gene_name) | dld1.fc$Gene_name %in% target_genes,
]

sw480_snhg.fc <- sw480.fc[
  grepl('^SNHG', sw480.fc$Gene_name) | sw480.fc$Gene_name %in% target_genes,
]


gene_order <- c(
  "SNHG25", "SNHG6", "SNHG3", "SNHG29", "SNHG8", "SNHG5", "SNHG32", "SNHG9", "SNHG1",
  "SNHG19", "SNHG17", "SNHG20", "SNHG21", "SNHG15", "SNHG4", "SNHG12", "SNHG16", "SNHG10",
  "SNHG7", "SNHG14", "SNHG26", "SNHG18", "SNHG30", "SNHG11", "SNHG31", "SNHG22", "SNHG27", "SNHG28",
  "GAS5", "ZFAS1"
)


# DLD1 ordering
dld1_snhg.fc$Gene_name <- factor(dld1_snhg.fc$Gene_name, levels = gene_order)
dld1_snhg.fc <- dld1_snhg.fc[order(dld1_snhg.fc$Gene_name), ]

# SW480 ordering
sw480_snhg.fc$Gene_name <- factor(sw480_snhg.fc$Gene_name, levels = gene_order)
sw480_snhg.fc <- sw480_snhg.fc[order(sw480_snhg.fc$Gene_name), ]

# Combine
combined_fc <- merge(dld1_snhg.fc, sw480_snhg.fc, by = "Gene_name")
colnames(combined_fc) <- c("Gene_name", "DLD1", "SW480")

combined_plot.df <- combined_fc[, -1]
rownames(combined_plot.df) <- combined_fc$Gene_name

combined_plot.df <- combined_plot.df[gene_order[gene_order %in% rownames(combined_plot.df)], , drop = FALSE]

max_val <- max(abs(combined_plot.df), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 100)


pdf('~/Projects/BRL/SNHG16/UPF1_target_analysis/Final_figures/new_UPF1_target_heatmap_with_GAS5.pdf', 
    width = 3, height = 10)
pheatmap(
  combined_plot.df,
  breaks = breaks,
  color = colorRampPalette(c("blue", "white", "red"))(length(breaks)),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  cellwidth = 15,
  cellheight = 15,
  legend_breaks = seq(round(-max_val), round(max_val), by = 1)
)
dev.off()














