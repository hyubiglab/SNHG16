save.image('~/Programs/R/workspace/SNHGs_expression_with_latest_genome.RData')
load('~/Programs/R/workspace/SNHGs_expression_with_latest_genome.RData')

# Reading in the latest FLAIR tpm matrix with CCD841 and DLD1
flair_all_genes.tpm <- read.table(file = '/home/minwook/Projects/BRL/SnoRNA_host/With_BIGtranscriptome_final/New_aligned_LR_seq/Quantification/FLAIR/second_trial_with_SW480/flair.quantify.tpm.tsv',
                                  header = T, sep = '\t')
colnames(flair_all_genes.tpm) <- c('Isoform_id', 'CCD841', 'DLD1', 'SW480', 'MDM_no_use')
rownames(flair_all_genes.tpm) <- flair_all_genes.tpm$Isoform_id
flair_all_genes.tpm <- data.frame(flair_all_genes.tpm)
head(flair_all_genes.tpm)


# Extracting all SNHGs isoform ids from ref_genome
BIG_transcriptome.gtf <- data.frame(read.table('~/Ref_genome/Human/hg19/BIG_transcriptome_new_ver_250107_updated/Gene_length_fixed_FINAL_USE_THIS/GENCODEv34_addliftover_BIGv3_filt.sorted.v2-2.trxname.marked.SNHGs_LUCAT1_updated.gtf_final.gtf', 
                                           header = F, sep = '\t'))
gtf_colname.list <- c('Chr', 'Source', 'Type', 'Start', 'End', 'Strandedness', 'Strand', 'Extra', 'Info')
colnames(BIG_transcriptome.gtf) <- gtf_colname.list
head(BIG_transcriptome.gtf)

target_gene_names <- c("GAS5", "ZFAS1")

SNHG_genes <- BIG_transcriptome.gtf %>%
  dplyr::filter(Type == "transcript") %>%
  dplyr::mutate(
    gene_name = sub(".*gene_name ([^;]+);.*", "\\1", Info),
    gene_id   = sub(".*gene_id ([^;]+);.*", "\\1", Info)
  ) %>%
  dplyr::filter(grepl("^SNHG", gene_name) | gene_name %in% target_gene_names) %>%
  dplyr::select(gene_name, gene_id) %>% 
  dplyr::distinct()

SNHG_transcripts <- BIG_transcriptome.gtf %>%
  dplyr::filter(Type == "transcript") %>%
  dplyr::mutate(gene_id = sub(".*gene_id ([^;]+);.*", "\\1", Info),
         transcript_id = sub(".*transcript_id ([^;]+);.*", "\\1", Info)) %>% 
  dplyr::filter(gene_id %in% SNHG_genes$gene_id) %>%
  dplyr::select(gene_id, transcript_id) %>%
  dplyr::distinct()


# Filtering genes from TPM matrix
head(flair_all_genes.tpm)
SNHGs_tpm <- subset(flair_all_genes.tpm, subset = Isoform_id %in% SNHG_transcripts$transcript_id)
head(SNHGs_tpm)

rownames(SNHG_genes) <- SNHG_genes$gene_id
rownames(SNHG_transcripts) <- SNHG_transcripts$transcript_id

Isoform_list <- rownames(SNHGs_tpm)
for(id in Isoform_list){
  temp_gene_id <- SNHG_transcripts[id, 'gene_id']
  temp_gene_name <- SNHG_genes[temp_gene_id, 'gene_name']
  SNHGs_tpm[id, 'Gene_name'] <- temp_gene_name
  rm(temp_gene_name, temp_gene_id)
}


# Subsetting by each SNHGs for total TPM calculation -> Calculating sum of TPMs
# Automated
gene_names <- unique(SNHGs_tpm$Gene_name)
processed_data <- list()

for(gene in gene_names){
  gene_tpm <- subset(SNHGs_tpm, subset = Gene_name == gene)
  gene_tpm <- gene_tpm[ , -c(1, ncol(gene_tpm))]
  gene_tpm['Total', ] <- colSums(gene_tpm)
  processed_data[[gene]] <- gene_tpm
}

for(gene in gene_names){
  assign(paste0(gene, '.tpm'), processed_data[[gene]])
}
# head(SNHG1.tpm)
# tail(SNHG1.tpm)


# Bringing all the 'Total's and scaling for heatmap.
Total_tpms <- data.frame()
for(gene in gene_names){
  Total_tpms[gene, 'DLD1'] <- processed_data[[gene]]['Total', 'DLD1']
  Total_tpms[gene, 'SW480'] <- processed_data[[gene]]['Total', 'SW480']
}
# Removing no need genes, manually.
genes_to_force_include <- c("SNHG27", "SNHG28", "GAS5", "ZFAS1")

for (g in genes_to_force_include) {
  if (!g %in% rownames(Total_tpms)) {
    Total_tpms[g, ] <- 0
  }
}

Total_tpms$Gene_name <- rownames(Total_tpms)

# Keep SNHGs + GAS5 + ZFAS1
Total_tpms <- Total_tpms[
  grepl("^SNHG", Total_tpms$Gene_name) | Total_tpms$Gene_name %in% c("GAS5", "ZFAS1"),
]

Total_tpms <- Total_tpms[, c("DLD1", "SW480")]

Total_tpms_scaled <- data.frame(scale(Total_tpms))
Total_tpms_scaled <- Total_tpms_scaled[order(Total_tpms_scaled$DLD1, decreasing = T), ]

gene_order <- rownames(Total_tpms_scaled)


# With TPM values over 100, change to 100 for better visualization
gene_list <- rownames(Total_tpms)
for(gene in gene_list){
  if(Total_tpms[gene, 'DLD1'] > 100){
    Total_tpms[gene, 'DLD1'] <- 100
  }
  if(Total_tpms[gene, 'SW480'] > 100){
    Total_tpms[gene, 'SW480'] <- 100
  }
}
Total_tpms$Gene_name <- rownames(Total_tpms)
Total_tpms$Gene_name <- rownames(Total_tpms)

Total_tpms <- Total_tpms[
  order(match(Total_tpms$Gene_name, gene_order)),
]

Total_tpms <- Total_tpms[, c("DLD1", "SW480")]

# Add GAS5 and ZFAS1 to the bottom if they were not included in gene_order
extra_order <- c("GAS5", "ZFAS1")
final_order <- c(gene_order, extra_order)
final_order <- final_order[final_order %in% rownames(Total_tpms)]

Total_tpms <- Total_tpms[final_order, , drop = FALSE]

Total_tpms_scaled <- data.frame(scale(Total_tpms))

# Visualization
library(pheatmap)
breaks <- seq(2, -2, length.out = 200)
SNHGs_expression_heatmap <- pheatmap(Total_tpms_scaled, breaks = breaks,
                                     color = colorRampPalette(c("red", "white", "blue"))(length(breaks)),
         cluster_rows = F, cluster_cols = F,
         scale = 'none', 
         cellwidth = 10,
         cellheight = 10) #Format changed with illustrator.

SNHGs_expression_heatmap_tpm <- pheatmap(Total_tpms,
                                     cluster_rows = F, cluster_cols = F,
                                     scale = 'none', 
                                     cellwidth = 10,
                                     cellheight = 10)

pdf(file = '~/Projects/BRL/SnoRNA_host/With_BIGtranscriptome_final/New_aligned_LR_seq/Final_figures/SNHGs_heatmap.pdf', 
    width = 10, height = 10)
SNHGs_expression_heatmap
SNHGs_expression_heatmap_tpm
dev.off()


# Changing the color coding
breaks <- seq(100, 0, length.out = 10000)
SNHGs_expression_heatmap_tpm_changed <- pheatmap(Total_tpms, breaks = breaks,
                                         color = colorRampPalette(c('white', 'orange', 'red'))(length(breaks)),
                                         cluster_rows = F, cluster_cols = F,
                                         scale = 'none', 
                                         cellwidth = 10,
                                         cellheight = 10)
pdf(file = '~/Projects/BRL/SnoRNA_host/With_BIGtranscriptome_final/New_aligned_LR_seq/Final_figures/SNHGs_heatmap_with_SW480.pdf', 
    width = 10, height = 10)
SNHGs_expression_heatmap_tpm_changed
dev.off()

















