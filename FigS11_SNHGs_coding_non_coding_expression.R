# save.image('~/Programs/R/workspace/SNHGs_coding_non_coding_comparison_with_latest_genome.RData')
# load('~/Programs/R/workspace/SNHGs_coding_non_coding_comparison_with_latest_genome.RData')
# Bring out coding and non-coding isoform ids (CCD841)
ccd841_gene_ids <- read.table('/lustre/export/home/hansoll/projects/polysome_fraction/analyzed/BIFUNCTIONAL/Part2/outputs/1.CCD841_gene_classes.txt',
                              sep = '\t', header = T)

# First filtering all snhgs
ccd841_snhg_ids <- ccd841_gene_ids[grepl('^SNHG', ccd841_gene_ids$gene_name), ]

# Dividing all SNHGs by the gene name
gene_names <- unique(ccd841_snhg_ids$gene_name)
processed_data <- list()

for(gene in gene_names){
  snhgs_ids <- subset(ccd841_snhg_ids, subset = gene_name == gene)
  processed_data[[gene]] <- snhgs_ids
}

for(gene in gene_names){
  assign(paste0(gene, '_isoforms'), processed_data[[gene]])
}

# Comparing the expressions by transcripts
# Dividing the coding isoforms' tpm and non-codings' tpm values
# First, fix the transcripts' names
bifunctional_gene_list <- c('SNHG1', 'SNHG15', 'SNHG21', 'SNHG17', 'SNHG16', 'SNHG8')
for(gene in bifunctional_gene_list){
  tmp_tpm_file <- get(paste0(gene, '.tpm'))
  trxs_before_change <- rownames(tmp_tpm_file)
  trxs_before_change <- trxs_before_change[trxs_before_change != 'Total']
  trxs_after_change <- sub('_.*$', '', trxs_before_change)
  rownames(tmp_tpm_file) <- c(trxs_after_change, 'Total')
  assign(paste0(gene, '.tpm'), tmp_tpm_file)
}

# Dividing and adding the tpms of coding and non-coding isoforms
for(gene in bifunctional_gene_list){
  tmp_tpm_file <- get(paste0(gene, '.tpm'))
  tmp_isoform_file <- tryCatch(get(paste0(gene, '_isoforms')), 
                               error = function(e) NULL)
  
  if(is.null(tmp_isoform_file)){
    next
  }

  tmp_tpm_file$Transcript_id <- rownames(tmp_tpm_file)
  rownames(tmp_isoform_file) <- tmp_isoform_file$transcript_id

  tmp_coding_form_file <- subset(tmp_isoform_file, subset = tmp_isoform_file$coding_type == 'sORF' | tmp_isoform_file$coding_type == 'coding')
  tmp_coding_isoforms <- rownames(tmp_coding_form_file)

  tmp_coding_tpm <- subset(tmp_tpm_file, subset = Transcript_id %in% tmp_coding_isoforms)
  tmp_non_coding_tpm <- subset(tmp_tpm_file, subset = !(Transcript_id %in% tmp_coding_isoforms) & Transcript_id != 'Total')
  
  tmp_coding_tpm <- tmp_coding_tpm[ , -4]
  tmp_non_coding_tpm <- tmp_non_coding_tpm[ , -4]
  
  tmp_coding_tpm['Total', ] <- colSums(tmp_coding_tpm)
  tmp_non_coding_tpm['Total', ] <- colSums(tmp_non_coding_tpm)
  
  assign(paste0(gene, '_coding.tpm'), tmp_coding_tpm)
  assign(paste0(gene, '_non_coding.tpm'), tmp_non_coding_tpm)
}

# Visualization; Heatmaps with grey if non exists
# SNHG16 needs transcript ids updated
Total_tpms$Gene_name <- rownames(Total_tpms)
Total_tpms <- Total_tpms[grepl('^SNHG', Total_tpms$Gene_name), ]

gene_list <- rownames(Total_tpms)

ccd841_coding_non_coding <- data.frame()
for(gene in gene_list){
  tmp_coding_tpm <- tryCatch(get(paste0(gene, '_coding.tpm')),
                             error = function(e) NULL)
  if(is.null(tmp_coding_tpm)){
    ccd841_coding_non_coding[gene, 'Coding_tpm'] <- NA
    ccd841_coding_non_coding[gene, 'Non_coding_tpm'] <- Total_tpms[gene, 'CCD841']
  } else{
    ccd841_coding_non_coding[gene, 'Coding_tpm'] <- tmp_coding_tpm['Total', 'CCD841']
    ccd841_coding_non_coding[gene, 'Non_coding_tpm'] <- Total_tpms[gene, 'CCD841'] - tmp_coding_tpm['Total', 'CCD841']
  }
}

# Manually changing SNHG16's data due to different isoform names
ccd841_coding_non_coding['SNHG16', 'Coding_tpm'] <- 18.7694571
ccd841_coding_non_coding['SNHG16', 'Non_coding_tpm'] <- 42.3034686 - 18.7694571
ccd841_coding_non_coding$Gene_name <- rownames(ccd841_coding_non_coding)

# Visualization
ccd841_coding_non_coding <- ccd841_coding_non_coding[order(match(ccd841_coding_non_coding$Gene_name, gene_order)), ]
ccd841_visualization <- data.frame(ccd841_coding_non_coding$Coding_tpm, ccd841_coding_non_coding$Non_coding_tpm)
colnames(ccd841_visualization) <- c('Coding', 'Non_coding')
gene_list <- rownames(ccd841_visualization)

# Converting to > 100 values
ccd841_visualization[1:5, 'Non_coding'] <- 100

breaks <- seq(100, 0, length.out = 10000)
ccd841_coding_non_coding_heatmap <- pheatmap(ccd841_visualization, breaks = breaks,
                                             color = colorRampPalette(c('red', 'orange', 'white'))(length(breaks)),
                                             cluster_rows = F, cluster_cols = F,
                                             scale = 'none',
                                             cellwidth = 10,
                                             cellheight = 10) #Format changed with illustrator.

pdf(file = '~/Projects/BRL/SnoRNA_host/With_BIGtranscriptome_final/New_aligned_LR_seq/Final_figures/CCD841_SNHGs_coding_noncoding_heatmap.pdf', 
    width = 10, height = 10)
ccd841_coding_non_coding_heatmap
dev.off()



# DLD1
# Bring out coding and non-coding isoform ids (CCD841)
dld1_gene_ids <- read.table('/lustre/export/home/hansoll/projects/polysome_fraction/analyzed/BIFUNCTIONAL/Part2/outputs/1.DLD1_gene_classes.txt',
                              sep = '\t', header = T)

# First filtering all snhgs
dld1_snhg_ids <- dld1_gene_ids[grepl('^SNHG', dld1_gene_ids$gene_name), ]

# Dividing all SNHGs by the gene name
gene_names <- unique(dld1_snhg_ids$gene_name)
processed_data <- list()

for(gene in gene_names){
  snhgs_ids <- subset(dld1_snhg_ids, subset = gene_name == gene)
  processed_data[[gene]] <- snhgs_ids
}

for(gene in gene_names){
  assign(paste0(gene, '_isoforms'), processed_data[[gene]])
}

# Comparing the expressions by transcripts
# Dividing the coding isoforms' tpm and non-codings' tpm values
# First, fix the transcripts' names
bifunctional_gene_list <- c('SNHG1', 'SNHG15', 'SNHG21', 'SNHG17', 'SNHG16', 'SNHG8')
for(gene in bifunctional_gene_list){
  tmp_tpm_file <- get(paste0(gene, '.tpm'))
  trxs_before_change <- rownames(tmp_tpm_file)
  trxs_before_change <- trxs_before_change[trxs_before_change != 'Total']
  trxs_after_change <- sub('_.*$', '', trxs_before_change)
  rownames(tmp_tpm_file) <- c(trxs_after_change, 'Total')
  assign(paste0(gene, '.tpm'), tmp_tpm_file)
}

# Dividing and adding the tpms of coding and non-coding isoforms
for(gene in bifunctional_gene_list){
  tmp_tpm_file <- get(paste0(gene, '.tpm'))
  tmp_isoform_file <- tryCatch(get(paste0(gene, '_isoforms')), 
                               error = function(e) NULL)
  
  if(is.null(tmp_isoform_file)){
    next
  }
  
  tmp_tpm_file$Transcript_id <- rownames(tmp_tpm_file)
  rownames(tmp_isoform_file) <- tmp_isoform_file$transcript_id
  
  tmp_coding_form_file <- subset(tmp_isoform_file, subset = tmp_isoform_file$coding_type == 'sORF' | tmp_isoform_file$coding_type == 'coding')
  tmp_coding_isoforms <- rownames(tmp_coding_form_file)
  
  tmp_coding_tpm <- subset(tmp_tpm_file, subset = Transcript_id %in% tmp_coding_isoforms)
  tmp_non_coding_tpm <- subset(tmp_tpm_file, subset = !(Transcript_id %in% tmp_coding_isoforms) & Transcript_id != 'Total')
  
  tmp_coding_tpm <- tmp_coding_tpm[ , -4]
  tmp_non_coding_tpm <- tmp_non_coding_tpm[ , -4]
  
  tmp_coding_tpm['Total', ] <- colSums(tmp_coding_tpm)
  tmp_non_coding_tpm['Total', ] <- colSums(tmp_non_coding_tpm)
  
  assign(paste0(gene, '_coding.tpm'), tmp_coding_tpm)
  assign(paste0(gene, '_non_coding.tpm'), tmp_non_coding_tpm)
}

# Visualization; Heatmaps with grey if non exists
# SNHG16 needs transcript ids updated
Total_tpms$Gene_name <- rownames(Total_tpms)
Total_tpms <- Total_tpms[grepl('^SNHG', Total_tpms$Gene_name), ]

gene_list <- rownames(Total_tpms)

dld1_coding_non_coding <- data.frame()
for(gene in gene_list){
  tmp_coding_tpm <- tryCatch(get(paste0(gene, '_coding.tpm')),
                             error = function(e) NULL)
  if(is.null(tmp_coding_tpm)){
    dld1_coding_non_coding[gene, 'Coding_tpm'] <- NA
    dld1_coding_non_coding[gene, 'Non_coding_tpm'] <- Total_tpms[gene, 'DLD1']
  } else{
    dld1_coding_non_coding[gene, 'Coding_tpm'] <- tmp_coding_tpm['Total', 'DLD1']
    dld1_coding_non_coding[gene, 'Non_coding_tpm'] <- Total_tpms[gene, 'DLD1'] - tmp_coding_tpm['Total', 'DLD1']
  }
}

# Manually changing SNHG16's data due to different isoform names
dld1_coding_non_coding['SNHG16', 'Coding_tpm'] <- 7.748066
dld1_coding_non_coding['SNHG16', 'Non_coding_tpm'] <- 15.331279 - 7.748066
dld1_coding_non_coding$Gene_name <- rownames(dld1_coding_non_coding)

# Visualization
dld1_coding_non_coding <- dld1_coding_non_coding[order(match(dld1_coding_non_coding$Gene_name, gene_order)), ]
dld1_visualization <- data.frame(dld1_coding_non_coding$Coding_tpm, dld1_coding_non_coding$Non_coding_tpm)
colnames(dld1_visualization) <- c('Coding', 'Non_coding')

dld1_visualization[1:7, 'Non_coding'] <- 100

breaks <- seq(100, 0, length.out = 10000)
dld1_coding_non_coding_heatmap <- pheatmap(dld1_visualization, breaks = breaks,
                                           color = colorRampPalette(c('red', 'orange', 'white'))(length(breaks)),
                                             cluster_rows = F, cluster_cols = F,
                                             scale = 'none', 
                                             cellwidth = 10,
                                             cellheight = 10) #Format changed with illustrator.

pdf(file = '~/Projects/BRL/SnoRNA_host/With_BIGtranscriptome_final/New_aligned_LR_seq/Final_figures/DLD1_SNHGs_coding_noncoding_heatmap.pdf', 
    width = 10, height = 10)
dld1_coding_non_coding_heatmap
dev.off()



































