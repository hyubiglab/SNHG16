save.image('~/Programs/R/workspace/SNHGs.RData')
# TPM already calculated
all_isoform_tpm.tsv <- data.frame(read.table('~/Projects/BRL/SnoRNA_host/SNHGs_quantification/Total/quantified_.tpm.tsv', 
                                               header = TRUE, sep = '\t'))
head(all_isoform_tpm.tsv)
rownames(all_isoform_tpm.tsv) <- all_isoform_tpm.tsv$ids
colname_list <- c('Transcript_id', 'CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                  'SRR12004858', 'SRR12004859')
colnames(all_isoform_tpm.tsv) <- colname_list

# Read GTF
SNHG_included.gtf <- data.frame(read.table('~/Projects/BRL/SnoRNA_host/SNHGs_assembled/GENCODEv34_addliftover_BIGv3_filt.sorted.v2-2.trxname.marked.SNHG16_nanopore_SNHG_updated.gtf', 
                                           header = F, sep = '\t'))
gtf_colname.list <- c('Chr', 'Source', 'Type', 'Start', 'End', 'Strandedness', 'Strand', 'Extra', 'Info')
colnames(SNHG_included.gtf) <- gtf_colname.list


# SNHGs list up
SNHG_genes <- SNHG_included.gtf %>%
  filter(grepl("SNHG", Info)) %>%
  mutate(gene_name = sub(".*gene_name ([^;]+);.*", "\\1", Info),
         gene_id = sub(".*gene_id ([^;]+);.*", "\\1", Info)) %>%
  select(gene_name, gene_id) %>% 
  distinct()

# Transcript id list up
SNHG_transcripts <- SNHG_included.gtf %>%
  filter(Type == "transcript") %>%
  mutate(gene_id = sub(".*gene_id ([^;]+);.*", "\\1", Info),
         transcript_id = sub(".*transcript_id ([^;]+);.*", "\\1", Info)) %>% 
  filter(gene_id %in% SNHG_genes$gene_id) %>%
  select(gene_id, transcript_id) %>%
  distinct()  

head(SNHG_transcripts)


# Filtering
head(all_isoform_tpm.tsv)
SNHGs_tpm <- subset(all_isoform_tpm.tsv, subset = Transcript_id %in% SNHG_transcripts$transcript_id)
head(SNHGs_tpm)


# Percentage calculation (SNHG1)
SNHG_genes.filtered <- subset(SNHG_genes, subset = grepl("^S", gene_name))
rownames(SNHG_genes.filtered) <- SNHG_genes.filtered$gene_name
SNHG_transcripts.filtered <- subset(SNHG_transcripts, subset = grepl("^E", gene_id))

SNHG_transcripts.filtered <- SNHG_transcripts.filtered %>%
  left_join(SNHG_genes.filtered, by = "gene_id")

SNHG1 <- subset(SNHG_transcripts.filtered, subset = gene_name == 'SNHG1')
SNHG1_tpm <- subset(SNHGs_tpm, subset = Transcript_id %in% SNHG1$transcript_id)
SNHG1_tpm <- SNHG1_tpm[ , -1]
SNHG1_tpm['Total', ] <- colSums(SNHG1_tpm)

SNHG1_percentage <- data.frame()
sample_list <- c('CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                  'SRR12004858', 'SRR12004859')
trancripts <- rownames(SNHG1_tpm[1:8, ])
for(trans in trancripts){
  for(samp in sample_list){
    SNHG1_percentage[trans, samp] <- SNHG1_tpm[trans, samp] / SNHG1_tpm['Total', samp] * 100
  }
}


# Percentage calculation (SNHG15)
SNHG15 <- subset(SNHG_transcripts.filtered, subset = gene_name == 'SNHG15')
SNHG15_tpm <- subset(SNHGs_tpm, subset = Transcript_id %in% SNHG15$transcript_id)
SNHG15_tpm <- SNHG15_tpm[ , -1]
SNHG15_tpm['Total', ] <- colSums(SNHG15_tpm)

SNHG15_percentage <- data.frame()
sample_list <- c('CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                 'SRR12004858', 'SRR12004859')
trancripts <- rownames(SNHG15_tpm[1:6, ])
for(trans in trancripts){
  for(samp in sample_list){
    SNHG15_percentage[trans, samp] <- SNHG15_tpm[trans, samp] / SNHG15_tpm['Total', samp] * 100
  }
}


# Percentage calculation (SNHG17)
SNHG17 <- subset(SNHG_transcripts.filtered, subset = gene_name == 'SNHG17')
SNHG17_tpm <- subset(SNHGs_tpm, subset = Transcript_id %in% SNHG17$transcript_id)
SNHG17_tpm <- SNHG17_tpm[ , -1]
SNHG17_tpm['Total', ] <- colSums(SNHG17_tpm)

SNHG17_percentage <- data.frame()
sample_list <- c('CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                 'SRR12004858', 'SRR12004859')
trancripts <- rownames(SNHG17_tpm[1:14, ])
for(trans in trancripts){
  for(samp in sample_list){
    SNHG17_percentage[trans, samp] <- SNHG17_tpm[trans, samp] / SNHG17_tpm['Total', samp] * 100
  }
}


# Percentage calculation (SNHG21)
SNHG21 <- subset(SNHG_transcripts.filtered, subset = gene_name == 'SNHG21')
SNHG21_tpm <- subset(SNHGs_tpm, subset = Transcript_id %in% SNHG21$transcript_id)
SNHG21_tpm <- SNHG21_tpm[ , -1]
SNHG21_tpm['Total', ] <- colSums(SNHG21_tpm)

SNHG21_percentage <- data.frame()
sample_list <- c('CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                 'SRR12004858', 'SRR12004859')
trancripts <- rownames(SNHG21_tpm[1:3, ])
for(trans in trancripts){
  for(samp in sample_list){
    SNHG21_percentage[trans, samp] <- SNHG21_tpm[trans, samp] / SNHG21_tpm['Total', samp] * 100
  }
}


# Percentage calculation (SNHG8)
SNHG8 <- subset(SNHG_transcripts.filtered, subset = gene_name == 'SNHG8')
SNHG8_tpm <- subset(SNHGs_tpm, subset = Transcript_id %in% SNHG8$transcript_id)
SNHG8_tpm <- SNHG8_tpm[ , -1]
SNHG8_tpm['Total', ] <- colSums(SNHG8_tpm)

SNHG8_percentage <- data.frame()
sample_list <- c('CCD841', 'DLD1', 'SRR12004855', 'SRR12004856', 'SRR12004857',
                 'SRR12004858', 'SRR12004859')
trancripts <- rownames(SNHG8_tpm[1:6, ])
for(trans in trancripts){
  for(samp in sample_list){
    SNHG8_percentage[trans, samp] <- SNHG8_tpm[trans, samp] / SNHG8_tpm['Total', samp] * 100
  }
}




CCD841_total <- sum(SNHGs_tpm$CCD841)
DLD1_total <- sum(SNHGs_tpm$DLD1)
SRR12004855_total <- sum(SNHGs_tpm$SRR12004855)
SRR12004856_total <- sum(SNHGs_tpm$SRR12004856)
SRR12004857_total <- sum(SNHGs_tpm$SRR12004857)
SRR12004858_total <- sum(SNHGs_tpm$SRR12004858)
SRR12004859_total <- sum(SNHGs_tpm$SRR12004859)

CCD841.tpm <- data.frame()
DLD1.tpm <- data.frame()
SRR12004855.tpm <- data.frame()
SRR12004856.tpm <- data.frame()
SRR12004857.tpm <- data.frame()
SRR12004858.tpm <- data.frame()
SRR12004859.tpm <- data.frame()












for(sample in sample_list){
  LUCAT1_subset.count[sample, 'Total'] <- sum(LUCAT1_subset.count[sample, 1:12])
}

LUCAT1_portion <- data.frame()

for(sample in sample_list){
  for(isoform in isoform_list){
    LUCAT1_portion[sample, isoform] <- LUCAT1.tpm[isoform, sample] / LUCAT1.tpm['Total', sample] * 100
  }
}

write.table(LUCAT1_portion, '~/Projects/BRL/LUCAT1/Isoform_Quantification/FLAIR/isoform_percent.tsv', 
            sep = '\t', row.names = T, col.names = T)


# Calculating TPM for LUCAT1 gene
sample_list
for(sample in sample_list){
  LUCAT1.tpm['Total', sample] <- sum(LUCAT1.tpm[1:12, sample])
}

tpm.portion






















