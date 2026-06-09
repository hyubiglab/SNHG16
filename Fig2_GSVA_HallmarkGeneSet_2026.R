#GSVA using hallmark gene sets (JAK_STAT 추가) 

library('GSEABase')
library('GSVA')
library('pheatmap')
#library('ComplexHeatmap')

merged <- read.delim('merged.readCount.cbAdj.tpm.txt', header = TRUE) # Input tpm count matrix

sum_data <- data.frame(ID = merged$ID, Symbol = merged$Symbol)
sum_data$'DLD1(w/M0)' <- (merged$DLD1_co.rep1 + merged$DLD1_co.rep2 + merged$DLD1_co.rep3) / 3
sum_data$'HT29(w/M0)' <- (merged$HT29_co.rep2 + merged$HT29_co.rep3) / 2
sum_data$'SW480(w/M0)' <- merged$SW480_co.rep4
sum_data$'DLD1' <- (merged$DLD1.rep1 + merged$DLD1.rep2 + merged$DLD1.rep3) / 3
sum_data$'HT29' <- (merged$HT29.rep2 + merged$HT29.rep3) / 2
sum_data$'SW480' <- merged$SW480.rep4

data <- sum_data
data[1:4, 1:4]

ids <- data$ID
symbols <- data$Symbol
data$ID = NULL
data$Symbol = NULL
colnames(data)

# Changing to continuous data for gaussian approach (log2TPM)
data <- log2(data)

label <- data.frame(matrix(nrow = ncol(data), ncol = 1))
rownames(label) <- colnames(data)
colnames(label) <- 'Sample'
label$Sample <- colnames(data)
head(label)
identical(rownames(label), colnames(data))

symbols_df <- as.data.frame(symbols)

pData <- AnnotatedDataFrame(data=label)
rownames(data) <- make.unique(symbols_df$symbols)
exprSet <- ExpressionSet(as.matrix(data), phenoData=pData); exprSet

table(duplicated(symbols_df)) # FALSE

ex <- exprs(exprSet)

### Get gene sets ###
genesetGmt <- 'h.all.v2023.2.Hs.symbols.gmt'
geneSet <- read.delim(genesetGmt, header = F, row.names = 1); geneSet[,1:4]

jak_stat <- read.delim('KEGG_JAK_STAT_SIGNALING_PATHWAY.v2026.1.Hs.gmt', header = FALSE)
rownames(jak_stat) <- jak_stat$V1; jak_stat$V1 <- NULL
rownames(jak_stat) <- unlist(lapply(rownames(jak_stat), function(x) unlist(strsplit(as.character(x), split = 'KEGG_'))[2]))

rownames(geneSet) <- unlist(lapply(rownames(geneSet), function(x) unlist(strsplit(as.character(x), split = 'HALLMARK_'))[2]))
geneSet[1:4, 1:4]

## Changing 'IL6_JAK_STAT_SIGNALING' to 'JAK_STAT_SIGNALING_PATHWAY' in HALLMARK GENESET
#geneSet[24, 1:ncol(jak_stat)] <- jak_stat[1, ]
#rownames(geneSet)[24] <- 'JAK_STAT_SIGNALING'

geneSet[51, 1:ncol(jak_stat)] <- jak_stat[1,]

geneSet[51, is.na(geneSet[51, ])] <- ""

rownames(geneSet)[51] <- 'JAK_STAT_SIGNALING'

gsInfo <- data.frame(row.names = rownames(geneSet), Source = geneSet$V2); head(gsInfo)
geneSet$V2 <- NULL
geneSet <- data.frame(t(geneSet))
head(geneSet)

gsInfo <- AnnotatedDataFrame(data=gsInfo)
identical(rownames(pData(gsInfo)), colnames(geneSet))
signatureSet <- ExpressionSet(as.matrix(geneSet), phenoData=gsInfo)

head(exprs(signatureSet))

gs_list <- list()
for (i in c(1:ncol(exprs(signatureSet)))) {
  gs_list[[i]] <- GeneSet(exprs(signatureSet)[, i][exprs(signatureSet)[, i] != ''], setName = as.character(colnames(signatureSet)[i]))
}
mySignatures <- GeneSetCollection(object = gs_list) #
#saveRDS(mySignatures, 'tmp/mySignatures.Rds')
remove(gs_list)


### Calculate enrichment scores ###
### GSVA ###

#es.gsva <- gsva(expr = exprSet, gset.idx.list = mySignatures, method='gsva', verbose=T, kcdf = 'Gaussian', parallel.sz =1)

param <- gsvaParam(
  exprData = exprSet,
  geneSets = mySignatures,
  kcdf = "Gaussian"
)
es.gsva <- gsva(param)
#saveRDS(es.gsva, 'tmp/es.gsva.Rds')

gsvaScores <- as.data.frame(t(exprs(es.gsva)))
head(gsvaScores)

#write.table(gsvaScores, 'gsva/gsvaScores.txt', quote = F, col.names = NA, sep = '\t')


### ssGSEA ###

#es.ss <- gsva(exprSet, gset.idx.list = mySignatures, 
#              method="ssgsea", verbose=FALSE, kcdf = 'Gaussian', parallel.sz =1)
#es.ss <- gsva(param)

es.ss <- gsva(ssgseaParam(exprSet, mySignatures), verbose = FALSE)

#saveRDS(es.ss, 'tmp/es.ss.Rds')

ssgseaScores <- as.data.frame(t(exprs(es.ss)))
head(ssgseaScores)

#write.table(ssgseaScores, 'ssgsea/ssgseaScores.txt', quote = F, col.names = NA, sep = '\t')


callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


##rownames sorting!!
t_gsvaScores <- t(gsvaScores)
gsva
col_ann1 <- data.frame(rownames(gsvaScores))
colnames(col_ann1) <- 'Condition'
rownames(col_ann1) <- colnames(t_gsvaScores)
col_ann1
test_matrix <- data.matrix(col_ann1, rownames.force=NA)


#row_order <- c('HEME_METABOLISM', 'JAK_STAT_SIGNALING', 'COAGULATION', 'MYOGENESIS', 'KRAS_SIGNALING_DN', 'BILE_ACID_METABOLISM',
#               'NOTCH_SIGNALING', 'INFLAMMATORY_RESPONSE', 'TNFA_SIGNALING_VIA_NFKB', 'TGF_BETA_SIGNALING', 'HYPOXIA', 'EPITHELIAL_MESENCHYMAL_TRANSITION',
#               'ANGIOGENESIS', 'KRAS_SIGNALING_UP', 'INTERFERON_GAMMA_RESPONSE', 'COMPLEMENT', 'ALLOGRAFT_REJECTION', 'INTERFERON_ALPHA_RESPONSE',
#               'APICAL_JUNCTION', 'PROTEIN_SECRETION', 'PI3K_AKT_MTOR_SIGNALING', 'ANDROGEN_RESPONSE', 'IL2_STAT5_SIGNALING', 'APOPTOSIS',
#               'MYC_TARGETS_V2', 'G2M_CHECKPOINT', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'APICAL_SURFACE', 'XENOBIOTIC_METABOLISM', 'FATTY_ACID_METABOLISM', 
#               'ADIPOGENESIS', 'CHOLESTEROL_HOMEOSTASIS', 'PANCREAS_BETA_CELLS', 'OXIDATIVE_PHOSPHORYLATION', 'UNFOLDED_PROTEIN_RESPONSE', 
#               'UV_RESPONSE_UP', 'DNA_REPAIR', 'REACTIVE_OXYGEN_SPECIES_PATHWAY', 'MTORC1_SIGNALING', 'WNT_BETA_CATENIN_SIGNALING',
#               'HEDGEHOG_SIGNALING', 'P53_PATHWAY', 'SPERMATOGENESIS', 'PEROXISOME', 'ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE',
#               'GLYCOLYSIS', 'UV_RESPONSE_DN', 'MITOTIC_SPINDLE')

#row_order <- rev(row_order)
#t_gsvaScores <- t_gsvaScores[row_order, ]

#gsvaScores$condition <- NULL

#gsvaScores <- t(gsvaScores)
#colnames(gsvaScores)

#pdf('/Users/muje/Desktop/lab/BRL_Project/Coculture/GSVA/GSVAScores_HallmarkGeneSet_JAKSTAT.pdf')
#pheatmap(t(t_gsvaScores), cluster_row = F, cluster_cols = F, 
#         scale = 'none', treeheight_row = 10, treeheight_col = 10, 
#         cellwidth = 7, cellheight = 10, fontsize_row = 10, fontsize_col = 8, border_color = NA,
#         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', 
#         clustering_method = 'complete', clustering_callback = callback,
#         legend = T,
#         legend_breaks = c(-0.6, -0.3, 0.0, 0.3, 0.6), 
#         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
#         color = colorRampPalette(c('#3a5fcd', 'white', '#ee0000'))(n=14))
#dev.off()


###########################################

#cellular_component <- c('APICAL_JUNCTION', 'APICAL_SURFACE', 'PEROXISOME') #apical~은 development, peroxisome은 metabolic로 이동

development <- c('ADIPOGENESIS', 'ANGIOGENESIS', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 
                 'MYOGENESIS', 'SPERMATOGENESIS', 'PANCREAS_BETA_CELLS',
                 'APICAL_JUNCTION', 'APICAL_SURFACE')

DNA_damage <- c('DNA_REPAIR', 'UV_RESPONSE_DN', 'UV_RESPONSE_UP')

inflammation <- c('ALLOGRAFT_REJECTION', 'COAGULATION', 'COMPLEMENT', 'INTERFERON_ALPHA_RESPONSE',
                  'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE')

metabolic <- c('BILE_ACID_METABOLISM', 'CHOLESTEROL_HOMEOSTASIS', 'FATTY_ACID_METABOLISM', 'PEROXISOME',
               'GLYCOLYSIS', 'HEME_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION', 'XENOBIOTIC_METABOLISM')

stress_response <- c('APOPTOSIS', 'HYPOXIA', 'PROTEIN_SECRETION', 
                     'UNFOLDED_PROTEIN_RESPONSE', 'REACTIVE_OXYGEN_SPECIES_PATHWAY')

proliferation <- c('E2F_TARGETS', 'G2M_CHECKPOINT', 
                   'MYC_TARGETS_V1', 'MYC_TARGETS_V2', 'P53_PATHWAY', 'MITOTIC_SPINDLE')

signaling <- c('ANDROGEN_RESPONSE', 'ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE', 
               'IL2_STAT5_SIGNALING', 'KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN',
               'MTORC1_SIGNALING', 'NOTCH_SIGNALING', 'PI3K_AKT_MTOR_SIGNALING', 'JAK_STAT_SIGNALING', 
               'HEDGEHOG_SIGNALING', 'TGF_BETA_SIGNALING', 'TNFA_SIGNALING_VIA_NFKB', 'WNT_BETA_CATENIN_SIGNALING')


categories_list <- list(
#  Cellular_Component = cellular_component,
  Development = development,
  DNA_damage = DNA_damage,
  Inflammation = inflammation,
  Metabolic = metabolic,
  Stress_response = stress_response,
  Proliferation = proliferation,
  Signaling = signaling
)

ordered_names <- unlist(categories_list)

annotation_row <- data.frame(
  Category = rep(names(categories_list), sapply(categories_list, length))
)
rownames(annotation_row) <- ordered_names

common_rows <- intersect(ordered_names, rownames(t_gsvaScores))
t_gsvaScores_ordered <- t_gsvaScores[common_rows, ]
annotation_row_final <- annotation_row[common_rows, , drop = FALSE]

row_gaps <- cumsum(sapply(categories_list, length))


actual_categories <- annotation_row_final$Category
category_lengths <- rle(as.character(actual_categories))$lengths
row_gaps_fixed <- cumsum(category_lengths)
row_gaps_fixed <- row_gaps_fixed[row_gaps_fixed < nrow(t_gsvaScores_ordered)]

ann_colors <- list(
  Category = c(
    Development = "#f65162",
    DNA_damage = "#f59244",
    Immune = "#e1c855",
    Metabolic = "#8ac25e",
    Stress_response = "#7bd2d1",
    Proliferation = "#6ea1d1",
    Signaling = "#bf9ecb"
  )
)
annotation_row_final$Category <- as.character(annotation_row_final$Category)

sorted_names <- unlist(lapply(names(categories_list), function(cat) {
  genes <- categories_list[[cat]]
  genes <- intersect(genes, rownames(t_gsvaScores))
  variances <- apply(t_gsvaScores[genes, ], 1, var)
  genes[order(variances, decreasing = TRUE)]
}))
common_rows <- intersect(sorted_names, rownames(t_gsvaScores))
t_gsvaScores_ordered <- t_gsvaScores[common_rows, ]
annotation_row_final <- annotation_row[common_rows, , drop = FALSE]

pdf('GSVAScores_HallmarkGeneSet_JAKSTATadded.pdf', width = 10, height = 12)
pheatmap(t(t_gsvaScores_ordered), cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = annotation_row_final, gaps_col = row_gaps_fixed, scale = 'none', 
         cellwidth = 7, cellheight = 10, fontsize_row = 10, fontsize_col = 8, 
         border_color = NA, legend = TRUE, 
         annotation_colors = ann_colors,
         legend_breaks = c(-0.6, -0.3, 0.0, 0.3, 0.6), 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14))

dev.off()


###########################################

# 'PI3K_AKT_MTOR_SIGNALING', 'MTORC1_SIGNALING', 'KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN', 'TGF_BETA_SIGNALING' > proliferation
# 'IL2_STAT5_SIGNALING', 'TNFA_SIGNALING_VIA_NFKB', 'JAK_STAT_SIGNALING' > inflammation
# 'NOTCH_SIGNALING', 'HEDGEHOG_SIGNALING', 'WNT_BETA_CATENIN_SIGNALING' > development
# development와 metastasis를 구분
# 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'ANGIOGENESIS', 'APICAL_JUNCTION', 'APICAL_SURFACE' > metastasis

cell_state <- c('ADIPOGENESIS', 'MYOGENESIS', 'SPERMATOGENESIS', 'PANCREAS_BETA_CELLS',
                'NOTCH_SIGNALING', 'HEDGEHOG_SIGNALING',
                'ANDROGEN_RESPONSE', 'ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE')

#DNA_damage <- c('DNA_REPAIR', 'UV_RESPONSE_DN', 'UV_RESPONSE_UP')

metastasis <- c('EPITHELIAL_MESENCHYMAL_TRANSITION', 'ANGIOGENESIS',
                'APICAL_JUNCTION', 'APICAL_SURFACE')

inflammation <- c('ALLOGRAFT_REJECTION', 'COAGULATION', 'COMPLEMENT', 'INTERFERON_ALPHA_RESPONSE',
                  'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING', 'INFLAMMATORY_RESPONSE',
                  'IL2_STAT5_SIGNALING', 'TNFA_SIGNALING_VIA_NFKB',  'JAK_STAT_SIGNALING')

metabolic <- c('BILE_ACID_METABOLISM', 'CHOLESTEROL_HOMEOSTASIS', 'FATTY_ACID_METABOLISM', 'PEROXISOME',
               'GLYCOLYSIS', 'HEME_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION', 'XENOBIOTIC_METABOLISM')

stress_response <- c('APOPTOSIS', 'HYPOXIA', 'PROTEIN_SECRETION', 
                     'UNFOLDED_PROTEIN_RESPONSE', 'REACTIVE_OXYGEN_SPECIES_PATHWAY',
                     'DNA_REPAIR', 'UV_RESPONSE_DN', 'UV_RESPONSE_UP')

proliferation <- c('E2F_TARGETS', 'G2M_CHECKPOINT', 'P53_PATHWAY', 
                   'MYC_TARGETS_V1', 'MYC_TARGETS_V2', 'MITOTIC_SPINDLE')

oncogenic_signaling <- c('PI3K_AKT_MTOR_SIGNALING', 'MTORC1_SIGNALING', 
                         'KRAS_SIGNALING_UP', 'KRAS_SIGNALING_DN', 'TGF_BETA_SIGNALING')

#development_signaling <- c('NOTCH_SIGNALING', 'HEDGEHOG_SIGNALING', 'WNT_BETA_CATENIN_SIGNALING')

#hormone_response <- c('ANDROGEN_RESPONSE', 'ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE')


categories_list <- list(
  #  Cellular_Component = cellular_component,
  Cell_state = cell_state,
  #DNA_damage = DNA_damage,
  Metastasis = metastasis,
  Inflammation = inflammation,
  Metabolic = metabolic,
  Stress_response = stress_response,
  Proliferation = proliferation,
  Oncogenic_signaling  = oncogenic_signaling
  #Development_signaling = development_signaling,
  #Hormone_response = hormone_response
)

ordered_names <- unlist(categories_list)

annotation_row <- data.frame(
  Category = rep(names(categories_list), sapply(categories_list, length))
)
rownames(annotation_row) <- ordered_names

common_rows <- intersect(ordered_names, rownames(t_gsvaScores))
t_gsvaScores_ordered <- t_gsvaScores[common_rows, ]
annotation_row_final <- annotation_row[common_rows, , drop = FALSE]

row_gaps <- cumsum(sapply(categories_list, length))


actual_categories <- annotation_row_final$Category
category_lengths <- rle(as.character(actual_categories))$lengths
row_gaps_fixed <- cumsum(category_lengths)
row_gaps_fixed <- row_gaps_fixed[row_gaps_fixed < nrow(t_gsvaScores_ordered)]


ann_colors <- list(
  Category = c(
    Cell_state = "#f65162",
    Metastasis = "#f59244",
    Inflammation = "#e1c855",
    Metabolic = "#8ac25e",
    Stress_response = "#7bd2d1",
    Proliferation = "#6ea1d1",
    Oncogenic_signaling = "#bf9ecb"
  )
)

annotation_row_final$Category <- as.character(annotation_row_final$Category)

sorted_names <- unlist(lapply(names(categories_list), function(cat) {
  genes <- categories_list[[cat]]
  genes <- intersect(genes, rownames(t_gsvaScores))
  variances <- apply(t_gsvaScores[genes, ], 1, var)
  genes[order(variances, decreasing = TRUE)]
}))
common_rows <- intersect(sorted_names, rownames(t_gsvaScores))
t_gsvaScores_ordered <- t_gsvaScores[common_rows, ]
annotation_row_final <- annotation_row[common_rows, , drop = FALSE]

pdf('GSVAScores_HallmarkGeneSet_newcategory.pdf', width = 10, height = 12)
pheatmap(t(t_gsvaScores_ordered), cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_col = annotation_row_final, gaps_col = row_gaps_fixed, scale = 'none', 
         cellwidth = 7, cellheight = 10, fontsize_row = 10, fontsize_col = 8, 
         border_color = NA, legend = TRUE, 
         annotation_colors = ann_colors,
         legend_breaks = c(-0.6, -0.3, 0.0, 0.3, 0.6), 
         breaks = c(-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
         color = colorRampPalette(c("#3a5fcd", 'white', "#ee0000"))(14))

dev.off()

