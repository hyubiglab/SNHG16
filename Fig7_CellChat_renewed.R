###### Heatmap visualization ######
library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# 1) load
cellchat <- readRDS("~/Projects/BRL/SNHG16/allcells/cell_types/tmp/cellchat.Rds")

sender_cell   <- "Macrophages"
receiver_cell <- "Malignant cells"

# 2) direct LR pairs: Macrophages -> Malignant cells
df.direct <- subsetCommunication(
  object      = cellchat,
  slot.name   = "net",
  sources.use = sender_cell,
  targets.use = receiver_cell
)

df.direct <- df.direct[df.direct$prob > 0, , drop = FALSE]

if (nrow(df.direct) == 0) {
  stop("No inferred LR pairs found for Macrophages -> Malignant cells.")
}

selected_pairs <- unique(df.direct$interaction_name)

pair_order_df <- df.direct[order(-df.direct$prob, df.direct$pval), ]
pair_order <- unique(pair_order_df$interaction_name)

pair_annot <- pair_order_df[!duplicated(pair_order_df$interaction_name),
                            c("interaction_name", "ligand", "receptor", "pathway_name")]
pair_labels <- paste0(
  pair_annot$ligand, " -> ", pair_annot$receptor, " (", pair_annot$pathway_name, ")"
)
names(pair_labels) <- pair_annot$interaction_name

# 3) all LR communications across all cell groups
df.all <- subsetCommunication(
  object    = cellchat,
  slot.name = "net"
)

df.all <- df.all[df.all$prob > 0, , drop = FALSE]
df.sel <- df.all[df.all$interaction_name %in% selected_pairs, , drop = FALSE]

if (nrow(df.sel) == 0) {
  stop("Selected LR pairs were not found in the full communication table.")
}

# 4) 모든 cell group 이름
cell_groups <- levels(cellchat@idents)

# 5) outgoing / incoming matrix
out.df <- aggregate(prob ~ interaction_name + source, data = df.sel, sum)
mat.out.ori <- xtabs(prob ~ interaction_name + source, data = out.df)

in.df <- aggregate(prob ~ interaction_name + target, data = df.sel, sum)
mat.in.ori <- xtabs(prob ~ interaction_name + target, data = in.df)

add_missing_cols <- function(mat, wanted_cols) {
  missing_cols <- setdiff(wanted_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    add <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      dimnames = list(rownames(mat), missing_cols)
    )
    mat <- cbind(mat, add)
  }
  mat[, wanted_cols, drop = FALSE]
}

add_missing_rows <- function(mat, wanted_rows) {
  missing_rows <- setdiff(wanted_rows, rownames(mat))
  if (length(missing_rows) > 0) {
    add <- matrix(
      0,
      nrow = length(missing_rows),
      ncol = ncol(mat),
      dimnames = list(missing_rows, colnames(mat))
    )
    mat <- rbind(mat, add)
  }
  mat[wanted_rows, , drop = FALSE]
}

mat.out.ori <- add_missing_cols(mat.out.ori, cell_groups)
mat.in.ori  <- add_missing_cols(mat.in.ori,  cell_groups)

mat.out.ori <- add_missing_rows(mat.out.ori, pair_order)
mat.in.ori  <- add_missing_rows(mat.in.ori,  pair_order)

rownames(mat.out.ori) <- pair_labels[rownames(mat.out.ori)]
rownames(mat.in.ori)  <- pair_labels[rownames(mat.in.ori)]

# 6) CellChat style row scaling
row_scale_by_max <- function(mat) {
  maxv <- apply(mat, 1, max, na.rm = TRUE)
  maxv[maxv == 0] <- 1
  sweep(mat, 1, maxv, "/")
}

mat.out <- row_scale_by_max(mat.out.ori)
mat.in  <- row_scale_by_max(mat.in.ori)

# 7) CellChat style color
color.heatmap.use <- colorRampPalette(brewer.pal(n = 9, name = "BuGn"))(100)

# cell color bar용 색
cell_colors <- setNames(
  CellChat::scPalette(length(cell_groups)),
  cell_groups
)

# column annotation (아래 색 띠)
col_df <- data.frame(group = colnames(mat.out))
rownames(col_df) <- colnames(mat.out)

col_annotation <- HeatmapAnnotation(
  df = col_df,
  col = list(group = cell_colors),
  which = "column",
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.2, "cm")
)

# CellChat source와 비슷한 right bar 계산
make_right_bar <- function(mat_ori) {
  pSum <- rowSums(mat_ori)
  pSum.original <- pSum
  pSum <- -1 / log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(seq_along(idx1), position)]
  }
  pSum
}

# outgoing annotations
ha2.out <- HeatmapAnnotation(
  Strength = anno_barplot(
    colSums(mat.out.ori),
    border = FALSE,
    gp = gpar(fill = cell_colors[colnames(mat.out.ori)], col = cell_colors[colnames(mat.out.ori)])
  ),
  show_annotation_name = FALSE
)

ha1.out <- rowAnnotation(
  Strength = anno_barplot(make_right_bar(mat.out.ori), border = FALSE),
  show_annotation_name = FALSE
)

# incoming annotations
ha2.in <- HeatmapAnnotation(
  Strength = anno_barplot(
    colSums(mat.in.ori),
    border = FALSE,
    gp = gpar(fill = cell_colors[colnames(mat.in.ori)], col = cell_colors[colnames(mat.in.ori)])
  ),
  show_annotation_name = FALSE
)

ha1.in <- rowAnnotation(
  Strength = anno_barplot(make_right_bar(mat.in.ori), border = FALSE),
  show_annotation_name = FALSE
)

pdf(
  "~/Projects/BRL/SNHG16/allcells/cell_types/pathways/LRpairs_MacrophageOut_MalignantIn_allCells_heatmap_withBars.pdf",
  width = max(10, 0.55 * ncol(mat.out) + 8),
  height = max(6, 0.22 * nrow(mat.out) + 2)
)

ht.out <- Heatmap(
  mat.out,
  col = color.heatmap.use,
  na_col = "white",
  name = "Relative strength",
  bottom_annotation = col_annotation,
  top_annotation = ha2.out,
  right_annotation = ha1.out,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_rot = 0,
  column_names_rot = 90,
  column_title = "Outgoing signaling patterns"
)

ht.in <- Heatmap(
  mat.in,
  col = color.heatmap.use,
  na_col = "white",
  name = "Relative strength",
  bottom_annotation = col_annotation,
  top_annotation = ha2.in,
  right_annotation = ha1.in,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_rot = 0,
  column_names_rot = 90,
  column_title = "Incoming signaling patterns"
)

draw(ht.out + ht.in, heatmap_legend_side = "right")
dev.off()



###### JAK/STAT3 narrowing (meanScore version for transcript-like features; robust matching) ######
outdir <- "/home/minwook/Projects/BRL/SNHG16"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1) load
obj <- readRDS(file.path(outdir, "crc.mtcut.norm.harmony.Rds"))

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Cell_type_v2"

# sanity check
print(table(obj$Cell_type_v2))
print(head(obj@meta.data[, c("Patient", "Class", "Condition", "Cell_type_v2")]))

# 2) helper: Seurat v4/v5 compatible data extraction
get_norm_mat <- function(seu, assay = "RNA") {
  tryCatch({
    GetAssayData(seu, assay = assay, layer = "data")
  }, error = function(e) {
    GetAssayData(seu, assay = assay, slot = "data")
  })
}

# 3) helper: regex escape
escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# 4) helper: detect gene-symbol column in meta.features
get_meta_symbol_col <- function(seu, assay = "RNA") {
  mf <- seu[[assay]]@meta.features
  candidate_cols <- intersect(
    c("gene_symbol", "gene_name", "symbol", "gene", "Gene", "external_gene_name"),
    colnames(mf)
  )
  if (length(candidate_cols) == 0) return(NULL)
  candidate_cols[1]
}

# 5) helper: feature -> symbol map for export / DEG annotation
#    meta.features symbol이 있으면 그대로 사용, 없으면 fallback 정리
make_feature_map <- function(seu, assay = "RNA") {
  features <- rownames(seu[[assay]])
  mf <- seu[[assay]]@meta.features
  symbol_col <- get_meta_symbol_col(seu, assay = assay)
  
  if (!is.null(symbol_col)) {
    symbols <- as.character(mf[[symbol_col]])
  } else {
    symbols <- features
    symbols <- sub("\\.tc$", "", symbols)
    symbols <- sub("-[0-9]+$", "", symbols)
  }
  
  out <- data.frame(
    feature = features,
    symbol = symbols,
    stringsAsFactors = FALSE
  )
  
  out <- out[!is.na(out$symbol) & out$symbol != "", , drop = FALSE]
  out
}

# 6) helper: robust symbol aggregation for transcript-like features
#    priority:
#    (1) meta.features exact symbol match
#    (2) fallback prefix regex match: ^STAT3($|-|\\.)
aggregate_symbol_expr_robust <- function(seu, assay = "RNA", symbols) {
  mat <- get_norm_mat(seu, assay = assay)
  features <- rownames(mat)
  mf <- seu[[assay]]@meta.features
  symbol_col <- get_meta_symbol_col(seu, assay = assay)
  
  symbols <- unique(symbols)
  out_list <- list()
  match_list <- list()
  
  for (sym in symbols) {
    feats <- character(0)
    
    # 1) exact match from meta.features symbol column
    if (!is.null(symbol_col)) {
      meta_symbols <- as.character(mf[[symbol_col]])
      feats <- features[toupper(meta_symbols) == toupper(sym)]
    }
    
    # 2) fallback prefix regex on rownames
    if (length(feats) == 0) {
      pat <- paste0("^", escape_regex(sym), "(?:$|-|\\.)")
      feats <- grep(
        pat,
        features,
        value = TRUE,
        ignore.case = TRUE,
        perl = TRUE
      )
    }
    
    feats <- intersect(unique(feats), rownames(mat))
    if (length(feats) == 0) next
    
    vals <- if (length(feats) == 1) {
      as.numeric(mat[feats, ])
    } else {
      Matrix::colMeans(mat[feats, , drop = FALSE])
    }
    
    names(vals) <- colnames(mat)
    out_list[[sym]] <- vals
    
    match_list[[sym]] <- data.frame(
      symbol = sym,
      n_features = length(feats),
      features = paste(feats, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(out_list) == 0) {
    out_mat <- matrix(
      numeric(0),
      nrow = 0,
      ncol = ncol(mat),
      dimnames = list(NULL, colnames(mat))
    )
    match_df <- data.frame(
      symbol = character(0),
      n_features = integer(0),
      features = character(0),
      stringsAsFactors = FALSE
    )
    return(list(mat = out_mat, match_table = match_df))
  }
  
  out_mat <- do.call(rbind, out_list)
  
  if (is.null(dim(out_mat))) {
    out_mat <- matrix(
      out_mat,
      nrow = 1,
      dimnames = list(names(out_list), colnames(mat))
    )
  } else {
    rownames(out_mat) <- names(out_list)
    colnames(out_mat) <- colnames(mat)
  }
  
  match_df <- do.call(rbind, match_list)
  rownames(match_df) <- NULL
  
  list(mat = out_mat, match_table = match_df)
}

# 7) helper: quick diagnostic for transcript-like naming
inspect_gene_patterns <- function(seu, genes, assay = "RNA", n_show = 10) {
  features <- rownames(seu[[assay]])
  out <- lapply(genes, function(g) {
    pat <- paste0("^", escape_regex(g), "(?:$|-|\\.)")
    hits <- grep(pat, features, value = TRUE, ignore.case = TRUE, perl = TRUE)
    data.frame(
      query_gene = g,
      matched_feature = head(hits, n_show),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(out)
  if (nrow(out) == 0) {
    out <- data.frame(
      query_gene = genes,
      matched_feature = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  out
}

# 8) export feature map
feature_map <- make_feature_map(obj, assay = "RNA")
write.table(
  feature_map,
  file.path(outdir, "RNA_feature_to_symbol_map.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# optional diagnostic
diag_genes <- c("STAT3", "JAK1", "JAK2", "IL6ST", "SOCS3", "MYC", "ICAM1", "VEGFA")
diag_tbl <- inspect_gene_patterns(obj, diag_genes, assay = "RNA", n_show = 15)
write.table(
  diag_tbl,
  file.path(outdir, "RNA_transcriptLike_feature_diagnostic.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 9) Hallmark IL6/JAK/STAT3 gene set
msig <- msigdbr::msigdbr(species = "Homo sapiens")
hallmark_jak.raw <- unique(
  msig$gene_symbol[msig$gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING"]
)

hallmark_res <- aggregate_symbol_expr_robust(
  seu = obj,
  assay = "RNA",
  symbols = hallmark_jak.raw
)
hallmark_gene_mat <- hallmark_res$mat

cat("Matched hallmark genes:", nrow(hallmark_gene_mat), "\n")
print(head(rownames(hallmark_gene_mat), 30))

write.table(
  hallmark_res$match_table,
  file.path(outdir, "JAKSTAT3_hallmark_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(hallmark_gene_mat) < 5) {
  unmatched <- setdiff(hallmark_jak.raw, rownames(hallmark_gene_mat))
  write.table(
    data.frame(symbol = unmatched, stringsAsFactors = FALSE),
    file.path(outdir, "JAKSTAT3_hallmark_unmatched_symbols.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  stop("Too few matched hallmark genes even after robust matching. Check diagnostic output and RNA assay content.")
}

# 10) cell-level JAK/STAT3 mean score
obj$JAKSTAT3_meanScore <- colMeans(hallmark_gene_mat)
score_col <- "JAKSTAT3_meanScore"

# 11) subsets
mal <- subset(obj, subset = Cell_type_v2 == "Malignant cells")
mac <- subset(obj, subset = Cell_type_v2 == "Macrophages")

DefaultAssay(mal) <- "RNA"
DefaultAssay(mac) <- "RNA"

mal_feature_map <- make_feature_map(mal, assay = "RNA")
mac_feature_map <- make_feature_map(mac, assay = "RNA")

# 12) overview plots
red_use <- if ("umap" %in% Reductions(obj)) "umap" else Reductions(obj)[1]

p1 <- FeaturePlot(
  obj,
  features = score_col,
  reduction = red_use,
  cols = c("grey90", "red")
) + ggtitle("JAK/STAT3 mean score")

p2 <- VlnPlot(
  obj,
  features = score_col,
  group.by = "Cell_type_v2",
  pt.size = 0
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(outdir, "JAKSTAT3_overview.allcells.meanScore.pdf"), width = 16, height = 6)
print(p1 + p2)
dev.off()

p3 <- FeaturePlot(
  mal,
  features = score_col,
  reduction = red_use,
  cols = c("grey90", "red")
) + ggtitle("Malignant only: JAK/STAT3 mean score")

pdf(file.path(outdir, "JAKSTAT3_overview.malignant.meanScore.pdf"), width = 7, height = 6)
print(p3)
dev.off()

# 13) receptor-high vs low analysis in malignant cells
receptors_raw <- c("BSG", "ADGRE5", "CD44", "P4HB", "PLAUR", "TNFRSF1A")

receptor_res <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = receptors_raw
)
receptor_mat <- receptor_res$mat

write.table(
  receptor_res$match_table,
  file.path(outdir, "malignant_receptor_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(receptor_mat) == 0) {
  stop("No receptor genes matched in malignant RNA assay.")
}

receptor_test_list <- list()

for (g in rownames(receptor_mat)) {
  vals <- as.numeric(receptor_mat[g, ])
  names(vals) <- colnames(mal)
  
  status <- ifelse(vals > 0, "High", "Low")
  split_method <- "detected_gt0"
  
  if (length(unique(status)) < 2) {
    if (length(unique(vals)) < 2) next
    cutoff <- median(vals, na.rm = TRUE)
    status <- ifelse(vals > cutoff, "High", "Low")
    split_method <- "median_split"
  }
  
  df <- data.frame(
    cell = colnames(mal),
    receptor = g,
    expr = vals,
    status = status,
    split_method = split_method,
    JAKSTAT3 = mal@meta.data[colnames(mal), score_col],
    stringsAsFactors = FALSE
  )
  
  if (length(unique(df$status)) < 2) next
  
  wt <- wilcox.test(JAKSTAT3 ~ status, data = df)
  
  receptor_test_list[[g]] <- data.frame(
    receptor = g,
    split_method = split_method,
    n_high = sum(df$status == "High"),
    n_low = sum(df$status == "Low"),
    median_high = median(df$JAKSTAT3[df$status == "High"], na.rm = TRUE),
    median_low = median(df$JAKSTAT3[df$status == "Low"], na.rm = TRUE),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
  
  pp <- ggplot(df, aes(status, JAKSTAT3, fill = status)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    theme_classic() +
    ggtitle(paste0(g, " in malignant cells")) +
    labs(
      subtitle = paste0(split_method, " | Wilcoxon p = ", signif(wt$p.value, 3)),
      x = NULL,
      y = "JAK/STAT3 mean score"
    )
  
  ggsave(
    filename = file.path(outdir, paste0("receptor_", g, "_vs_JAKSTAT3.malignant.meanScore.pdf")),
    plot = pp, width = 4.5, height = 4.5
  )
}

receptor_test_res <- bind_rows(receptor_test_list)

if (nrow(receptor_test_res) > 0) {
  receptor_test_res <- receptor_test_res %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value)
}

write.table(
  receptor_test_res,
  file.path(outdir, "receptor_vs_JAKSTAT3.malignant.meanScore.stats.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(receptor_test_res)

# 14) macrophage ligand mean score
ligands_raw <- c("PPIA", "CD55", "LGALS9", "PLAU", "SPP1", "TNF")

ligand_res <- aggregate_symbol_expr_robust(
  seu = obj,
  assay = "RNA",
  symbols = ligands_raw
)
ligand_mat_all <- ligand_res$mat

write.table(
  ligand_res$match_table,
  file.path(outdir, "macrophage_ligand_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(ligand_mat_all) == 0) {
  stop("No ligand genes matched in the RNA assay.")
}

obj$MacrophageLigand_meanScore <- colMeans(ligand_mat_all)
lig_score_col <- "MacrophageLigand_meanScore"

# refresh subsets after metadata update
mal <- subset(obj, subset = Cell_type_v2 == "Malignant cells")
mac <- subset(obj, subset = Cell_type_v2 == "Macrophages")

DefaultAssay(mal) <- "RNA"
DefaultAssay(mac) <- "RNA"

mal_feature_map <- make_feature_map(mal, assay = "RNA")
mac_feature_map <- make_feature_map(mac, assay = "RNA")

# 15) patient-level summary:
# macrophage ligand score vs malignant JAK/STAT3 score
mac_df <- as.data.frame(mac@meta.data)
mac_df$cell <- rownames(mac_df)
mac_df$lig_score_use <- mac_df[[lig_score_col]]

mal_df <- as.data.frame(mal@meta.data)
mal_df$cell <- rownames(mal_df)
mal_df$jak_score_use <- mal_df[[score_col]]

mac_patient <- mac_df %>%
  dplyr::group_by(Patient, Class, Condition) %>%
  dplyr::summarise(
    macrophage_ligand_score = mean(lig_score_use, na.rm = TRUE),
    n_macrophage = length(lig_score_use),
    .groups = "drop"
  )

mal_patient <- mal_df %>%
  dplyr::group_by(Patient, Class, Condition) %>%
  dplyr::summarise(
    malignant_JAKSTAT3_score = mean(jak_score_use, na.rm = TRUE),
    n_malignant = length(jak_score_use),
    .groups = "drop"
  )

patient_pair <- inner_join(
  mac_patient,
  mal_patient,
  by = c("Patient", "Class", "Condition")
)

write.table(
  patient_pair,
  file.path(outdir, "macrophageLigand_vs_malignantJAKSTAT3.patientTable.meanScore.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(patient_pair) >= 3) {
  sp <- suppressWarnings(cor.test(
    patient_pair$macrophage_ligand_score,
    patient_pair$malignant_JAKSTAT3_score,
    method = "spearman"
  ))
  
  p_corr <- ggplot(
    patient_pair,
    aes(x = macrophage_ligand_score, y = malignant_JAKSTAT3_score, label = Patient)
  ) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(vjust = -0.6, size = 3) +
    theme_classic() +
    labs(
      x = "Macrophage ligand mean score",
      y = "Malignant JAK/STAT3 mean score",
      subtitle = paste0(
        "Spearman rho = ", round(sp$estimate, 3),
        ", p = ", signif(sp$p.value, 3)
      )
    )
  
  ggsave(
    file.path(outdir, "macrophageLigand_vs_malignantJAKSTAT3.corr.meanScore.pdf"),
    p_corr, width = 5.5, height = 5
  )
}

# 16) pair-by-pair patient-level analysis
pair_list <- list(
  PPIA_BSG = c("PPIA", "BSG"),
  CD55_ADGRE5 = c("CD55", "ADGRE5"),
  LGALS9_CD44 = c("LGALS9", "CD44"),
  LGALS9_P4HB = c("LGALS9", "P4HB"),
  PLAU_PLAUR = c("PLAU", "PLAUR"),
  SPP1_CD44 = c("SPP1", "CD44"),
  TNF_TNFRSF1A = c("TNF", "TNFRSF1A")
)

all_ligands <- unique(vapply(pair_list, `[`, character(1), 1))
all_receptors <- unique(vapply(pair_list, `[`, character(1), 2))

mac_symbol_res <- aggregate_symbol_expr_robust(
  seu = mac,
  assay = "RNA",
  symbols = all_ligands
)
mac_symbol_mat <- mac_symbol_res$mat

mal_symbol_res <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = all_receptors
)
mal_symbol_mat <- mal_symbol_res$mat

write.table(
  mac_symbol_res$match_table,
  file.path(outdir, "macrophage_pairLigands_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  mal_symbol_res$match_table,
  file.path(outdir, "malignant_pairReceptors_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pair_res <- list()

for (nm in names(pair_list)) {
  lg <- pair_list[[nm]][1]
  rc <- pair_list[[nm]][2]
  
  if (!(lg %in% rownames(mac_symbol_mat))) next
  if (!(rc %in% rownames(mal_symbol_mat))) next
  
  mac_tmp <- data.frame(
    Patient = mac$Patient,
    ligand_expr = as.numeric(mac_symbol_mat[lg, ]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      ligand_expr = mean(ligand_expr, na.rm = TRUE),
      .groups = "drop"
    )
  
  mal_tmp <- data.frame(
    Patient = mal$Patient,
    receptor_expr = as.numeric(mal_symbol_mat[rc, ]),
    JAKSTAT3 = mal@meta.data[, score_col],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      receptor_expr = mean(receptor_expr, na.rm = TRUE),
      JAKSTAT3 = mean(JAKSTAT3, na.rm = TRUE),
      .groups = "drop"
    )
  
  tmp <- inner_join(mac_tmp, mal_tmp, by = "Patient")
  if (nrow(tmp) < 3) next
  
  cor_lig_jak <- suppressWarnings(cor.test(
    tmp$ligand_expr, tmp$JAKSTAT3, method = "spearman"
  ))
  cor_rec_jak <- suppressWarnings(cor.test(
    tmp$receptor_expr, tmp$JAKSTAT3, method = "spearman"
  ))
  
  pair_res[[nm]] <- data.frame(
    pair = nm,
    ligand = lg,
    receptor = rc,
    n_patient = nrow(tmp),
    rho_ligand_vs_JAKSTAT3 = unname(cor_lig_jak$estimate),
    p_ligand_vs_JAKSTAT3 = cor_lig_jak$p.value,
    rho_receptor_vs_JAKSTAT3 = unname(cor_rec_jak$estimate),
    p_receptor_vs_JAKSTAT3 = cor_rec_jak$p.value,
    stringsAsFactors = FALSE
  )
}

pair_res_df <- bind_rows(pair_res)

if (nrow(pair_res_df) > 0) {
  pair_res_df$FDR_ligand <- p.adjust(pair_res_df$p_ligand_vs_JAKSTAT3, method = "BH")
  pair_res_df$FDR_receptor <- p.adjust(pair_res_df$p_receptor_vs_JAKSTAT3, method = "BH")
}

write.table(
  pair_res_df,
  file.path(outdir, "LRpair_vs_JAKSTAT3.patientCorrelation.meanScore.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(pair_res_df)

# 17) malignant DE between high vs low JAK/STAT3 mean score

jak_vals <- mal@meta.data[[score_col]]

q_hi <- quantile(jak_vals, 0.75, na.rm = TRUE)
q_lo <- quantile(jak_vals, 0.25, na.rm = TRUE)

jak_group <- rep(NA_character_, length(jak_vals))
jak_group[jak_vals >= q_hi] <- "High"
jak_group[jak_vals <= q_lo] <- "Low"

mal@meta.data$JAKSTAT3_group <- jak_group

print(table(mal@meta.data$JAKSTAT3_group, useNA = "ifany"))

cells_use <- rownames(mal@meta.data)[!is.na(mal@meta.data$JAKSTAT3_group)]
mal2 <- subset(mal, cells = cells_use)

Idents(mal2) <- "JAKSTAT3_group"

deg_jak <- FindMarkers(
  mal2,
  ident.1 = "High",
  ident.2 = "Low",
  assay = "RNA",
  logfc.threshold = 0.1,
  min.pct = 0.1
)

deg_jak$feature <- rownames(deg_jak)

mal2_feature_map <- make_feature_map(mal2, assay = "RNA")
feature_to_symbol <- setNames(mal2_feature_map$symbol, mal2_feature_map$feature)
deg_jak$gene_symbol <- unname(feature_to_symbol[deg_jak$feature])

write.table(
  deg_jak,
  file.path(outdir, "malignant_JAKSTAT3High_vs_Low.DEG.meanScore.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 18) save object with scores
saveRDS(obj, file.path(outdir, "crc.mtcut.norm.harmony.JAKSTAT3_meanScore.Rds"))





###### Redefining high and low ######
receptors_raw <- c("BSG", "ADGRE5", "CD44", "P4HB", "PLAUR", "TNFRSF1A")

receptor_res <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = receptors_raw
)
receptor_mat <- receptor_res$mat

write.table(
  receptor_res$match_table,
  file.path(outdir, "malignant_receptor_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(receptor_mat) == 0) {
  stop("No receptor genes matched in malignant RNA assay.")
}

# p-value label helper
format_p_label <- function(p) {
  if (is.na(p)) {
    return("p=NA")
  } else if (p == 0) {
    return("p<2.2e-308")
  } else if (p < 1e-300) {
    return("p<1e-300")
  } else {
    return(paste0("p=", format(p, digits = 2, scientific = TRUE)))
  }
}

receptor_test_list <- list()
receptor_df_list <- list()

receptor_order <- receptors_raw[receptors_raw %in% rownames(receptor_mat)]

for (g in receptor_order) {
  vals <- as.numeric(receptor_mat[g, ])
  names(vals) <- colnames(mal)
  
  q1 <- as.numeric(quantile(vals, 0.25, na.rm = TRUE))
  q3 <- as.numeric(quantile(vals, 0.75, na.rm = TRUE))
  
  # middle 50% 제외
  keep <- vals <= q1 | vals >= q3
  
  # usable cell이 너무 적으면 skip
  if (sum(keep, na.rm = TRUE) < 2) next
  
  vals_use <- vals[keep]
  cells_use <- names(vals)[keep]
  
  status <- ifelse(vals_use >= q3, "High", "Low")
  split_method <- "quartile_split_Q1Q3"
  
  df <- data.frame(
    cell = cells_use,
    receptor = g,
    expr = vals_use,
    status = status,
    split_method = split_method,
    JAKSTAT3 = mal@meta.data[cells_use, score_col],
    stringsAsFactors = FALSE
  )
  
  # q1 == q3 등으로 한 군만 생기면 skip
  if (length(unique(df$status)) < 2) next
  
  wt <- wilcox.test(JAKSTAT3 ~ status, data = df, exact = FALSE)
  
  median_high <- median(df$JAKSTAT3[df$status == "High"], na.rm = TRUE)
  median_low  <- median(df$JAKSTAT3[df$status == "Low"], na.rm = TRUE)
  delta_median <- median_high - median_low
  
  n_high <- sum(df$status == "High")
  n_low  <- sum(df$status == "Low")
  
  receptor_test_list[[g]] <- data.frame(
    receptor = g,
    split_method = split_method,
    q1_expr = q1,
    q3_expr = q3,
    n_high = n_high,
    n_low = n_low,
    median_high = median_high,
    median_low = median_low,
    delta_median = delta_median,
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
  
  receptor_df_list[[g]] <- df
  
  pp <- ggplot(df, aes(status, JAKSTAT3, fill = status)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    theme_classic() +
    ggtitle(paste0(g, " in malignant cells")) +
    labs(
      subtitle = paste0(
        split_method,
        " | ",
        format_p_label(wt$p.value),
        " | Δmed=",
        round(delta_median, 3),
        " | nL=",
        n_low,
        ", nH=",
        n_high
      ),
      x = NULL,
      y = "JAK/STAT3 mean score"
    )
  
  ggsave(
    filename = file.path(outdir, paste0("receptor_", g, "_vs_JAKSTAT3.malignant.meanScore.pdf")),
    plot = pp, width = 4.5, height = 4.5
  )
}

receptor_test_res <- bind_rows(receptor_test_list)

if (nrow(receptor_test_res) > 0) {
  receptor_test_res <- receptor_test_res %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value)
}

write.table(
  receptor_test_res,
  file.path(outdir, "receptor_vs_JAKSTAT3.malignant.meanScore.stats.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(receptor_test_res)

receptor_df_all <- bind_rows(receptor_df_list)

if (nrow(receptor_df_all) > 0) {
  receptor_df_all$receptor <- factor(receptor_df_all$receptor, levels = receptor_order)
  receptor_df_all$status <- factor(receptor_df_all$status, levels = c("Low", "High"))
  
  label_df <- receptor_df_all %>%
    dplyr::group_by(receptor) %>%
    dplyr::summarise(
      y_max = max(JAKSTAT3, na.rm = TRUE),
      y_min = min(JAKSTAT3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      receptor_test_res %>%
        dplyr::select(receptor, p_value, FDR, split_method, delta_median, n_low, n_high),
      by = "receptor"
    ) %>%
    dplyr::mutate(
      y_span = pmax(y_max - y_min, 0.05),
      y_pos = y_max + y_span * 0.18,
      label = paste0(
        vapply(p_value, format_p_label, character(1)),
        "\nΔmed=", round(delta_median, 3),
        "\nLow n=", n_low, ", High n=", n_high
      )
    )
  
  p_receptor_onepanel <- ggplot(
    receptor_df_all,
    aes(x = receptor, y = JAKSTAT3, fill = status)
  ) +
    geom_violin(
      position = position_dodge(width = 0.8),
      scale = "width",
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_text(
      data = label_df,
      aes(x = receptor, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3,
      lineheight = 0.9
    ) +
    theme_classic() +
    labs(
      title = "Receptor-high vs low JAK/STAT3 score in malignant cells",
      subtitle = "Low = bottom quartile, High = top quartile; middle 50% excluded",
      x = NULL,
      y = "JAK/STAT3 mean score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    filename = file.path(outdir, "receptor_all_vs_JAKSTAT3.malignant.meanScore.onepanel.pdf"),
    plot = p_receptor_onepanel,
    width = 9,
    height = 5
  )
  
  print(p_receptor_onepanel)
}



###### Final fig formatting ######
if (nrow(receptor_df_all) > 0) {
  receptor_df_all$receptor <- factor(receptor_df_all$receptor, levels = receptor_order)
  receptor_df_all$status <- factor(receptor_df_all$status, levels = c("Low", "High"))
  
  # p-value -> significance stars
  p_to_star <- function(p) {
    if (is.na(p)) {
      "NA"
    } else if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      "ns"
    }
  }
  
  label_df <- receptor_df_all %>%
    dplyr::group_by(receptor) %>%
    dplyr::summarise(
      y_max = max(JAKSTAT3, na.rm = TRUE),
      y_min = min(JAKSTAT3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      receptor_test_res %>%
        dplyr::select(receptor, p_value),
      by = "receptor"
    ) %>%
    dplyr::mutate(
      receptor_num = match(receptor, receptor_order),
      y_span = pmax(y_max - y_min, 0.05),
      y_tick = y_max + y_span * 0.05,
      y_bracket = y_max + y_span * 0.12,
      y_text = y_max + y_span * 0.18,
      x_left = receptor_num - 0.20,
      x_right = receptor_num + 0.20,
      label = vapply(p_value, p_to_star, character(1))
    )
  
  p_receptor_onepanel_star <- ggplot(
    receptor_df_all,
    aes(x = receptor, y = JAKSTAT3, fill = status)
  ) +
    geom_violin(
      position = position_dodge(width = 0.8),
      scale = "width",
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_segment(
      data = label_df,
      aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_segment(
      data = label_df,
      aes(x = x_left, xend = x_left, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_segment(
      data = label_df,
      aes(x = x_right, xend = x_right, y = y_tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_text(
      data = label_df,
      aes(x = receptor_num, y = y_text, label = label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = "JAK/STAT3 mean score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    filename = file.path(outdir, "receptor_all_vs_JAKSTAT3.malignant.meanScore.onepanel.stars_bracket.pdf"),
    plot = p_receptor_onepanel_star,
    width = 7,
    height = 2
  )
  
  print(p_receptor_onepanel_star)
}



###### Progeny faster ver ######
# Receptor high / low group
receptors_of_interest <- c("BSG", "ADGRE5", "CD44", "P4HB", "PLAUR", "TNFRSF1A")

receptor_mat <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = receptors_of_interest
)$mat

if (nrow(receptor_mat) < 2) {
  stop("Too few matched receptors in malignant cells.")
}

mal$receptor_combined_score <- colMeans(receptor_mat)

q_hi <- quantile(mal$receptor_combined_score, 0.75, na.rm = TRUE)
q_lo <- quantile(mal$receptor_combined_score, 0.25, na.rm = TRUE)

group <- rep(NA_character_, length(colnames(mal)))
names(group) <- colnames(mal)
group[mal$receptor_combined_score >= q_hi] <- "High"
group[mal$receptor_combined_score <= q_lo] <- "Low"
mal$receptor_group <- group[colnames(mal)]

cat("Receptor group counts:\n")
print(table(mal$receptor_group, useNA = "ifany"))

# Fast gene-level matrix for PROGENy
# Fast gene-level matrix for PROGENy
mal_mat <- get_norm_mat(mal, assay = "RNA")

# PROGENy model
model500 <- progeny::getModel("Human", top = 500)
rownames(model500) <- toupper(rownames(model500))

pathways_use <- c("JAK-STAT", "NFkB", "MAPK", "PI3K", "TNFa")
pathways_use <- intersect(pathways_use, colnames(model500))

if (length(pathways_use) == 0) {
  stop("Requested pathways not found in PROGENy model.")
}

# only requested pathway targets
target_symbols <- unique(rownames(model500[, pathways_use, drop = FALSE]))
target_symbols <- unique(toupper(target_symbols))

cat("Requested PROGENy target symbols:", length(target_symbols), "\n")

# robust aggregation directly from Seurat object
gene_res <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = target_symbols
)

gene_mat <- gene_res$mat

if (is.null(gene_mat) || length(gene_mat) == 0 || is.null(dim(gene_mat)) || nrow(gene_mat) == 0) {
  stop("No PROGENy target genes matched in malignant RNA assay after robust aggregation.")
}

rownames(gene_mat) <- toupper(rownames(gene_mat))
colnames(gene_mat) <- colnames(mal_mat)

cat("Gene-level matrix:", nrow(gene_mat), "genes x", ncol(gene_mat), "cells\n")

# overlap check
common_genes <- intersect(rownames(gene_mat), rownames(model500))
cat("Common genes with PROGENy model:", length(common_genes), "\n")

if (length(common_genes) < 20) {
  write.table(
    gene_res$match_table,
    file.path(outdir, "PROGENy_target_matched_features.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  stop("Too few overlapping genes between gene_mat and PROGENy model.")
}

# Fast PROGENy scoring
# result: cells x pathways
expr_use  <- gene_mat[common_genes, , drop = FALSE]
model_use <- as.matrix(model500[common_genes, pathways_use, drop = FALSE])

progeny_scores <- t(expr_use) %*% model_use
progeny_scores <- as.matrix(progeny_scores)

cat("Available PROGENy pathways:\n")
print(colnames(progeny_scores))

print(apply(progeny_scores, 2, function(x) {
  c(min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    sd  = sd(x,  na.rm = TRUE),
    nonzero = sum(x != 0, na.rm = TRUE))
}))


# Pathway selection
plot_df <- do.call(rbind, lapply(pathways_use, function(pw) {
  data.frame(
    pathway = pw,
    cell = rownames(progeny_scores),
    score = as.numeric(progeny_scores[, pw]),
    group = mal$receptor_group[rownames(progeny_scores)],
    stringsAsFactors = FALSE
  )
}))

plot_df <- plot_df[!is.na(plot_df$group), , drop = FALSE]


# Comparison
progeny_res <- lapply(unique(plot_df$pathway), function(pw) {
  df <- plot_df[plot_df$pathway == pw, , drop = FALSE]
  if (length(unique(df$group)) < 2) return(NULL)
  
  wt <- wilcox.test(score ~ group, data = df, exact = FALSE)
  
  data.frame(
    pathway = pw,
    median_high = median(df$score[df$group == "High"], na.rm = TRUE),
    median_low  = median(df$score[df$group == "Low"],  na.rm = TRUE),
    delta_median = median(df$score[df$group == "High"], na.rm = TRUE) -
      median(df$score[df$group == "Low"], na.rm = TRUE),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
})

progeny_res <- dplyr::bind_rows(progeny_res)
progeny_res$FDR <- p.adjust(progeny_res$p_value, method = "BH")
progeny_res <- progeny_res %>%
  dplyr::arrange(FDR, dplyr::desc(delta_median))

print(progeny_res)

write.table(
  progeny_res,
  file.path(outdir, "PROGENy_receptorHigh_vs_Low.malignant.stats.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Visualization
plot_df$pathway <- factor(plot_df$pathway, levels = progeny_res$pathway)
plot_df$group <- factor(plot_df$group, levels = c("Low", "High"))

label_df <- plot_df %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    y_max = max(score, na.rm = TRUE),
    y_min = min(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    progeny_res %>% dplyr::select(pathway, p_value, FDR, delta_median),
    by = "pathway"
  ) %>%
  dplyr::mutate(
    y_span = pmax(y_max - y_min, 0.05),
    y_pos = y_max + y_span * 0.18,
    label = paste0(
      "p=", format(p_value, digits = 2, scientific = TRUE)
    )
  )

p_progeny <- ggplot(plot_df, aes(x = pathway, y = score, fill = group)) +
  geom_violin(
    position = position_dodge(width = 0.8),
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = label_df,
    aes(x = pathway, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 3,
    lineheight = 0.9
  ) +
  theme_classic() +
  labs(
    title = "PROGENy pathway activity in malignant cells",
    subtitle = "Low = bottom quartile, High = top quartile of 6-receptor combined score",
    x = NULL,
    y = "PROGENy activity score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_progeny)

ggsave(
  file.path(outdir, "PROGENy_receptorHigh_vs_Low.malignant.violin.pdf"),
  p_progeny,
  width = 6,
  height = 4
)



###### Providing context to Progeny scoring ######
library(dplyr)
library(ggplot2)

# 1. Define receptors and pathways
receptors_raw <- c("BSG", "ADGRE5", "CD44", "P4HB", "PLAUR", "TNFRSF1A")

# If pathways_use already exists, this line is optional
pathways_use <- intersect(
  c("JAK-STAT", "NFkB", "MAPK", "PI3K", "TNFa"),
  colnames(progeny_scores)
)

# 2. Make receptor expression matrix in malignant cells
# This assumes aggregate_symbol_expr_robust() already exists from your previous code.
receptor_res <- aggregate_symbol_expr_robust(
  seu = mal,
  assay = "RNA",
  symbols = receptors_raw
)

receptor_mat <- receptor_res$mat

write.table(
  receptor_res$match_table,
  file.path(outdir, "PROGENy_individualReceptor_matched_features.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(receptor_mat) == 0) {
  stop("No receptor genes matched in malignant RNA assay.")
}

# 3. Match cells between receptor_mat and progeny_scores
common_cells <- intersect(colnames(receptor_mat), rownames(progeny_scores))

receptor_mat_use <- receptor_mat[, common_cells, drop = FALSE]
progeny_use <- progeny_scores[common_cells, pathways_use, drop = FALSE]

cat("Common malignant cells:", length(common_cells), "\n")
cat("Receptors used:\n")
print(rownames(receptor_mat_use))
cat("Pathways used:\n")
print(colnames(progeny_use))

receptor_progeny_list <- list()

for (rc in rownames(receptor_mat_use)) {
  
  vals <- as.numeric(receptor_mat_use[rc, ])
  names(vals) <- colnames(receptor_mat_use)
  
  # Quartile split: Low = bottom 25%, High = top 25%
  q1 <- as.numeric(quantile(vals, 0.25, na.rm = TRUE))
  q3 <- as.numeric(quantile(vals, 0.75, na.rm = TRUE))
  
  keep <- vals <= q1 | vals >= q3
  
  if (sum(keep, na.rm = TRUE) < 10) next
  
  cells_use <- names(vals)[keep]
  group_use <- ifelse(vals[cells_use] >= q3, "High", "Low")
  
  if (length(unique(group_use)) < 2) next
  
  for (pw in colnames(progeny_use)) {
    
    df.tmp <- data.frame(
      cell = cells_use,
      receptor = rc,
      pathway = pw,
      receptor_expr = vals[cells_use],
      group = group_use,
      progeny_score = as.numeric(progeny_use[cells_use, pw]),
      stringsAsFactors = FALSE
    )
    
    df.tmp$group <- factor(df.tmp$group, levels = c("Low", "High"))
    
    wt <- wilcox.test(progeny_score ~ group, data = df.tmp, exact = FALSE)
    
    median_high <- median(df.tmp$progeny_score[df.tmp$group == "High"], na.rm = TRUE)
    median_low  <- median(df.tmp$progeny_score[df.tmp$group == "Low"],  na.rm = TRUE)
    delta_median <- median_high - median_low
    
    receptor_progeny_list[[paste(rc, pw, sep = "_")]] <- data.frame(
      receptor = rc,
      pathway = pw,
      split_method = "quartile_split_Q1Q3",
      q1_expr = q1,
      q3_expr = q3,
      n_low = sum(df.tmp$group == "Low"),
      n_high = sum(df.tmp$group == "High"),
      median_high = median_high,
      median_low = median_low,
      delta_median = delta_median,
      high_gt_low = delta_median > 0,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

receptor_progeny_res <- bind_rows(receptor_progeny_list)

# Multiple-testing correction across all receptor-pathway tests
receptor_progeny_res <- receptor_progeny_res %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(pathway, dplyr::desc(delta_median))

write.table(
  receptor_progeny_res,
  file.path(outdir, "PROGENy_individualReceptor_highLow_byPathway.stats.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(receptor_progeny_res)


# Visualization
denominator <- length(receptors_raw)  # should be 6

pathway_n6_summary <- receptor_progeny_res %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    n_positive = sum(high_gt_low, na.rm = TRUE),
    n_total = denominator,
    count_label = paste0(n_positive, "/", n_total),
    n_tested = dplyr::n_distinct(receptor),
    positive_receptors = paste(receptor[high_gt_low], collapse = ", "),
    negative_or_equal_receptors = paste(receptor[!high_gt_low], collapse = ", "),
    mean_delta_median = mean(delta_median, na.rm = TRUE),
    median_delta_median = median(delta_median, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_positive), dplyr::desc(mean_delta_median))

print(pathway_n6_summary)

write.table(
  pathway_n6_summary,
  file.path(outdir, "PROGENy_pathway_positive_receptor_count_n6.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pathway_effect_summary <- receptor_progeny_res %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    n_receptors = dplyr::n_distinct(receptor),
    n_positive = sum(delta_median > 0, na.rm = TRUE),
    count_label = paste0(n_positive, "/", n_receptors),
    mean_delta_median = mean(delta_median, na.rm = TRUE),
    median_delta_median = median(delta_median, na.rm = TRUE),
    sum_delta_positive = sum(delta_median[delta_median > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_positive),
    dplyr::desc(mean_delta_median)
  )

print(pathway_effect_summary)














