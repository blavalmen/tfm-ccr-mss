# ============================================================
# PIPELINE PRINCIPAL — Estratificación inmunológica CCR MSS
# ============================================================
# TFM Máster Bioinformática y Bioestadística — UOC 2025-2026
# Autora: Blanca Valenzuela Méndez
# Tutora: Andrea Moreno Manuel
# ============================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(MCPcounter)
library(caret)
library(glmnet)
library(randomForest)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggrepel)

# Crear estructura de directorios
dir.create("data_raw",       showWarnings = FALSE, recursive = TRUE)
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results",        showWarnings = FALSE, recursive = TRUE)
dir.create("figures",        showWarnings = FALSE, recursive = TRUE)

# ── PASO 1: DESCARGA TCGA-COAD/READ ─────────────────────────
file_se <- "data_raw/tcga_coad_read_counts_se.rds"

if (!file.exists(file_se)) {
  query_exp <- GDCquery(
    project       = c("TCGA-COAD", "TCGA-READ"),
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query = query_exp)
  se_exp <- GDCprepare(query = query_exp)
  saveRDS(se_exp, file_se)
} else {
  se_exp <- readRDS(file_se)
}

# ── PASO 2: FILTRADO TUMORES PRIMARIOS ───────────────────────
metadata  <- as.data.frame(colData(se_exp))
se_tumor  <- se_exp[, metadata$shortLetterCode == "TP"]

# ── PASO 3: DEDUPLICACIÓN POR PACIENTE ───────────────────────
metadata_tumor             <- as.data.frame(colData(se_tumor))
metadata_tumor$patient_id  <- substr(rownames(metadata_tumor), 1, 12)
metadata_tumor$sample_vial <- substr(rownames(metadata_tumor), 14, 16)
metadata_tumor    <- metadata_tumor[order(metadata_tumor$patient_id,
                                          metadata_tumor$sample_vial), ]
selected_barcodes <- metadata_tumor$barcode[!duplicated(metadata_tumor$patient_id)]
se_tumor_unique   <- se_tumor[, selected_barcodes]

# ── PASO 4: NORMALIZACIÓN TMM + logCPM ───────────────────────
counts <- assay(se_tumor_unique)
dge    <- DGEList(counts = counts)
keep   <- filterByExpr(dge)
dge    <- dge[keep, , keep.lib.sizes = FALSE]
dge    <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)
saveRDS(logCPM, "data_processed/tcga_coad_read_primary_tumor_unique_logCPM.rds")

# ── PASO 5: SELECCIÓN COHORTE MSS ────────────────────────────
metadata_msi <- as.data.frame(colData(se_tumor_unique))
mss_idx      <- metadata_msi$paper_MSI_status == "MSS" &
                !is.na(metadata_msi$paper_MSI_status)
se_mss       <- se_tumor_unique[, mss_idx]
logCPM_mss   <- logCPM[, colnames(se_mss)]
saveRDS(se_mss,     "data_processed/tcga_crc_mss_primary_unique_se.rds")
saveRDS(logCPM_mss, "data_processed/tcga_crc_mss_primary_unique_logCPM.rds")

# ── PASO 6: CONVERSIÓN ENSEMBL → HUGO ────────────────────────
ensembl_ids  <- gsub("\\..*", "", rownames(logCPM_mss))
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                       column = "SYMBOL", keytype = "ENSEMBL",
                       multiVals = "first")
valid_idx         <- !is.na(gene_symbols)
logCPM_mss_hugo   <- logCPM_mss[valid_idx, ]
rownames(logCPM_mss_hugo) <- gene_symbols[valid_idx]
logCPM_mss_hugo   <- logCPM_mss_hugo[!duplicated(rownames(logCPM_mss_hugo)), ]
saveRDS(logCPM_mss_hugo, "data_processed/tcga_crc_mss_logCPM_HUGO_symbols.rds")

# ── PASO 7: DECONVOLUCIÓN MCP-counter ────────────────────────
mcp_results <- MCPcounter.estimate(logCPM_mss_hugo,
                                   featuresType = "HUGO_symbols")
mcp_scores  <- t(mcp_results)
saveRDS(mcp_scores, "data_processed/mcp_counter_scores_mss.rds")

# ── PASO 8: CLUSTERING JERÁRQUICO (Ward.D2, k=3) ─────────────
mcp_scaled  <- scale(mcp_scores)
dist_matrix <- dist(mcp_scaled, method = "euclidean")
hc          <- hclust(dist_matrix, method = "ward.D2")
clusters    <- cutree(hc, k = 3)

# Asignar etiquetas biológicas por score global de infiltración
cluster_df        <- as.data.frame(mcp_scaled)
cluster_df$cluster <- clusters
cluster_summary   <- aggregate(. ~ cluster, data = cluster_df, FUN = mean)
cluster_summary$tme_score <- rowMeans(cluster_summary[, -1])
ordered_clusters  <- cluster_summary$cluster[order(cluster_summary$tme_score)]
labels_map <- setNames(
  c("Immune desert", "Immune intermediate", "Immune-enriched subtype"),
  ordered_clusters
)
mcp_clustered             <- as.data.frame(mcp_scores)
mcp_clustered$cluster_bio <- factor(
  labels_map[as.character(clusters)],
  levels = c("Immune desert", "Immune intermediate", "Immune-enriched subtype")
)
saveRDS(mcp_clustered, "data_processed/mcp_counter_clusters_annotated.rds")

# ── PASO 9: PCA ──────────────────────────────────────────────
cluster_colors <- c(
  "Immune desert"           = "#377EB8",
  "Immune intermediate"     = "#999999",
  "Immune-enriched subtype" = "#E41A1C"
)

pca_mcp <- prcomp(mcp_scaled, center = TRUE, scale. = FALSE)
pca_df  <- data.frame(pca_mcp$x, cluster_bio = mcp_clustered$cluster_bio)
var_exp <- round(100 * summary(pca_mcp)$importance[2, 1:2], 1)

cat("Loadings PC1 y PC2:\n")
print(round(pca_mcp$rotation[, 1:2], 3))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cluster_bio)) +
  geom_point(size = 2.5, alpha = 0.75) +
  stat_ellipse(aes(group = cluster_bio), type = "norm",
               linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = cluster_colors, name = "Subtipo") +
  labs(title = "PCA de perfiles inmunológicos en CCR MSS",
       x = paste0("PC1 (", var_exp[1], "% varianza)"),
       y = paste0("PC2 (", var_exp[2], "% varianza)")) +
  theme_minimal(base_size = 12)

ggsave("figures/pca_immune_subtypes.png", p_pca, width = 7, height = 6, dpi = 300)

# ── PASO 10: HEATMAP ─────────────────────────────────────────
pop_order <- c("Cytotoxic lymphocytes", "CD8 T cells", "T cells",
               "NK cells", "B lineage", "Myeloid dendritic cells",
               "Monocytic lineage", "Neutrophils", "Endothelial cells",
               "Fibroblasts")

col_order <- order(mcp_clustered$cluster_bio)
mcp_ord   <- as.data.frame(mcp_scaled)[col_order, ]

ann_col <- data.frame(
  Subtipo = mcp_clustered$cluster_bio[col_order],
  row.names = rownames(mcp_ord)
)

group_sizes <- table(factor(mcp_clustered$cluster_bio[col_order],
                            levels = c("Immune desert", "Immune intermediate",
                                       "Immune-enriched subtype")))
gaps_col <- cumsum(as.numeric(group_sizes))[-length(group_sizes)]

mat <- t(as.matrix(mcp_ord[, pop_order]))

png("figures/heatmap_tcga_subtypes.png", width = 3000, height = 1500, res = 300)
pheatmap(
  mat,
  annotation_col    = ann_col,
  annotation_colors = list(Subtipo = cluster_colors),
  color             = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  gaps_col          = gaps_col,
  show_colnames     = FALSE,
  border_color      = NA,
  fontsize_row      = 11,
  main              = "Microambiente tumoral por subtipo inmunológico (TCGA, CCR MSS)"
)
dev.off()

# ── PASO 11: BOXPLOTS MCP-counter ────────────────────────────
mcp_long <- mcp_clustered %>%
  select(all_of(pop_order), cluster_bio) %>%
  pivot_longer(cols = -cluster_bio,
               names_to = "Poblacion", values_to = "Score") %>%
  mutate(Poblacion = factor(Poblacion, levels = pop_order))

p_box <- ggplot(mcp_long, aes(x = cluster_bio, y = Score, fill = cluster_bio)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ Poblacion, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = cluster_colors, name = "Subtipo") +
  scale_x_discrete(labels = NULL) +
  labs(title = "Distribución de poblaciones celulares por subtipo (CCR MSS)",
       y = "Score MCP-counter", x = NULL) +
  theme_minimal(base_size = 11)

ggsave("figures/boxplots_mcp_subtypes.png", p_box,
       width = 14, height = 7, dpi = 300)

# ── PASO 12: MODELOS SUPERVISADOS ────────────────────────────
logCPM_hugo <- readRDS("data_processed/tcga_crc_mss_logCPM_HUGO_symbols.rds")
common      <- intersect(colnames(logCPM_hugo), rownames(mcp_clustered))
y   <- make.names(factor(mcp_clustered[common, ]$cluster_bio))
X   <- t(logCPM_hugo[, common])
top500 <- names(sort(apply(X, 2, var), decreasing = TRUE))[1:500]
X   <- X[, top500]

set.seed(123)
ctrl <- trainControl(method = "repeatedcv", number = 5,
                     repeats = 3, classProbs = TRUE)
model_enet <- train(x = X, y = y, method = "glmnet",
                    trControl = ctrl, tuneLength = 10)
model_rf   <- train(x = X, y = y, method = "rf",
                    trControl = ctrl, tuneLength = 5)

saveRDS(model_enet, "data_processed/model_elastic_net.rds")
saveRDS(model_rf,   "data_processed/model_random_forest.rds")

cat("\nRendimiento Elastic Net:\n")
print(model_enet$results[which.max(model_enet$results$Accuracy), ])
cat("\nRendimiento Random Forest:\n")
print(model_rf$results[which.max(model_rf$results$Accuracy), ])

# ── PASO 13: GENES DISCRIMINANTES ────────────────────────────
coef_list  <- coef(model_enet$finalModel, s = model_enet$bestTune$lambda)
genes_list <- lapply(coef_list, function(cm) {
  df      <- as.data.frame(as.matrix(cm))
  df$gene <- rownames(df)
  df$coef <- df[, 1]
  df      <- df[df$coef != 0 & df$gene != "(Intercept)", ]
  df[order(abs(df$coef), decreasing = TRUE), c("gene", "coef")]
})
saveRDS(genes_list, "data_processed/elastic_net_selected_genes_by_subtype.rds")

cat("\nTop genes discriminantes Immune-enriched:\n")
print(head(genes_list[["Immune.enriched.subtype"]], 10))

# ── PASO 14: EXPRESIÓN DIFERENCIAL Y ENRIQUECIMIENTO ─────────
group  <- make.names(factor(mcp_clustered$cluster_bio))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(factor(group))
fit2   <- eBayes(contrasts.fit(
  lmFit(logCPM_mss_hugo, design),
  makeContrasts(Immune.enriched.subtype - Immune.desert, levels = design)
))
deg    <- topTable(fit2, number = Inf, adjust.method = "BH")
write.csv(deg, "results/DEG_immune_enriched_vs_desert.csv")

sig_up <- rownames(deg[deg$adj.P.Val < 0.05 & deg$logFC > 1, ])
cat(paste0("\nGenes sobreexpresados en Immune-enriched (FDR<0.05, logFC>1): ",
           length(sig_up), "\n"))

entrez <- bitr(sig_up, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

ego   <- simplify(enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                            ont = "BP", pAdjustMethod = "BH", readable = TRUE))
ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                    pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

saveRDS(ego,   "data_processed/ego_enriched_vs_desert.rds")
saveRDS(ekegg, "data_processed/ekegg_enriched_vs_desert.rds")

cat("\nTop terminos GO:\n")
print(head(as.data.frame(ego)[, c("Description", "GeneRatio", "p.adjust")], 6))

cat("\nTop rutas KEGG:\n")
print(head(as.data.frame(ekegg)[, c("Description", "GeneRatio", "p.adjust")], 6))

cat("\n✓ Pipeline completado. Continuar con 02_survival_analysis.R\n")
