# ============================================================
# VALIDACIÓN EXTERNA EN GSE39582 — versión corregida
# ============================================================

library(GEOquery)
library(MCPcounter)
library(caret)
library(glmnet)
library(survival)
library(survminer)
library(dplyr)

# ── PASO 1: CARGAR (ya descargado) ───────────────────────────
cat("Cargando GSE39582...\n")
gse <- getGEO("GSE39582", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]

expr  <- exprs(gse)
pheno <- pData(gse)
feat  <- fData(gse)

# ── PASO 2: EXTRAER VARIABLES CLÍNICAS ──────────────────────
# Los nombres correctos según el output:
pheno$mmr_status <- pheno[["mmr.status:ch1"]]
pheno$os_time    <- as.numeric(pheno[["os.delay (months):ch1"]])  # ya en meses
pheno$os_event   <- as.numeric(pheno[["os.event:ch1"]])
pheno$tnm_stage  <- pheno[["tnm.stage:ch1"]]
pheno$dataset    <- pheno[["dataset:ch1"]]

cat("Distribución MMR status:\n")
print(table(pheno$mmr_status, useNA = "ifany"))

# ── PASO 3: SELECCIONAR COHORTE MSS (pMMR) ──────────────────
# En GSE39582 el estado se llama pMMR/dMMR
mss_idx <- !is.na(pheno$mmr_status) &
           grepl("pMMR|MSS|proficient", pheno$mmr_status, ignore.case = TRUE) &
           pheno$dataset != "Non Tumoral"

cat(sprintf("\nMuestras pMMR/MSS seleccionadas: %d\n", sum(mss_idx)))
print(table(pheno$mmr_status[mss_idx]))

expr_mss  <- expr[, mss_idx]
pheno_mss <- pheno[mss_idx, ]

# ── PASO 4: MAPEAR SONDAS A SÍMBOLOS HUGO ───────────────────
gene_symbols <- feat[["Gene symbol"]]
valid <- !is.na(gene_symbols) & gene_symbols != "" &
         !grepl("///", gene_symbols)

expr_v <- expr_mss[valid, ]
syms_v <- gene_symbols[valid]

# Colapsar por mayor varianza
vars_probe <- apply(expr_v, 1, var)
df_tmp     <- data.frame(gene = syms_v, var = vars_probe,
                          stringsAsFactors = FALSE)
df_tmp     <- df_tmp[order(df_tmp$var, decreasing = TRUE), ]
keep_rows  <- !duplicated(df_tmp$gene)
expr_hugo  <- expr_v[rownames(df_tmp)[keep_rows], ]
rownames(expr_hugo) <- df_tmp$gene[keep_rows]

cat(sprintf("Genes únicos: %d\n", nrow(expr_hugo)))

# ── PASO 5: MCP-COUNTER ─────────────────────────────────────
cat("Aplicando MCP-counter...\n")
mcp_val <- MCPcounter.estimate(
  expression   = expr_hugo,
  featuresType = "HUGO_symbols"
)
cat(sprintf("MCP-counter: %d poblaciones x %d muestras\n",
            nrow(mcp_val), ncol(mcp_val)))

saveRDS(mcp_val, "data_processed/mcp_counter_GSE39582_MSS.rds")

# ── PASO 6: PROYECTAR SUBTIPOS CON ELASTIC NET ──────────────
cat("Proyectando subtipos con modelo Elastic Net...\n")

model_enet   <- readRDS("data_processed/model_elastic_net.rds")
logCPM_tcga  <- readRDS("data_processed/tcga_crc_mss_logCPM_HUGO_symbols.rds")

# Genes top500 del entrenamiento
vars_tcga <- apply(logCPM_tcga, 1, var)
top500     <- names(sort(vars_tcga, decreasing = TRUE))[1:500]

genes_comun   <- intersect(top500, rownames(expr_hugo))
genes_missing <- setdiff(top500, genes_comun)
cat(sprintf("Genes en común: %d | Ausentes (se imputan a 0): %d\n",
            length(genes_comun), length(genes_missing)))

# Construir matriz de validación con exactamente los 500 genes en orden
X_val <- matrix(0, nrow = ncol(expr_hugo), ncol = length(top500),
                dimnames = list(colnames(expr_hugo), top500))
X_val[, genes_comun] <- t(expr_hugo[genes_comun, ])
X_val <- as.data.frame(X_val)

# Hacer nombres válidos para R (igual que en entrenamiento)
colnames(X_val) <- make.names(colnames(X_val))

pred_subtypes <- predict(model_enet, newdata = X_val)
pred_prob     <- predict(model_enet, newdata = X_val, type = "prob")

cat("\nDistribución subtipos predichos en GSE39582 MSS:\n")
print(table(pred_subtypes))

pheno_mss$subtipo_pred <- as.character(pred_subtypes)
saveRDS(pheno_mss, "data_processed/pheno_GSE39582_MSS_subtypes.rds")

# ── PASO 7: KAPLAN-MEIER EN COHORTE DE VALIDACIÓN ───────────
surv_val <- pheno_mss %>%
  filter(!is.na(os_time), !is.na(os_event),
         !is.na(subtipo_pred), os_time > 0)

cat(sprintf("\nMuestras con supervivencia completa: %d\n", nrow(surv_val)))
cat("Eventos (muertes):", sum(surv_val$os_event, na.rm = TRUE), "\n")
print(table(surv_val$subtipo_pred))

# Etiquetas limpias para la figura
surv_val$subtipo_label <- recode(surv_val$subtipo_pred,
  "Immune.desert"           = "Immune desert",
  "Immune.intermediate"     = "Immune intermediate",
  "Immune.enriched.subtype" = "Immune-enriched subtype"
)

colores <- c(
  "Immune desert"           = "#377EB8",
  "Immune intermediate"     = "#999999",
  "Immune-enriched subtype" = "#E41A1C"
)

km_val <- survfit(
  Surv(os_time, os_event) ~ subtipo_label,
  data = surv_val
)

p_km_val <- ggsurvplot(
  km_val,
  data              = surv_val,
  palette           = colores,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  xlab              = "Tiempo (meses)",
  ylab              = "Supervivencia global",
  title             = "Validación externa: supervivencia por subtipo inmunológico (GSE39582, pMMR/MSS)",
  legend.title      = "Subtipo predicho",
  ggtheme           = theme_minimal(base_size = 13)
)

ggsave(
  filename = "figures/KM_OS_GSE39582_validation.png",
  plot     = print(p_km_val),
  width = 8, height = 7, dpi = 300
)

# Log-rank y Cox
lr_val   <- survdiff(Surv(os_time, os_event) ~ subtipo_label, data = surv_val)
p_lr_val <- 1 - pchisq(lr_val$chisq, df = length(lr_val$n) - 1)
cat(sprintf("\nLog-rank GSE39582 — p-valor: %.4f\n", p_lr_val))

cat("\nMedianas supervivencia (meses):\n")
print(surv_median(km_val))

surv_val$subtipo_label <- relevel(factor(surv_val$subtipo_label),
                                   ref = "Immune desert")
cox_val <- coxph(Surv(os_time, os_event) ~ subtipo_label, data = surv_val)
cat("\nCox univariante GSE39582:\n")
print(summary(cox_val))

# ── PASO 8: VALIDAR COHERENCIA MCP-COUNTER ──────────────────
mcp_df        <- as.data.frame(t(mcp_val))
mcp_df$sample <- rownames(mcp_df)
mcp_df$subtipo <- pheno_mss$subtipo_pred[match(mcp_df$sample, rownames(pheno_mss))]

cat("\nMedianas MCP-counter por subtipo predicho (GSE39582):\n")
print(aggregate(. ~ subtipo,
                data = mcp_df[, !colnames(mcp_df) %in% "sample"],
                FUN  = median))

cat("\n✓ Validación completada.\n")
cat("Figura guardada en: figures/KM_OS_GSE39582_validation.png\n")
