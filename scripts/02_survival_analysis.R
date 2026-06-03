# ============================================================
# ANÁLISIS DE SUPERVIVENCIA POR SUBTIPO INMUNOLÓGICO - CCR MSS
# ============================================================
# Requiere: survival, survminer, dplyr
# Ejecutar desde el directorio raíz del proyecto TFM
# ============================================================

library(survival)
library(survminer)
library(dplyr)

# ── PASO 1: CARGAR LOS SUBTIPOS YA DEFINIDOS ────────────────
# Carga el objeto con los clusters que generaste en el pipeline
mcp_clustered <- readRDS("data_processed/mcp_counter_clusters_annotated.rds")
subtypes <- data.frame(
  barcode     = rownames(mcp_clustered),
  subtipo     = mcp_clustered$cluster_bio,
  stringsAsFactors = FALSE
)

# Estandarizar barcode a 12 caracteres (TCGA-XX-XXXX)
subtypes$patient_id <- substr(subtypes$barcode, 1, 12)

# ── PASO 2: CARGAR DATOS CLÍNICOS ───────────────────────────
# Carga el SummarizedExperiment con metadatos clínicos
se_mss <- readRDS("data_processed/tcga_crc_mss_primary_unique_se.rds")
clinical <- as.data.frame(colData(se_mss))
clinical$patient_id <- substr(rownames(clinical), 1, 12)

# Ver qué variables de supervivencia hay disponibles
surv_vars <- c("days_to_death", "days_to_last_follow_up",
               "vital_status", "paper_OS", "paper_OS.time",
               "paper_PFI", "paper_PFI.time",
               "days_to_last_followup")

available <- intersect(surv_vars, colnames(clinical))
cat("Variables de supervivencia disponibles:\n")
print(available)
cat("\nNombres de todas las columnas clínicas disponibles:\n")
print(colnames(clinical))

# ── PASO 3: CONSTRUIR TABLA DE SUPERVIVENCIA ────────────────
surv_data <- clinical %>%
  select(patient_id, any_of(surv_vars)) %>%
  left_join(subtypes, by = "patient_id")

# Construir tiempo y evento para OS
# Ajusta los nombres según lo que imprima el paso anterior
if (all(c("days_to_death", "days_to_last_follow_up", "vital_status") %in% colnames(surv_data))) {
  surv_data <- surv_data %>%
    mutate(
      OS_time = ifelse(vital_status == "Dead",
                       as.numeric(days_to_death),
                       as.numeric(days_to_last_follow_up)) / 30.44,
      OS_event = ifelse(vital_status == "Dead", 1, 0)
    )
} else if (all(c("paper_OS.time", "paper_OS") %in% colnames(surv_data))) {
  surv_data <- surv_data %>%
    rename(OS_time = paper_OS.time, OS_event = paper_OS) %>%
    mutate(OS_time = as.numeric(OS_time) / 30.44)
}

# Filtrar NAs y tiempos negativos
surv_clean <- surv_data %>%
  filter(!is.na(OS_time), !is.na(OS_event), !is.na(subtipo), OS_time > 0)

cat("\nMuestras con datos de supervivencia completos por subtipo:\n")
print(table(surv_clean$subtipo))

# ── PASO 4: KAPLAN-MEIER (OS) ───────────────────────────────
km_fit <- survfit(
  Surv(OS_time, OS_event) ~ subtipo,
  data = surv_clean
)

colores <- c(
  "Immune desert"           = "#377EB8",
  "Immune intermediate"     = "#999999",
  "Immune-enriched subtype" = "#E41A1C"
)

p_km <- ggsurvplot(
  km_fit,
  data              = surv_clean,
  palette           = colores,
  pval              = TRUE,
  pval.method       = TRUE,
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  xlab              = "Tiempo (meses)",
  ylab              = "Supervivencia global",
  title             = "Supervivencia global por subtipo inmunológico (CCR MSS, TCGA)",
  legend.title      = "Subtipo",
  ggtheme           = theme_minimal(base_size = 13)
)

ggsave(
  filename = "figures/KM_OS_immune_subtypes.png",
  plot     = print(p_km),
  width = 8, height = 7, dpi = 300
)

# ── PASO 5: LOG-RANK Y MEDIANAS ──────────────────────────────
logrank   <- survdiff(Surv(OS_time, OS_event) ~ subtipo, data = surv_clean)
p_logrank <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
cat(sprintf("\nLog-rank test OS — p-valor global: %.4f\n", p_logrank))

cat("\nMedianas de supervivencia (meses):\n")
print(surv_median(km_fit))

# ── PASO 6: COX UNIVARIANTE ──────────────────────────────────
surv_clean$subtipo <- relevel(factor(surv_clean$subtipo), ref = "Immune desert")
cox_uni <- coxph(Surv(OS_time, OS_event) ~ subtipo, data = surv_clean)
cat("\nCox univariante (referencia: Immune desert):\n")
print(summary(cox_uni))

# ── PASO 7: PFI (si disponible) ──────────────────────────────
if (all(c("paper_PFI.time", "paper_PFI") %in% colnames(surv_data))) {

  surv_pfi <- surv_data %>%
    mutate(
      PFI_time  = as.numeric(paper_PFI.time) / 30.44,
      PFI_event = as.numeric(paper_PFI)
    ) %>%
    filter(!is.na(PFI_time), !is.na(PFI_event), !is.na(subtipo), PFI_time > 0)

  km_pfi <- survfit(Surv(PFI_time, PFI_event) ~ subtipo, data = surv_pfi)

  p_pfi <- ggsurvplot(
    km_pfi,
    data              = surv_pfi,
    palette           = colores,
    pval              = TRUE,
    pval.method       = TRUE,
    conf.int          = FALSE,
    risk.table        = TRUE,
    risk.table.height = 0.28,
    xlab              = "Tiempo (meses)",
    ylab              = "Intervalo libre de progresión",
    title             = "PFI por subtipo inmunológico (CCR MSS, TCGA)",
    legend.title      = "Subtipo",
    ggtheme           = theme_minimal(base_size = 13)
  )

  ggsave(
    filename = "figures/KM_PFI_immune_subtypes.png",
    plot     = print(p_pfi),
    width = 8, height = 7, dpi = 300
  )

  logrank_pfi <- survdiff(Surv(PFI_time, PFI_event) ~ subtipo, data = surv_pfi)
  p_pfi_val   <- 1 - pchisq(logrank_pfi$chisq, df = length(logrank_pfi$n) - 1)
  cat(sprintf("\nLog-rank PFI — p-valor global: %.4f\n", p_pfi_val))

  cat("\nMedianas PFI (meses):\n")
  print(surv_median(km_pfi))
}

cat("\n✓ Análisis completado. Figuras guardadas en figures/\n")
