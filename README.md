# Immunological Stratification of MSS Colorectal Cancer

**TFM — Máster en Bioinformática y Bioestadística (UOC)**  
**Autora:** Blanca Valenzuela Méndez  
**Tutora:** Andrea Moreno Manuel  
**Curso:** 2025-2026

## Descripción

Pipeline bioinformático para la estratificación inmunológica del cáncer colorrectal microsatélite estable (CCR MSS) mediante deconvolución transcriptómica y aprendizaje automático.

Se identificaron tres subtipos inmunológicos — Immune desert, Immune intermediate e Immune-enriched subtype — a partir de datos RNA-seq de TCGA-COAD/READ, con validación biológica externa en GSE39582.

## Estructura del repositorio

```
├── scripts/
│   ├── 01_pipeline_main.R         # Pipeline completo: descarga, normalización,
│   │                              # deconvolución, clustering, modelos supervisados,
│   │                              # expresión diferencial y enriquecimiento funcional
│   ├── 02_survival_analysis.R     # Análisis de supervivencia en TCGA
│   └── 03_validation_GSE39582.R   # Validación externa en GSE39582
├── figures/                       # Figuras generadas (no incluidas en el repo)
├── data_processed/                # Objetos RDS intermedios (no incluidos)
└── README.md
```

## Requisitos

R versión 4.3.1 o superior.

### Paquetes Bioconductor

```r
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "edgeR",
  "limma",
  "org.Hs.eg.db",
  "clusterProfiler"
))
```

### Paquetes CRAN

```r
install.packages(c(
  "MCPcounter",
  "caret",
  "glmnet",
  "randomForest",
  "ggplot2",
  "pheatmap",
  "survival",
  "survminer",
  "GEOquery",
  "dplyr",
  "tidyr",
  "ggrepel",
  "reshape2"
))
```

## Uso

Ejecutar los scripts en orden desde el directorio raíz del proyecto:

```r
source("scripts/01_pipeline_main.R")
source("scripts/02_survival_analysis.R")
source("scripts/03_validation_GSE39582.R")
```

Los scripts crean automáticamente las carpetas `data_raw/`, `data_processed/`, `results/` y `figures/`.

## Datos

- **TCGA-COAD/READ**: descargados automáticamente mediante `TCGAbiolinks`
- **GSE39582**: descargado automáticamente mediante `GEOquery`

No se incluyen datos en el repositorio por restricciones de tamaño. Los objetos intermedios (`.rds`) se generan al ejecutar el pipeline.

## Referencia

Valenzuela Méndez B. Estratificación inmunológica del cáncer colorrectal microsatélite estable mediante deconvolución inmune y aprendizaje automático. TFM, Máster en Bioinformática y Bioestadística, UOC. 2025-2026.
