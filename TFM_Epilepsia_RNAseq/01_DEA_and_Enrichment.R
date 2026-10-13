# =========================================================================================
# TFM: Caracterización de la desregulación transcriptómica en la epilepsia farmacorresistente
# Autor: Sonia María López Marante (estudio secundario)
# Autor estudio original: Dwivedi AK, Mahesh A, Sanfeliu A, Larkin J, Siwicki RA, Sweeney KJ, et al. High-resolution multimodal profiling of human epileptic brain activity via explanted depth electrodes. JCI Insight. 2024 [citado 31 marzo 2026]; 10(1):e184518. 
# Disponible en: https://pmc.ncbi.nlm.nih.gov/articles/PMC11721296/
# Datos: GSE268714 (Dwivedi et al. 2024) - Muestras bulk RNA-seq de electrodos sEEG
# =========================================================================================

# Cargar librerías necesarias
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(org.Hs.eg.db)
library(limma)

# =========================================================================================
# 1. PREPARACIÓN DE DATOS Y METADATOS
# =========================================================================================
# Lectura de la matriz de conteos y metadatos
count_matrix <- read.delim("GSE268714_EPILEPSY.BULK.RNAseq.Count.txt.gz", row.names=1)
metadata <- read.csv("metadatos_pacientes.csv", row.names=1)

# Factorización y establecimiento de referencias (tejido sano "Cold" como referencia)
metadata$Condicion <- as.factor(metadata$Condicion)
metadata$Condicion <- relevel(metadata$Condicion, ref = "Cold")
metadata$Paciente <- as.factor(metadata$Paciente)

# Sincronización de matriz y metadatos
count_matrix_filtrada <- count_matrix[, rownames(metadata)]

# =========================================================================================
# 2. ANÁLISIS DE EXPRESIÓN DIFERENCIAL (DESeq2)
# =========================================================================================
# Creación del objeto DESeqDataSet con diseño intra-sujeto
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtrada,
                              colData = metadata,
                              design = ~ Paciente + Condicion)

# Filtrado de genes con baja expresión (umbral: >= 20 lecturas)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]

# Ejecución del modelo DESeq2
dds <- DESeq(dds)

# Aplicación de LFC Shrinkage (tipo 'normal') para atenuar falsos positivos
res_shrink <- lfcShrink(dds, coef="Condicion_Onset_vs_Cold", type="normal")
summary(res_shrink)

# =========================================================================================
# 3. VISUALIZACIÓN EXPLORATORIA (PCA)
# =========================================================================================
# Transformación estabilizadora de varianza (VST)
vsd <- vst(dds, blind = FALSE)

# Generación del gráfico PCA
pca_plot <- plotPCA(vsd, intgroup = c("Condicion", "Paciente"))
print(pca_plot)

# =========================================================================================
# 4. VISUALIZACIÓN DE EXPRESIÓN DIFERENCIAL (Volcano Plot)
# =========================================================================================
df_shrink <- as.data.frame(res_shrink) %>% drop_na()
df_shrink$Gen <- rownames(df_shrink)

# Clasificación según significancia (FDR < 0.05 y |log2FC| > 1)
df_shrink <- df_shrink %>%
  mutate(Categoria = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Sobreexpresado (Up)",
    padj < 0.05 & log2FoldChange < -1 ~ "Infraexpresado (Down)",
    TRUE ~ "No significativo"
  ))

# Top 10 genes significativos (excluyendo "LOC")
top_genes <- df_shrink %>%
  filter(Categoria != "No significativo", !grepl("^LOC", Gen)) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(10)

# Renderizado del Volcano Plot
plot_volcano_final <- ggplot(df_shrink, aes(x = log2FoldChange, y = -log10(padj), color = Categoria)) +
  geom_point(alpha = 0.8, size = 1) +
  scale_color_manual(values = c("Sobreexpresado (Up)" = "red", 
                                "Infraexpresado (Down)" = "navy", 
                                "No significativo" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 5)) + 
  geom_text_repel(data = top_genes, aes(label = Gen), size = 3.5, color = "black", 
                  box.padding = 0.5, max.overlaps = Inf, fontface = "italic") + 
  labs(title = "Expresión diferencial: Onset vs Cold", x = "Log2 Fold Change", y = "-Log10(P-valor ajustado)") +
  theme_minimal() + theme(legend.position = "bottom")

print(plot_volcano_final)

# Guardar resultados de DEG
genes_finales <- df_shrink %>% filter(Categoria != "No significativo")
write.csv(genes_finales, "Genes_Diferenciales.csv")

# =========================================================================================
# 5. PERFILES DE EXPRESIÓN (Heatmap)
# =========================================================================================
matriz_norm <- assay(vsd)
genes_heatmap <- top_genes$Gen
matriz_filtrada <- matriz_norm[genes_heatmap, ]

anotacion_columnas <- data.frame(Condicion = metadata$Condicion)
rownames(anotacion_columnas) <- rownames(metadata)
colores_anotacion <- list(Condicion = c("Onset" = "firebrick", "Cold" = "navy"))

pheatmap(matriz_filtrada, annotation_col = anotacion_columnas, annotation_colors = colores_anotacion,
         scale = "row", show_colnames = TRUE, treeheight_row = 0, fontsize_row = 8, 
         main = "Perfil de expresión de los 10 genes más significativos")

# =========================================================================================
# 6. ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL (GO)
# =========================================================================================
# Definición del universo génico (SYMBOL a ENTREZID)
genes_universo <- df_shrink %>% filter(!grepl("^LOC", Gen)) %>% pull(Gen)
universo_entrezid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = genes_universo, 
                                                 column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")))

# Genes Up y Down (Umbral ajustado: |log2FC| > 0.5 para mayor sensibilidad en pathways)
genes_up <- df_shrink %>% filter(padj < 0.05 & log2FoldChange > 0.5, !grepl("^LOC", Gen)) %>% pull(Gen)
genes_down <- df_shrink %>% filter(padj < 0.05 & log2FoldChange < -0.5, !grepl("^LOC", Gen)) %>% pull(Gen)

up_entrezid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = genes_up, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")))
down_entrezid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = genes_down, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")))

lista_genes <- list(UP = up_entrezid, DOWN = down_entrezid)

# Ejecución de Goana
res_goana <- goana(de = lista_genes, universe = universo_entrezid, species = "Hs")

# Extracción de Top 10 resultados (BP y CC)
top_BP_up <- topGO(res_goana, ontology = "BP", sort = "UP", number = 10)
top_CC_up <- topGO(res_goana, ontology = "CC", sort = "UP", number = 10)
top_BP_down <- topGO(res_goana, ontology = "BP", sort = "DOWN", number = 10)
top_CC_down <- topGO(res_goana, ontology = "CC", sort = "DOWN", number = 10)

# Guardar resultados de GO
write.csv(top_BP_up, "GO_ProcesosBio_UP.csv")
write.csv(top_CC_up, "GO_ComponenteCel_UP.csv")
write.csv(top_BP_down, "GO_ProcesosBio_DOWN.csv")
write.csv(top_CC_down, "GO_ComponenteCel_DOWN.csv")
