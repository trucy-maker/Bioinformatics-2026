# Pipeline de Análisis de Expresión Diferencial (DEG) y Vías (Pathways)
# Estudio transcriptómico basado en: 
# Fu, Nai Yang, et al. "EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival." Nature Cell Biology 17.4 (2015): 365-375.

# Cargar librerías necesarias
library(edgeR)
library(org.Mm.eg.db)
library(ggplot2)
library(limma)
library(pheatmap)
library(GO.db)

# ==============================================================================
# 1. LECTURA Y PREPARACIÓN DE LA MATRIZ DE CONTEOS
# ==============================================================================
seqdata <- read.delim(file = "data/results/GSE60450_Lactation-GenewiseCounts.txt", row.names = "EntrezGeneID")
seqdata <- seqdata[, 2:ncol(seqdata)]
colnames(seqdata) <- substr(x = colnames(seqdata), start = 1, stop = 7)

# Conversión de identificadores ENTREZID a SYMBOL
symbols <- mapIds(org.Mm.eg.db, keys=rownames(seqdata), column="SYMBOL", keytype="ENTREZID", multiVals="first")

# ==============================================================================
# 2. DEFINICIÓN DE METADATOS Y CREACIÓN DEL OBJETO DGEList
# ==============================================================================
sampleinfo <- read.delim(file = "data/raw/metadata.txt")
group <- paste(sampleinfo$CellType, sampleinfo$Status, sep = ".")
group <- factor(group)

matrizdge <- DGEList(counts = seqdata, group = group)

# Filtrado de genes de baja expresión
keep <- filterByExpr(matrizdge)
matrizdge_filtrado <- matrizdge[keep, , keep.lib.sizes=FALSE]

# ==============================================================================
# 3. NORMALIZACIÓN TMM Y EXPLORACIÓN VISUAL (MDS Plot)
# ==============================================================================
matrizdge_filtrado <- calcNormFactors(matrizdge_filtrado, method="TMM")
logCPM_postnorm <- cpm(matrizdge_filtrado, log=TRUE)

# Gráfico MDS para visualizar similitud entre muestras
colores_grupos <- as.numeric(matrizdge_filtrado$samples$group)
plotMDS(matrizdge_filtrado, 
        col=colores_grupos, 
        pch=16, 
        main="Gráfico MDS")
legend("topleft", 
       legend=levels(matrizdge_filtrado$samples$group), 
       col=1:6, 
       pch=16, 
       cex=0.8)

# ==============================================================================
# 4. MODELADO ESTADÍSTICO - ESTABLECIMIENTO DE CONTRASTES
# ==============================================================================
# Matriz de diseño
design <- model.matrix(~ 0 + matrizdge_filtrado$samples$group)
colnames(design) <- levels(matrizdge_filtrado$samples$group)

# Estimación de dispersiones y ajuste del modelo
matrizdge_filtrado <- estimateDisp(matrizdge_filtrado, design, robust=TRUE)
fit <- glmQLFit(matrizdge_filtrado, design, robust=TRUE)

# Contraste: luminal.lactate vs luminal.pregnant
contraste_lum <- makeContrasts(luminal.lactate - luminal.pregnant, levels=design)
res_luminal <- glmQLFTest(fit, contrast=contraste_lum)
res_corrected_luminal <- topTags(res_luminal, n=Inf)

# ==============================================================================
# 5. IDENTIFICACIÓN DE DEGs Y VISUALIZACIÓN (Volcano Plot)
# ==============================================================================
# Filtrado de DEGs (FDR <= 0.05 y logFC >= 1)
is.de_luminal <- decideTests(res_luminal, adjust.method="BH", p.value=0.05, lfc=1)
summary(is.de_luminal)

# Preparación de datos para Volcano Plot
data <- res_corrected_luminal$table
data$DE <- "NO"
data$DE[data$logFC >= 1 & data$FDR <= 0.05] <- "up"
data$DE[data$logFC <= -1 & data$FDR <= 0.05] <- "down"

degs_solo <- data[data$DE != "NO", ]
genes_top3 <- degs_solo[order(abs(degs_solo$logFC), decreasing = TRUE), ][1:3, ]

ggplot(data, aes(x = logFC, y = -log10(FDR), color = DE)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("NO" = "grey", "up" = "red", "down" = "blue")) +
  geom_text(data = genes_top3, aes(label = rownames(genes_top3)), vjust = -1, color = "black", show.legend = FALSE) +
  labs(title = "Luminal lactate vs pregnant", x = "logFC (Cambio)", y = "-log10 FDR")

# ==============================================================================
# 6. DIAGRAMA DE VENN Y HEATMAP
# ==============================================================================
# Contraste basal (para comparativa)
contraste_bas <- makeContrasts(basal.lactate - basal.pregnant, levels=design)
res_basal <- glmQLFTest(fit, contrast=contraste_bas)
is.de_basal <- decideTests(res_basal, adjust.method="BH", p.value=0.05, lfc=1)

# Diagrama de Venn
vennDiagram(cbind(Luminal = is.de_luminal, Basal = is.de_basal),
            include = c("up", "down"),
            counts.col = c("red", "blue"),
            circle.col = c("pink", "lightblue"),
            cex = c(0.7, 0.4, 1),
            main = "comparativa: luminal vs basal")

# Extracción de genes sobreexpresados compartidos para Heatmap
genesup_luminal <- rownames(is.de_luminal)[is.de_luminal == 1]
genesup_basal <- rownames(is.de_basal)[is.de_basal == 1]
genesup_compartidos <- intersect(genesup_luminal, genesup_basal)

matriz_genesup_compartidos <- logCPM_postnorm[genesup_compartidos, ]
pheatmap(matriz_genesup_compartidos, scale = "row", show_rownames = FALSE, main = "genes sobreexpresados compartidos")

# ==============================================================================
# 7. ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL (GO)
# ==============================================================================
genes_total <- rownames(logCPM_postnorm)
analisis_go <- goana(de = genesup_compartidos, universe = genes_total, species = "Mm")

# Top 10 Procesos Biológicos (BP) y Funciones Moleculares (MF)
top10_BP <- topGO(analisis_go, ontology = "BP", number = 10)
top10_MF <- topGO(analisis_go, ontology = "MF", number = 10)

print(top10_BP)
print(top10_MF)
