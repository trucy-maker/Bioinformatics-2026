# Análisis RNA-seq

Este proyecto contiene un pipeline completo de análisis de datos de ARN-seq para evaluar los perfiles de expresión génica entre poblaciones de células basales y luminales en glándulas mamarias de ratones en diferentes estadios fisiológicos (vírgenes, gestantes y lactantes).
Los datos de partida se basan en el estudio de referencia: *Fu, Nai Yang, et al. "EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival." Nature Cell Biology 17.4 (2015): 365-375.*

## Estructura

El análisis se divide en dos fases principales

### 1. Preprocesamiento (Bash) - `01_preprocessing.sh`
Flujo de trabajo para pasar de lecturas crudas (FASTQ) a matriz de conteos:
- **Control de calidad:** FastQC y eliminación de adaptadores con TrimGalore.
- **Alineamiento:** Mapeo contra el genoma de referencia de ratón (mm10) utilizando Hisat2.
- **Procesamiento de alineamientos:** Conversión y ordenación con SAMtools.
- **Evaluación post-alineamiento:** QualiMap para analizar la distribución genómica.
- **Recuento génico:** Cuantificación sobre exones utilizando HTSeq-count y anotaciones de GENCODE.

### 2. Expresión Diferencial y Enriquecimiento (R) - `02_DEG_and_Pathways.R`
Análisis estadístico de la matriz de conteos:
- **Filtrado y Normalización:** Uso de `edgeR` para normalización TMM y corrección de sesgos de composición (ej. sobreexpresión masiva de beta-caseína en lactantes).
- **Análisis de Expresión Diferencial (DEG):** Modelado GLM y contrastes estadísticos usando `limma`. Identificación de genes clave regulados al alza y a la baja (FDR <= 0.05, |logFC| >= 1).
- **Visualización de Datos:** Creación de gráficos MDS, Volcano Plots y Heatmaps (`ggplot2`, `pheatmap`) para ilustrar la segregación transcripcional por linaje celular.
- **Ontología Génica (GO):** Análisis de enriquecimiento funcional de procesos biológicos y funciones moleculares.

## 📊 Conclusiones Biológicas Clave
El análisis revela que la lactancia induce una reprogramación metabólica masiva independiente del linaje celular, destacando:
- Sobreexpresión de procesos relacionados con el metabolismo de los lípidos.
- Sobreexpresión de procesos orientados al metabolismo mineral y hormonal de la producción de la leche, y procesos de protección ante los daños oxidativos. Fuerte incremento en la actividad de la aldehído deshidrogenasa (ALDH2) como respuesta protectora ante el estrés oxidativo.

## 📄 Documentación Completa PDF
- `Actividad_1.pdf`: corresponde a `01_preprocessing.sh`
- `Actividad_2.pdf`: corresponde a `02_DEG_and_Pathways.R`
