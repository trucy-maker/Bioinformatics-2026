# Caracterización de la desregulación transcriptómica en la epilepsia farmacorresistente a partir de muestras de sEEG

Este proyecto conforma el análisis de mi TFM en Bioinformática. El objetivo principal es caracterizar el perfil transcriptómico de la zona epileptogénica (Onset) frente a tejido cerebral sano (Cold) en pacientes con epilepsia farmacorresistente (EFR) a partir de muestras procedentes de sEEG.
Los datos proceden de un análisis secundario del dataset **GEO GSE268714** (Dwivedi et al.)

## Metodología
El pipeline completo se ha desarrollado en **R** (`01_DEA_and_Enrichment.R`):
- **Preprocesamiento y DEA:** Uso de `DESeq2` aplicando un diseño experimental intra-sujeto (`~ Paciente + Condicion`) para suprimir la alta variabilidad biológica interindividual, aplicando LFC Shrinkage (estimador 'normal').
- **Visualización Analítica:** Análisis de Componentes Principales (PCA) y visualización de la expresión diferencial mediante Volcano Plots (`ggplot2`, `ggrepel`) y Heatmaps con agrupamiento jerárquico (`pheatmap`).
- **Ontología Génica:** Análisis de enriquecimiento funcional mediante `limma/goana` y la base de datos `org.Hs.eg.db` para mapeo a identificadores ENTREZID, identificando procesos biológicos (BP) y componentes celulares (CC) alterados.

## Resultados
El análisis identificó una firma clave de **10 genes biomarcadores** capaces de segregar clínicamente el tejido epiléptico del sano, revelando dos importantes procesos fisiopatológicos:
1. **Pérdida de la inhibición sináptica:** Infraexpresión significativa de elementos clave en las dendritas y membranas postsinápticas (PLEKHH2, HIPK4), lo que sugiere una incapacidad física para anclar receptores GABA-A. Esto justificaría la falta de respuesta a los fármacos antiepilépticos convencionales.
2. **Estrés celular y metabólico:** Sobreexpresión de procesos de catabolismo de ácidos grasos en peroxisomas (genes NUDT19, LINC02725, ERFE) como mecanismo de respuesta al estrés oxidativo derivado de las convulsiones.

## Documentación
La justificación clínica, revisión de la literatura, metodología detallada y discusión biológica de los resultados se encuentran en la memoria completa del proyecto: TFM_FINAL.pdf
