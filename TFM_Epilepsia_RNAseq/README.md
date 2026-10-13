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

## ⚖️ Licencia y Derechos de los Datos
- **Código y Memoria:** El código desarrollado en R para este pipeline y la memoria del TFM son de mi autoría. El código se distribuye bajo la licencia MIT.
- **Datos Originales:** Los datos crudos analizados en este proyecto son de acceso público y pertenecen al estudio original de: Dwivedi AK, Mahesh A, Sanfeliu A, Larkin J, Siwicki RA, Sweeney KJ, et al. High-resolution multimodal profiling of human epileptic brain activity via explanted depth electrodes. JCI Insight. 2024 [citado 31 marzo 2026]; 10(1):e184518. Disponible en: https://pmc.ncbi.nlm.nih.gov/articles/PMC11721296/
Los datos fueron obtenidos a través de la base de datos Gene Expression Omnibus (GEO) con el identificador **GSE268714**. Este repositorio no redistribuye los datos originales, únicamente contiene el código necesario para reproducir el análisis secundario.
