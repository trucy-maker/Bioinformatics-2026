#!/bin/bash
# Pipeline de preprocesamiento de RNA-seq (From reads to gene counts)
# Basado en: 
# Fu, Nai Yang, et al. "EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival." Nature Cell Biology 17.4 (2015): 365-375.

# 1. Control de calidad inicial con FastQC
fastqc data/raw/SRR1552444.fastq.gz -o data/processed/01_quality_control/

# 2. Recorte de adaptadores y baja calidad con TrimGalore
trim_galore data/raw/SRR1552444.fastq.gz -o data/processed/02_trimming/

# 3. Control de calidad post-filtrado
fastqc data/processed/02_trimming/SRR1552444_trimmed.fq.gz -o data/processed/01_quality_control/

# 4. Alineamiento con Hisat2 sobre el genoma de referencia de ratón (mm10)
# Se utiliza -k 1 para reportar un único alineamiento por lectura y evitar ambigüedad estadística
hisat2 -x data/raw/mm10/genome -U data/processed/02_trimming/SRR1552444_trimmed.fq.gz -S data/processed/03_alignment/SRR1552444_hisat2_k1.sam -k 1

# Alineamiento alternativo permitiendo hasta 5 mapeos por lectura para comparación
hisat2 -x data/raw/mm10/genome -U data/processed/02_trimming/SRR1552444_trimmed.fq.gz -S data/processed/03_alignment/SRR1552444_hisat2_k5.sam -k 5

# 5. Conversión, ordenación e indexado con SAMtools
# Se convierte el archivo SAM a BAM ordenado por coordenadas para optimizar el espacio y análisis
samtools sort -o data/processed/03_alignment/SRR1552444_hisat2_k1.sorted.bam data/processed/03_alignment/SRR1552444_hisat2_k1.sam
samtools index data/processed/03_alignment/SRR1552444_hisat2_k1.sorted.bam

# Comprobación rápida de FLAGS y valores MAPQ en el archivo SAM
samtools view data/processed/03_alignment/SRR1552444_hisat2_k1.sam | cut -f 2 | sort | uniq -c
samtools view data/processed/03_alignment/SRR1552444_hisat2_k1.sam | cut -f 5 | sort | uniq -c

# 6. Evaluación de la calidad del alineamiento con QualiMap
# Se evalúa la distribución genómica (exones, intrones) y la cobertura de los transcriptos
qualimap rnaseq -bam data/processed/03_alignment/SRR1552444_hisat2_k1.sorted.bam -gtf data/annotation/gencode.vM10.annotation.gtf -java-mem-size=4000M -outdir data/Results/qualimap_results

# 7. Exploración del archivo de anotaciones GTF
# Análisis de las bases de datos origen, features y distribución de genes por cromosoma
grep -v "^#" data/annotation/gencode.vM10.annotation.gtf | cut -f 2 | sort | uniq -c
grep -v "^#" data/annotation/gencode.vM10.annotation.gtf | cut -f 3 | sort | uniq -c
awk -F '\t' '$3=="gene" {print $1}' data/annotation/gencode.vM10.annotation.gtf | sort | uniq -c | sort -n
# Búsqueda del gen más largo
awk -F '\t' '$3=="gene" {print $5-$4, $9}' data/annotation/gencode.vM10.annotation.gtf | sort -nr | head -n 1

# 8. Generación de la matriz de recuentos génicos con HTSeq
# Se cuentan las lecturas que mapean en los exones
htseq-count -t exon -i gene_id -f bam -r pos -s no data/processed/03_alignment/SRR1552444_hisat2_k1.sorted.bam data/annotation/gencode.vM10.annotation.gtf > data/processed/03_alignment/SRR1552444_counts.tsv

# Comprobación de los recuentos finales
tail data/processed/03_alignment/SRR1552444_counts.tsv
grep -v "_" data/processed/03_alignment/SRR1552444_counts.tsv | awk '{sum+=$2} END {print sum}'
