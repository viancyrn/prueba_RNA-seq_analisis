# Bulk RNA-seq Analysis Pipeline

**Dataset:** GEO GSE201955  
**Organism:** *Homo sapiens*  
**Samples:** SRR19020908, SRR19020909, SRR19020910  
**Sequencing:** Single-end RNA-seq, 51 bp  
**Reference genome:** GRCh38  

**Languages:** [English](#english) | [Español](#español)

---

# English 🇬🇧 

## 1. Project overview

This repository is part of my training to create a reproducible bulk RNA-seq analysis workflow using publicly available sequencing data from GEO accession **GSE201955**, which includes RNA-seq data from human bronchial epithelial cells.
As of now, I have only selected 3 samples to test the workflow and understand what happens at each step.
The workflow covers the initial stages of RNA-seq analysis, from raw sequencing reads to gene-level count, not with the sole goal to run the tools, but also check the results, understand possible problems, and document the decisions I make along the way.

So far, I have completed:

1. Downloading the sequencing data.
2. Converting SRA files to FASTQ.
3. Checking raw read quality with FastQC.
4. Processing reads with fastp.
5. Comparing the original and optimized trimming results.
6. Aligning the reads to the human reference genome.
7. Checking the BAM files and alignment statistics.
8. Counting reads assigned to genes.
9. Combining the final quality-control results in MultiQC.

The next step will be to prepare the count matrix and continue the analysis in R.

**Differential gene expression analysis has not yet been performed.**

## 2. Dataset and samples

I am working with the following three samples:

| Sample | Sequencing type | Original read length |
|---|---|---|
| SRR19020908 | Single-end | 51 bp |
| SRR19020909 | Single-end | 51 bp |
| SRR19020910 | Single-end | 51 bp |

The sample IDs are listed in [`samples.txt`](samples.txt).

For the alignment, I used the human reference genome **GRCh38**.

For gene counting, I used the Ensembl GRCh38 release 115 annotation.

This is a small test dataset. Before moving to biological comparisons, I still need to check the sample metadata and experimental design.

## 3. Repository structure

The main files and folders in this project are:

```text
prueba_RNA-seq_analisis/
├── README.md
├── .gitignore
├── samples.txt
├── scripts/
│   ├── 01_download.sh
│   ├── 02_fastq.sh
│   ├── 03_fastp.sh
│   ├── 04_align.sh
│   └── 05_featurecounts.sh
├── docs/
│   └── fastqc_trimming_analysis.md
├── 00_raw/
├── 01_clean/
├── 02_bam/
├── 03_counts/
├── 04_multiqc/
│   ├── multiqc_config.yaml
│   └── multiqc_final_tail1.html
└── ref/
```

Some folders contain large files that are only stored locally.

I do not upload raw FASTQ files, BAM files, reference files, or count tables to GitHub. Instead, I keep the scripts, documentation, and selected reports needed to understand the analysis.

## 4. Tools and workflow

I ran the analysis in a Conda environment called `rnaseq_qc`.

These are the main tools I used:

| Tool | What I used it for |
|---|---|
| SRA Toolkit | Downloading and converting sequencing data |
| FastQC | Checking read quality |
| fastp | Filtering and trimming reads |
| Subread | Aligning reads to GRCh38 |
| SAMtools | Checking BAM files and alignment statistics |
| featureCounts | Counting reads assigned to genes |
| MultiQC | Combining results into QC reports |

### Workflow

```text
SRA download
     |
     v
FASTQ conversion
     |
     v
FastQC on raw reads
     |
     v
fastp trimming
     |
     v
FastQC after trimming
     |
     v
Alignment to GRCh38
     |
     v
BAM validation + flagstat
     |
     v
featureCounts
     |
     v
Final MultiQC report
     |
     v
Next: analysis in R
```

The scripts are available in the [`scripts/`](scripts/) folder.

They should be run from the main project directory.

## 5. Quality control and trimming

### 5.1 What I found in the raw reads

I first checked the raw FASTQ files with FastQC.

The main issues were:

- Lower sequencing quality at the last base (position 51).
- Higher `N` content at the last base.
- Uneven nucleotide composition at the beginning of the reads, mainly along the first 13 positions.

The problem at the last base was consistent across the three samples, so I decided to test whether removing that base would improve the results.

I did not remove the first 13 bases because the nucleotide composition pattern could be related to the RNA-seq library preparation. I did not want to trim those positions without a clear reason.

### 5.2 Original and optimized trimming

I first processed the reads using fastp using it's standard settings and saved the output as:

```text
_clean.fastq.gz
```

After reviewing the FastQC results, I tested a second approach using additional fastp settings:

```bash
--trim_tail1 1
```

This removes one base from the 3′ end of each read.

The optimized output files are named:

```text
_tail1.fastq.gz
```

I kept both versions so I could compare the results instead of replacing the original files.

### 5.3 How many reads were retained?

| Sample | Raw reads | Reads after optimized fastp | Retained |
|---|---:|---:|---:|
| SRR19020908 | 3,075,016 | 3,055,200 | 99.36% |
| SRR19020909 | 3,562,425 | 3,543,079 | 99.46% |
| SRR19020910 | 3,629,875 | 3,615,115 | 99.59% |

The Q30 percentages after filtering were approximately 99% for all three samples.

The optimized strategy also retained slightly more reads than the original fastp processing:

| Sample | Additional reads retained |
|---|---:|
| SRR19020908 | +929 |
| SRR19020909 | +1,049 |
| SRR19020910 | +1,100 |

The differences are small, but the optimized strategy gave better results for the problem I was trying to correct.

### 5.4 FastQC after optimized trimming

After trimming the last base, I ran FastQC again.

| FastQC module | Result |
|---|---|
| Per base sequence quality | PASS |
| Per sequence quality scores | PASS |
| Per base N content | PASS |
| Per base sequence content | FAIL |
| Per sequence GC content | PASS |
| Adapter Content | PASS |
| Overrepresented sequences | PASS |
| Sequence Length Distribution | WARNING |

The problems with quality and `N` content at the last base were resolved.

The **Per base sequence content** module still failed.

I decided not to apply additional trimming just to make every FastQC module pass. I think it is more important to understand why a warning appears than to remove data without knowing whether that is necessary.

The full analysis between the different fastp results are documented here:

[FastQC and trimming analysis](docs/fastqc_trimming_analysis.md)

## 6. Alignment to the reference genome

After selecting the optimized reads, I aligned them to the human reference genome GRCh38 using Subread.

The script is:

[`scripts/04_align.sh`](scripts/04_align.sh)

The new BAM files were saved with the `_tail1.bam` suffix so the original BAM files would not be overwritten so I could still use them for later comparisons.

### Alignment results

| Sample | Input reads | Mapped reads | Mapping rate |
|---|---:|---:|---:|
| SRR19020908 | 3,055,200 | 2,714,036 | 88.83% |
| SRR19020909 | 3,543,079 | 3,156,878 | 89.10% |
| SRR19020910 | 3,615,115 | 3,228,911 | 89.32% |

The mapping rates were similar across all three samples.

I also checked the BAM files using:

```bash
samtools quickcheck -v 02_bam/*_tail1.bam
```

All three files passed the check.

After that, I generated SAMtools `flagstat` reports to include the mapping statistics in MultiQC.

## 7. Gene counting with featureCounts

Once the BAM files were ready, I used featureCounts to count reads assigned to genes.

The script is:

[`scripts/05_featurecounts.sh`](scripts/05_featurecounts.sh)

I used the following annotation:

```text
ref/Homo_sapiens.GRCh38.115.gtf.gz
```

The main parameters were:

```bash
-T 2
-t exon
-g gene_id
-s 0
```

I compared the different strandedness settings in featureCounts to see which one assigned the most reads to genes.
The unstranded setting (`-s 0`) gave the highest number of assigned reads, so I decided to use it for the final gene counts. This result supports the choice of unstranded counting, although the number of assigned reads alone does not definitively confirm the library preparation method.

For other RNA-seq datasets, strandedness information may already be available in the study documentation. It is useful to check the GEO or SRA metadata, the original publication, and especially the library preparation methods or sequencing protocol.

If this information is not available or is unclear, strandedness can be investigated by comparing featureCounts results using `-s 0`, `-s 1`, and `-s 2`, or by using a tool such as RSeQC (`infer_experiment.py`).

It is important to check this before gene counting because using the wrong strandedness setting can affect the number of reads assigned to genes.

The optimized results were saved in:

```text
03_counts/gene_counts_tail1.txt
```

### Gene assignment results

### Why I repeated the analysis

I initially processed the reads using fastp and continued with alignment and gene counting. At that point, I had not yet performed a separate FastQC analysis.
After going over my workflow with my professor, he pointed out that I had skipped FastQC. He explained that although fastp provides quality statistics and reports, I should also examine the sequencing data independently before deciding whether additional trimming is necessary.

I therefore went back to the raw FASTQ files and ran FastQC. This helped me identify the lower quality and increased N content at the final base, as well as the nucleotide composition bias at the beginning of the reads.
Based on these results, I decided to test an optimized fastp strategy using `--trim_tail1 1`. I then repeated FastQC, alignment, and gene counting so I could compare both workflows.

Using FastQC together with fastp was useful because the tools provide complementary information. FastQC gives clear figures showing quality patterns across read positions, nucleotide composition, GC content, and other potential problems. These figures made it easier for me to identify what needed attention and to evaluate the results of fastp more carefully.

The comparison also helped me understand that a trimming decision should be based on the quality of the data, not only on whether a tool finishes successfully.

### Gene assignment results

| Sample | Assigned reads | Assignment rate |
|---|---:|---:|
| SRR19020908 | 2,093,666 | 68.5% |
| SRR19020909 | 2,563,479 | 72.4% |
| SRR19020910 | 2,572,873 | 71.2% |

I also compared the gene counts from the initial and optimized workflows.

| Sample | Original assigned reads | Optimized assigned reads | Difference |
|---|---:|---:|---:|
| SRR19020908 | 2,093,495 | 2,093,666 | +171 |
| SRR19020909 | 2,563,077 | 2,563,479 | +402 |
| SRR19020910 | 2,572,872 | 2,572,873 | +1 |

The differences were very small.

This suggests that trimming the final base improved the read quality without substantially changing the gene counting results.

## 8. MultiQC reports

I generated separate MultiQC reports because I wanted to keep the different processing stages clear.

The reports include:

1. `multiqc_before.html` — FastQC before trimming.
2. `multiqc_after.html` — FastQC after the original trimming.
3. `multiqc_tail1.html` — FastQC after optimized trimming.
4. `multiqc_fastp_comparison.html` — comparison of both fastp strategies.
5. `multiqc_final_tail1.html` — final report for the optimized pipeline.

The intermediate reports are stored locally.

### Final MultiQC report

For the final report, I combined:

- FastQC results from the optimized reads.
- fastp results from the optimized processing.
- SAMtools `flagstat` statistics.
- featureCounts assignment statistics.

At first, MultiQC displayed multiple rows for the same sample because different tools used slightly different sample names.

This was solved by using a configuration file to standardize the sample names.

Now the **General Statistics** table shows three rows, one for each sample, with results from the different tools combined.

**Final report:**

[Download the final MultiQC report](04_multiqc/multiqc_final_tail1.html)

**Configuration file:**

[View the MultiQC configuration](04_multiqc/multiqc_config.yaml)

**Note:** GitHub may not preview the HTML file because of its size. To view the interactive report, download the raw HTML file and open it in a browser.

## 9. What I learned from the results

The optimized pipeline worked for all three samples.

The main results were:

- More than 99% of the reads were retained after filtering.
- Q30 was approximately 99% for all samples.
- Removing the final base corrected the quality and `N` content issues.
- Mapping rates were between 88.83% and 89.32%.
- Gene assignment rates were between 68.5% and 72.4%.
- All optimized BAM files passed validation.
- The optimized and original gene counts were very similar.

## 10. Next steps

The next part of the project will focus on working with the gene counts in R.

I plan to:

1. Import the featureCounts output.
2. Prepare the count matrix.
3. Review the sample metadata and experimental groups.
4. Confirm library strandedness.
5. Explore the data using PCA and sample clustering.
6. Check whether the experimental design supports differential expression analysis with DESeq2.
7. Document the results and any decisions made during the analysis.

---

---

# Español 🇪🇸

## 1. Descripción general del proyecto

Este repositorio forma parte de mi proceso de aprendizaje para crear un flujo de trabajo reproducible de análisis de bulk RNA-seq utilizando datos públicos de secuenciación de GEO, correspondientes al acceso **GSE201955**, que incluye datos de RNA-seq de células epiteliales bronquiales humanas.

Por ahora, solamente he seleccionado 3 muestras para probar el flujo de trabajo y entender qué ocurre en cada etapa.

El análisis cubre las primeras etapas del procesamiento de RNA-seq, desde las lecturas originales hasta el conteo de lecturas por gen. Mi objetivo no es solamente ejecutar las herramientas, sino también revisar los resultados, entender los posibles problemas y documentar las decisiones que voy tomando durante el proceso.

Hasta ahora he completado:

1. La descarga de los datos de secuenciación.
2. La conversión de los archivos SRA a FASTQ.
3. La evaluación de la calidad de las lecturas originales con FastQC.
4. El procesamiento de las lecturas con fastp.
5. La comparación entre los resultados del recorte original y el optimizado.
6. El alineamiento de las lecturas al genoma de referencia humano.
7. La revisión de los archivos BAM y las estadísticas de alineamiento.
8. El conteo de las lecturas asignadas a genes.
9. La integración de los resultados finales de control de calidad con MultiQC.

El siguiente paso será preparar la matriz de conteos y continuar el análisis en R.

**Todavía no he realizado el análisis de expresión diferencial.**

## 2. Datos y muestras

Estoy trabajando con las siguientes tres muestras:

| Muestra | Tipo de secuenciación | Longitud original de las lecturas |
|---|---|---|
| SRR19020908 | Single-end | 51 pb |
| SRR19020909 | Single-end | 51 pb |
| SRR19020910 | Single-end | 51 pb |

Los identificadores de las muestras están en [`samples.txt`](samples.txt).

Para el alineamiento utilicé el genoma de referencia humano **GRCh38**.

Para el conteo de genes utilicé la anotación Ensembl GRCh38, versión 115.

Este es un conjunto pequeño de muestras para probar el análisis. Antes de realizar comparaciones biológicas, todavía necesito revisar los metadatos de las muestras y el diseño experimental.

## 3. Estructura del repositorio

Los principales archivos y carpetas de este proyecto son:

```text
prueba_RNA-seq_analisis/
├── README.md
├── .gitignore
├── samples.txt
├── scripts/
│   ├── 01_download.sh
│   ├── 02_fastq.sh
│   ├── 03_fastp.sh
│   ├── 04_align.sh
│   └── 05_featurecounts.sh
├── docs/
│   └── fastqc_trimming_analysis.md
├── 00_raw/
├── 01_clean/
├── 02_bam/
├── 03_counts/
├── 04_multiqc/
│   ├── multiqc_config.yaml
│   └── multiqc_final_tail1.html
└── ref/
```

Algunas carpetas contienen archivos grandes que solamente guardo de manera local.

No subo los archivos FASTQ originales, los archivos BAM, los archivos de referencia ni las tablas de conteos a GitHub. En su lugar, mantengo los scripts, la documentación y algunos informes seleccionados que permiten entender el análisis.

## 4. Herramientas y flujo de trabajo

Realicé el análisis en un entorno Conda llamado `rnaseq_qc`.

Estas son las principales herramientas que utilicé:

| Herramienta | Para qué la utilicé |
|---|---|
| SRA Toolkit | Descargar y convertir los datos de secuenciación |
| FastQC | Revisar la calidad de las lecturas |
| fastp | Filtrar y recortar las lecturas |
| Subread | Alinear las lecturas a GRCh38 |
| SAMtools | Revisar los archivos BAM y las estadísticas de alineamiento |
| featureCounts | Contar las lecturas asignadas a genes |
| MultiQC | Combinar los resultados en informes de control de calidad |

### Flujo de trabajo

```text
Descarga desde SRA
       |
       v
Conversión a FASTQ
       |
       v
FastQC de las lecturas originales
       |
       v
Recorte con fastp
       |
       v
FastQC después del recorte
       |
       v
Alineamiento a GRCh38
       |
       v
Validación BAM + flagstat
       |
       v
featureCounts
       |
       v
Informe final de MultiQC
       |
       v
Siguiente etapa: análisis en R
```

Los scripts están disponibles en la carpeta [`scripts/`](scripts/).

Deben ejecutarse desde el directorio principal del proyecto.

## 5. Control de calidad y recorte

### 5.1 Qué encontré en las lecturas originales

Primero revisé los archivos FASTQ originales con FastQC.

Los principales problemas fueron:

- Menor calidad de secuenciación en la última base (posición 51).
- Mayor contenido de `N` en la última base.
- Composición desigual de nucleótidos al inicio de las lecturas, principalmente durante las primeras 13 posiciones.

El problema de la última base aparecía en las tres muestras, por lo que decidí probar si eliminar esa base mejoraba los resultados.

No eliminé las primeras 13 bases porque el patrón de composición nucleotídica podría estar relacionado con la preparación de la biblioteca de RNA-seq. No quería recortar esas posiciones sin tener una razón clara.

### 5.2 Recorte original y optimizado

Primero procesé las lecturas utilizando fastp con su configuración estándar y guardé los archivos de salida con el sufijo:

```text
_clean.fastq.gz
```

Después de revisar los resultados de FastQC, probé una segunda estrategia añadiendo un parámetro a fastp:

```bash
--trim_tail1 1
```

Este parámetro elimina una base del extremo 3′ de cada lectura.

Los archivos optimizados tienen el sufijo:

```text
_tail1.fastq.gz
```

Conservé las dos versiones para poder comparar los resultados sin reemplazar los archivos originales.

### 5.3 ¿Cuántas lecturas se conservaron?

| Muestra | Lecturas originales | Lecturas después de fastp optimizado | Retención |
|---|---:|---:|---:|
| SRR19020908 | 3,075,016 | 3,055,200 | 99.36% |
| SRR19020909 | 3,562,425 | 3,543,079 | 99.46% |
| SRR19020910 | 3,629,875 | 3,615,115 | 99.59% |

Los porcentajes Q30 después del filtrado fueron aproximadamente del 99% en las tres muestras.

La estrategia optimizada también conservó un poco más de lecturas que el procesamiento original con fastp:

| Muestra | Lecturas adicionales conservadas |
|---|---:|
| SRR19020908 | +929 |
| SRR19020909 | +1,049 |
| SRR19020910 | +1,100 |

Las diferencias son pequeñas, pero la estrategia optimizada dio mejores resultados para el problema que estaba intentando corregir.

### 5.4 FastQC después del recorte optimizado

Después de eliminar la última base, ejecuté FastQC nuevamente.

| Módulo FastQC | Resultado |
|---|---|
| Per base sequence quality | PASS |
| Per sequence quality scores | PASS |
| Per base N content | PASS |
| Per base sequence content | FAIL |
| Per sequence GC content | PASS |
| Adapter Content | PASS |
| Overrepresented sequences | PASS |
| Sequence Length Distribution | WARNING |

Los problemas de calidad y contenido de `N` en la última base se corrigieron.

El módulo **Per base sequence content** todavía presentó un resultado FAIL.

Decidí no realizar un recorte adicional solamente para conseguir que todos los módulos de FastQC mostraran PASS. Considero más importante entender por qué aparece una advertencia que eliminar datos sin saber si realmente es necesario.

El análisis completo y la comparación entre los diferentes resultados de fastp están documentados aquí:

[FastQC and trimming analysis](docs/fastqc_trimming_analysis.md)

## 6. Alineamiento al genoma de referencia

Después de seleccionar las lecturas optimizadas, las alineé al genoma de referencia humano GRCh38 utilizando Subread.

El script es:

[`scripts/04_align.sh`](scripts/04_align.sh)

Guardé los nuevos archivos BAM con el sufijo `_tail1.bam` para no sobrescribir los archivos BAM originales y poder conservarlos para comparaciones posteriores.

### Resultados del alineamiento

| Muestra | Lecturas de entrada | Lecturas alineadas | Porcentaje de alineamiento |
|---|---:|---:|---:|
| SRR19020908 | 3,055,200 | 2,714,036 | 88.83% |
| SRR19020909 | 3,543,079 | 3,156,878 | 89.10% |
| SRR19020910 | 3,615,115 | 3,228,911 | 89.32% |

Los porcentajes de alineamiento fueron similares entre las tres muestras.

También revisé los archivos BAM utilizando:

```bash
samtools quickcheck -v 02_bam/*_tail1.bam
```

Los tres archivos pasaron la validación.

Después generé los informes `flagstat` de SAMtools para incluir las estadísticas de alineamiento en MultiQC.

## 7. Conteo de genes con featureCounts

Una vez que los archivos BAM estuvieron listos, utilicé featureCounts para contar las lecturas asignadas a genes.

El script es:

[`scripts/05_featurecounts.sh`](scripts/05_featurecounts.sh)

Utilicé la siguiente anotación:

```text
ref/Homo_sapiens.GRCh38.115.gtf.gz
```

Los parámetros principales fueron:

```bash
-T 2
-t exon
-g gene_id
-s 0
```

Comparé las diferentes opciones de especificidad de cadena en featureCounts para ver cuál permitía asignar más lecturas a los genes.

La opción no específica de cadena (`-s 0`) produjo el mayor número de lecturas asignadas, por lo que decidí utilizarla para los conteos finales. Este resultado respalda la elección de un conteo no específico de cadena, aunque el número de lecturas asignadas por sí solo no confirma de manera definitiva el método de preparación de la biblioteca.

En otros conjuntos de datos de RNA-seq, esta información puede estar disponible en la documentación del estudio. Es útil revisar los metadatos de GEO o SRA, la publicación original y, especialmente, los métodos de preparación de la biblioteca o el protocolo de secuenciación.

Si esta información no está disponible o no es clara, se puede investigar la especificidad de cadena comparando los resultados de featureCounts con `-s 0`, `-s 1` y `-s 2`, o utilizando una herramienta como RSeQC (`infer_experiment.py`).

Es importante comprobar este aspecto antes del conteo de genes porque utilizar una configuración incorrecta puede afectar la cantidad de lecturas asignadas a los genes.

Los resultados optimizados se guardaron en:

```text
03_counts/gene_counts_tail1.txt
```

### Por qué repetí el análisis

Inicialmente procesé las lecturas utilizando fastp y continué con el alineamiento y el conteo de genes. En ese momento todavía no había realizado un análisis independiente con FastQC.

Después de revisar mi flujo de trabajo con mi profesor, me señaló que había omitido FastQC. Me explicó que, aunque fastp genera estadísticas e informes de calidad, también era importante revisar los datos de secuenciación de manera independiente antes de decidir si era necesario realizar algún recorte adicional.

Por eso regresé a los archivos FASTQ originales y ejecuté FastQC. Esto me permitió identificar la menor calidad y el aumento del contenido de `N` en la última base, además del sesgo de composición nucleotídica al inicio de las lecturas.

Con base en estos resultados, decidí probar una estrategia optimizada de fastp utilizando `--trim_tail1 1`. Después repetí FastQC, el alineamiento y el conteo de genes para poder comparar ambos procedimientos.

Utilizar FastQC junto con fastp fue útil porque las dos herramientas proporcionan información complementaria. FastQC genera figuras claras sobre la calidad a lo largo de las lecturas, la composición de nucleótidos, el contenido de GC y otros posibles problemas. Estas figuras me ayudaron a identificar qué aspectos necesitaban atención y a evaluar con más detalle los resultados de fastp.

Esta comparación también me permitió entender que las decisiones de recorte deben basarse en la calidad de los datos y no solamente en que una herramienta haya terminado correctamente.

### Resultados de asignación a genes

| Muestra | Lecturas asignadas | Porcentaje de asignación |
|---|---:|---:|
| SRR19020908 | 2,093,666 | 68.5% |
| SRR19020909 | 2,563,479 | 72.4% |
| SRR19020910 | 2,572,873 | 71.2% |

También comparé los conteos de genes del procesamiento inicial con los del procesamiento optimizado.

| Muestra | Lecturas asignadas originales | Lecturas asignadas optimizadas | Diferencia |
|---|---:|---:|---:|
| SRR19020908 | 2,093,495 | 2,093,666 | +171 |
| SRR19020909 | 2,563,077 | 2,563,479 | +402 |
| SRR19020910 | 2,572,872 | 2,572,873 | +1 |

Las diferencias fueron muy pequeñas.

Esto sugiere que eliminar la última base mejoró la calidad de las lecturas sin cambiar de manera importante los resultados del conteo de genes.

## 8. Informes de MultiQC

Generé informes de MultiQC por separado porque quería mantener claras las diferentes etapas del procesamiento.

Los informes son:

1. `multiqc_before.html` — FastQC antes del recorte.
2. `multiqc_after.html` — FastQC después del recorte original.
3. `multiqc_tail1.html` — FastQC después del recorte optimizado.
4. `multiqc_fastp_comparison.html` — comparación de las dos estrategias de fastp.
5. `multiqc_final_tail1.html` — informe final del flujo optimizado.

Los informes intermedios se guardan localmente.

### Informe final de MultiQC

Para el informe final reuní:

- Los resultados de FastQC de las lecturas optimizadas.
- Los resultados de fastp del procesamiento optimizado.
- Las estadísticas de SAMtools `flagstat`.
- Las estadísticas de asignación de featureCounts.

Al principio, MultiQC mostraba varias filas para una misma muestra porque las diferentes herramientas utilizaban nombres ligeramente distintos.

Resolví esto mediante un archivo de configuración que unifica los nombres de las muestras.

Ahora la tabla **General Statistics** muestra tres filas, una por muestra, con los resultados de las diferentes herramientas integrados.

**Informe final:**

[Descargar el informe final de MultiQC](04_multiqc/multiqc_final_tail1.html)

**Archivo de configuración:**

[Ver la configuración de MultiQC](04_multiqc/multiqc_config.yaml)

**Nota:** GitHub puede no mostrar una vista previa del archivo HTML debido a su tamaño. Para ver el informe interactivo, hay que descargar el archivo HTML original y abrirlo en un navegador.

## 9. Qué aprendí de los resultados

El flujo optimizado funcionó para las tres muestras.

Los resultados principales fueron:

- Se conservó más del 99% de las lecturas después del filtrado.
- El porcentaje Q30 fue aproximadamente del 99% en las tres muestras.
- El recorte de la última base corrigió los problemas de calidad y contenido de `N`.
- Los porcentajes de alineamiento estuvieron entre 88.83% y 89.32%.
- Los porcentajes de asignación a genes estuvieron entre 68.5% y 72.4%.
- Todos los archivos BAM optimizados pasaron la validación.
- Los conteos de genes originales y optimizados fueron muy similares.

## 10. Próximos pasos

La siguiente parte del proyecto estará enfocada en trabajar con los conteos de genes en R.

Tengo previsto:

1. Importar los resultados de featureCounts.
2. Preparar la matriz de conteos.
3. Revisar los metadatos de las muestras y los grupos experimentales.
4. Confirmar la especificidad de cadena de la biblioteca.
5. Explorar los datos mediante PCA y agrupamiento de muestras.
6. Revisar si el diseño experimental permite realizar un análisis de expresión diferencial con DESeq2.
7. Documentar los resultados y las decisiones tomadas durante el análisis.


---

**Repository / Repositorio:**  
[prueba_RNA-seq_analisis](https://github.com/viancyrn/prueba_RNA-seq_analisis)
