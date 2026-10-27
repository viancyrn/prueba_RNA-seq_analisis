# Bulk RNA-seq Analysis Pipeline | Análisis de RNA-seq

**A reproducible workflow for preprocessing human bulk RNA-seq data.**  
**Un flujo de trabajo reproducible para el preprocesamiento de datos de RNA-seq humano.**

**Dataset:** GSE201955  
**Organism / Organismo:** *Homo sapiens*  
**Sequencing type / Tipo de secuenciación:** Bulk RNA-seq, single-end, 51 bp  
**Status / Estado:** In progress / En desarrollo

---

# 🇬🇧 English

## 1. Project overview

This repository contains a reproducible **bulk RNA-seq analysis pipeline** for publicly available sequencing data from human bronchial epithelial cells.

The dataset **GSE201955** includes samples from healthy controls and individuals with asthma.

The objective is to develop, automate, and validate a workflow covering data retrieval, quality assessment, preprocessing, read alignment, and gene-level quantification.

The pipeline is being tested with three SRA accessions. Initial alignment and gene counting have been completed, and the quality-control strategy has been refined through three stages of FastQC assessment.

> **Note:** The current results are for technical pipeline validation and quality assessment. They are not sufficient for differential gene expression analysis.

## 2. Dataset and test samples

**Public data source:** NCBI Gene Expression Omnibus (GEO)

| Identifier | Accession |
|---|---|
| GEO Series | [GSE201955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201955) |
| BioProject | [PRJNA833856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833856) |
| SRA Study | SRP373003 |

### Samples used for pipeline testing

| SRA accession | Description |
|---|---|
| `SRR19020908` | Healthy control (provisional classification) |
| `SRR19020909` | Asthma (provisional classification) |
| `SRR19020910` | Experimental group pending verification |

The third accession was added to verify that the scripts can automatically process multiple samples.

Final experimental group assignments must be verified against GEO metadata before differential gene expression analysis.

## 3. Pipeline workflow

The following diagram represents the intended workflow and includes the optimized quality-control strategy.

```text
SRA accessions (samples.txt)
           |
           v
1. Download SRA files
        prefetch
           |
           v
2. Convert SRA to FASTQ
       fasterq-dump
           |
           v
3. Initial quality assessment
     FastQC (raw reads)
           |
           v
4. Quality filtering and trimming
    fastp (--trim_tail1 1)
           |
           v
5. Post-trimming quality assessment
   FastQC (processed reads)
           |
           v
6. Reference genome preparation
         GRCh38
           |
           v
7. Read alignment
         Subread
           |
           v
8. Library strandedness assessment
   featureCounts (-s 0, 1, 2)
           |
           v
9. Gene-level quantification
     featureCounts (-s 0)
           |
           v
10. Integrated quality report
          MultiQC
           |
           v
11. Gene expression analysis
         (planned)
```

### Pipeline progress

| Step | Tool | Status |
|---|---|---|
| 1. SRA download | SRA Toolkit (`prefetch`) | ✅ Completed |
| 2. SRA-to-FASTQ conversion | SRA Toolkit (`fasterq-dump`) | ✅ Completed |
| 3. Initial quality assessment | FastQC | ✅ Completed |
| 4. Quality filtering and trimming optimization | `fastp` | ✅ Completed |
| 5. Post-trimming quality assessment | FastQC | ✅ Completed |
| 6. Reference genome preparation | GRCh38 / Subread | ✅ Completed |
| 7. Initial read alignment | Subread | ✅ Completed |
| 8. Initial strandedness assessment | featureCounts (`-s 0`, `-s 1`, `-s 2`) | ✅ Completed |
| 9. Initial gene-level quantification | featureCounts (`-s 0`) | ✅ Completed |
| 10. Integrated quality report | MultiQC | ⏳ Pending |
| 11. Gene expression analysis | To be determined | ⏳ Planned |

> **Note:** Initial alignment, strandedness assessment, and gene counting were performed before the final trimming strategy was selected. These steps have not yet been repeated using the optimized `_tail1.fastq.gz` files.

### Quality-control optimization

FastQC was performed at three stages:

1. Before trimming, using the original FASTQ files.
2. After the initial fastp run, using default trimming settings.
3. After optimized fastp trimming, using `--trim_tail1 1`.

The comparison of these stages is documented in [FastQC analysis and trimming optimization](docs/fastqc_trimming_analysis.md).

## 4. Repository structure

The repository contains scripts, sample metadata, documentation, and selected FastQC reports.

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
│   ├── fastqc_before/
│   ├── fastqc_after/
│   └── fastqc_tail1/
└── ref/
```

### Directory descriptions

| Directory | Purpose |
|---|---|
| `scripts/` | Bash scripts for downloading, preprocessing, alignment, and gene counting |
| `docs/` | Detailed analyses and documentation of quality-control decisions |
| `00_raw/` | Original SRA files and raw FASTQ reads |
| `01_clean/` | Processed FASTQ files and fastp HTML/JSON reports |
| `02_bam/` | Aligned reads in BAM format |
| `03_counts/` | Gene-level count matrix and featureCounts summaries |
| `04_multiqc/fastqc_before/` | FastQC reports for original reads |
| `04_multiqc/fastqc_after/` | FastQC reports after the initial fastp run |
| `04_multiqc/fastqc_tail1/` | FastQC reports after optimized trimming |
| `ref/` | Reference genome, annotation files, and alignment index |

> **Note:** Large sequencing files, reference files, and generated analysis outputs are stored locally. Selected FastQC HTML reports are explicitly included in Git for documentation, while large data files remain excluded.

## 5. Software and environment

The pipeline is developed on **macOS (Apple Silicon)** using a Conda environment named `rnaseq_qc`. FastQC was installed separately using Homebrew.

| Software | Version | Purpose | Status |
|---|---|---|---|
| SRA Toolkit | 3.4.1 | Download and FASTQ conversion | Used |
| fastp | 1.4.0 | Quality filtering and trimming | Used |
| FastQC | 0.13.0 | Quality assessment before and after trimming | Three assessment stages completed |
| Subread | 2.1.1 | Genome indexing and read alignment | Used |
| featureCounts | Included in Subread 2.1.1 | Gene-level quantification and strandedness assessment | Used |
| samtools | 1.24 | BAM validation and inspection | Used |
| MultiQC | 1.35 | Integrated quality reporting | Installed; report pending |
| Bash | System installation | Workflow automation | Used |

### Software usage and validation

- SRA Toolkit was used to download and convert three SRA accessions into FASTQ files.
- fastp was used for quality filtering and trimming. The final strategy includes `--trim_tail1 1`.
- FastQC was used to evaluate raw reads, reads after the initial fastp run, and reads after optimized trimming.
- Subread was used for initial alignment against the GRCh38 reference genome.
- samtools was used to validate the resulting BAM files with `samtools quickcheck`.
- featureCounts was used for initial gene-level quantification and comparison of strandedness settings (`-s 0`, `-s 1`, and `-s 2`).
- MultiQC is installed and will be used to integrate quality-control results.

### FastQC installation

FastQC was installed using Homebrew:

```bash
brew install fastqc
```

The installed version is `0.13.0`. The output directory is specified using `--outdir`.

For example, the original FASTQ files were analyzed using:

```bash
mkdir -p 04_multiqc/fastqc_before

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq
```

FastQC generates HTML and ZIP reports containing the quality-control results and figures.

## 6. How to run the preprocessing workflow

The following commands are executed from the project root directory.

### Step 1 — Activate the environment

```bash
conda activate rnaseq_qc
```

### Step 2 — Check the sample list

The `samples.txt` file contains one SRA accession per line:

```text
SRR19020908
SRR19020909
SRR19020910
```

### Step 3 — Download SRA files

```bash
bash scripts/01_download.sh
```

This script downloads missing SRA files and stores them in `00_raw/`.

### Step 4 — Convert SRA files to FASTQ

```bash
bash scripts/02_fastq.sh
```

This script converts the downloaded SRA files to single-end FASTQ files.

### Step 5 — Assess the quality of original reads

```bash
mkdir -p 04_multiqc/fastqc_before

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq
```

This assessment was completed for all three samples.

### Step 6 — Perform quality filtering and optimized trimming

```bash
bash scripts/03_fastp.sh
```

The current script applies quality filtering and removes the final nucleotide from each read using:

```bash
--trim_tail1 1
```

The processed FASTQ files use the suffix `_tail1.fastq.gz`.

The script also generates HTML and JSON quality-control reports.

The initial fastp run used default trimming settings and produced files with the suffix `_clean.fastq.gz`. Those results were retained for comparison during optimization.

### Step 7 — Assess quality after optimized trimming

```bash
mkdir -p 04_multiqc/fastqc_tail1

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_tail1 \
  01_clean/*_tail1.fastq.gz
```

FastQC assessment was completed for all three samples.

The detailed comparison of all three quality-control stages is available in [FastQC analysis and trimming optimization](docs/fastqc_trimming_analysis.md).

### Step 8 — Align reads

```bash
bash scripts/04_align.sh
```

Initial alignment using Subread was completed for all three samples. BAM files were stored in `02_bam/` and validated using samtools.

**Important:** The initial alignment was performed before the optimized trimming strategy was selected. The alignment script must be updated to use `_tail1.fastq.gz` before repeating this step.

### Step 9 — Quantify reads at gene level

```bash
bash scripts/05_featurecounts.sh
```

Initial quantification used exon annotations, the `gene_id` attribute, and the unstranded setting (`-s 0`).

The count matrix was saved as:

`03_counts/gene_counts.txt`

**Important:** Gene counting must be repeated after alignment with the final trimmed reads.

### Step 10 — Integrate quality-control reports

MultiQC reporting is pending.

The objective is to combine quality-control results into a single report for easier comparison between samples and processing stages.

## 7. Preliminary results

### 7.1. Quality filtering and optimized trimming with fastp

The three samples were processed successfully using the optimized fastp strategy.

| Sample | Reads before filtering | Reads after filtering | Reads retained |
|---|---:|---:|---:|
| `SRR19020908` | 3,075,016 | 3,055,200 | 99.36% |
| `SRR19020909` | 3,562,425 | 3,543,079 | 99.46% |
| `SRR19020910` | 3,629,875 | 3,615,115 | 99.59% |

**Main observations:**

- More than **99% of reads** were retained in each sample.
- Approximately **99% of bases achieved Q30 quality** after filtering.
- The optimized strategy removes one nucleotide from the 3' end of each read.
- Processed reads are stored as compressed FASTQ files.

### 7.2. Initial FastQC assessment

FastQC reports were generated and reviewed for all three original FASTQ files.

| FastQC module | SRR19020908 | SRR19020909 | SRR19020910 |
|---|---|---|---|
| Basic Statistics | PASS | PASS | PASS |
| Per base sequence quality | FAIL | FAIL | FAIL |
| Per sequence quality scores | PASS | PASS | PASS |
| Per base sequence content | FAIL | FAIL | FAIL |
| Per sequence GC content | PASS | PASS | PASS |
| Per base N content | FAIL | FAIL | FAIL |
| Sequence Length Distribution | PASS | PASS | PASS |
| Sequence Duplication Levels | PASS | PASS | PASS |
| Overrepresented sequences | WARNING | WARNING | WARNING |
| Adapter Content | PASS | PASS | PASS |

The three samples showed similar quality patterns.

The initial analysis identified reduced sequencing quality and increased N content at nucleotide position 51.

Uneven nucleotide composition was also observed during approximately the first 13 positions. However, this pattern does not automatically indicate poor sequencing quality in RNA-seq data.

### 7.3. FastQC assessment and trimming optimization

FastQC reports were compared at three stages:

1. Before trimming.
2. After the initial fastp run.
3. After optimized trimming with `--trim_tail1 1`.

The final nucleotide was removed because the original FastQC results showed reduced quality and increased N content at position 51.

The first 13 nucleotides were retained because uneven base composition can be associated with RNA-seq library preparation and does not necessarily justify removing a substantial part of short reads.

After optimized trimming, the FastQC modules **Per base sequence quality** and **Per base N content** passed for all three samples.

The **Per base sequence content** module continued to fail, but this result alone was not considered sufficient evidence to remove the first 13 bases.

**Detailed analysis and supporting FastQC reports:**

- [FastQC analysis and trimming decisions](docs/fastqc_trimming_analysis.md)
- [FastQC reports before trimming](04_multiqc/fastqc_before/)
- [FastQC reports after initial fastp trimming](04_multiqc/fastqc_after/)
- [FastQC reports after optimized trimming](04_multiqc/fastqc_tail1/)

Each report directory contains the HTML reports for all three samples, including the quality-control figures.

### 7.4. Initial alignment and gene-level quantification

Initial alignment was completed for all three samples, and the resulting BAM files passed validation with `samtools quickcheck`.

The initial featureCounts analysis produced the following results:

| Sample | Reads assigned to genes | Assignment rate |
|---|---:|---:|
| `SRR19020908` | 2,093,495 | 68.5% |
| `SRR19020909` | 2,563,077 | 72.4% |
| `SRR19020910` | 2,572,872 | 71.2% |

The strandedness settings `-s 0`, `-s 1`, and `-s 2` were compared. The unstranded setting (`-s 0`) was selected for initial gene quantification.

The initial count matrix is located at:

`03_counts/gene_counts.txt`

**These alignment and counting results were obtained before the final trimming strategy was selected and must be reassessed using the optimized reads.**

The current results are intended for technical pipeline validation and do not establish biological differences between experimental groups.

## 8. Reproducibility and data management

- The `samples.txt` file controls sample selection.
- Bash scripts automate the main processing stages.
- Samples are processed sequentially where necessary to limit memory usage.
- Original sequencing files are preserved during validation.
- Large sequencing files and reference files are excluded from Git tracking.
- Selected FastQC HTML reports are included to document quality-control decisions.
- Scripts and documentation are version-controlled using Git.
- Software versions and strandedness parameters are documented to support reproducibility.

### Next steps

1. Update the alignment script to use the final `_tail1.fastq.gz` files.
2. Repeat alignment and gene-level quantification using the optimized reads.
3. Generate an integrated quality-control report with MultiQC.
4. Verify GEO metadata and experimental group assignments.
5. Plan the gene expression analysis.

---

# 🇨🇴 Español

## 1. Descripción del proyecto

Este repositorio contiene un **pipeline reproducible de análisis de RNA-seq de tipo bulk** para procesar datos públicos de secuenciación de células epiteliales bronquiales humanas.

El conjunto de datos **GSE201955** incluye muestras de controles sanos y personas con asma.

El objetivo es desarrollar, automatizar y validar un flujo de trabajo que abarque la descarga de datos, la evaluación de calidad, el preprocesamiento, el alineamiento y la cuantificación de lecturas por gen.

El pipeline se está probando con tres accesiones SRA. El alineamiento y la cuantificación iniciales ya se completaron, y la estrategia de control de calidad se optimizó mediante tres etapas de evaluación con FastQC.

> **Nota:** Los resultados actuales corresponden a la validación técnica del pipeline y la evaluación de calidad. Todavía no constituyen un análisis de expresión diferencial.

## 2. Datos y muestras

**Fuente de datos:** NCBI Gene Expression Omnibus (GEO)

| Identificador | Accesión |
|---|---|
| GEO Series | [GSE201955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201955) |
| BioProject | [PRJNA833856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833856) |
| SRA Study | SRP373003 |

### Muestras utilizadas para validar el pipeline

| Accesión SRA | Descripción |
|---|---|
| `SRR19020908` | Control sano (clasificación provisional) |
| `SRR19020909` | Asma (clasificación provisional) |
| `SRR19020910` | Grupo experimental pendiente de verificación |

La tercera muestra se agregó para comprobar que los scripts pueden procesar automáticamente múltiples accesiones.

Las asignaciones definitivas de los grupos experimentales deberán verificarse con los metadatos de GEO antes del análisis de expresión diferencial.

## 3. Etapas del pipeline

El siguiente diagrama representa el flujo de trabajo previsto e incluye la estrategia optimizada de control de calidad.

```text
Accesiones SRA (samples.txt)
             |
             v
1. Descarga de archivos SRA
           prefetch
             |
             v
2. Conversión de SRA a FASTQ
         fasterq-dump
             |
             v
3. Evaluación inicial de calidad
    FastQC (lecturas originales)
             |
             v
4. Filtrado de calidad y trimming
     fastp (--trim_tail1 1)
             |
             v
5. Evaluación después del trimming
    FastQC (lecturas procesadas)
             |
             v
6. Preparación del genoma de referencia
             GRCh38
             |
             v
7. Alineamiento de lecturas
             Subread
             |
             v
8. Evaluación de strandedness
   featureCounts (-s 0, 1, 2)
             |
             v
9. Cuantificación de lecturas por gen
       featureCounts (-s 0)
             |
             v
10. Informe integrado de calidad
             MultiQC
             |
             v
11. Análisis de expresión génica
          (planificado)
```

### Progreso del pipeline

| Etapa | Herramienta | Estado |
|---|---|---|
| 1. Descarga de archivos SRA | SRA Toolkit (`prefetch`) | ✅ Completada |
| 2. Conversión de SRA a FASTQ | SRA Toolkit (`fasterq-dump`) | ✅ Completada |
| 3. Evaluación inicial de calidad | FastQC | ✅ Completada |
| 4. Filtrado y optimización del trimming | `fastp` | ✅ Completada |
| 5. Evaluación después del trimming | FastQC | ✅ Completada |
| 6. Preparación del genoma de referencia | GRCh38 / Subread | ✅ Completada |
| 7. Alineamiento inicial de lecturas | Subread | ✅ Completada |
| 8. Evaluación inicial de strandedness | featureCounts (`-s 0`, `-s 1`, `-s 2`) | ✅ Completada |
| 9. Cuantificación inicial por gen | featureCounts (`-s 0`) | ✅ Completada |
| 10. Informe integrado de calidad | MultiQC | ⏳ Pendiente |
| 11. Análisis de expresión génica | Por determinar | ⏳ Planificado |

> **Nota:** El alineamiento, la evaluación de strandedness y la cuantificación iniciales se realizaron antes de seleccionar la estrategia definitiva de trimming. Estas etapas todavía no se han repetido utilizando los archivos optimizados `_tail1.fastq.gz`.

### Optimización del control de calidad

Se realizó una evaluación con FastQC en tres etapas:

1. Antes del trimming, utilizando los archivos FASTQ originales.
2. Después de la ejecución inicial de fastp, utilizando la configuración predeterminada.
3. Después del trimming optimizado, utilizando `--trim_tail1 1`.

La comparación de estas etapas se encuentra en [FastQC analysis and trimming optimization](docs/fastqc_trimming_analysis.md).

## 4. Organización del proyecto

El repositorio contiene scripts, metadatos de muestras, documentación e informes FastQC seleccionados.

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
│   ├── fastqc_before/
│   ├── fastqc_after/
│   └── fastqc_tail1/
└── ref/
```

### Descripción de los directorios

| Directorio | Contenido |
|---|---|
| `scripts/` | Scripts Bash para descarga, preprocesamiento, alineamiento y cuantificación |
| `docs/` | Análisis detallados y documentación de las decisiones de calidad |
| `00_raw/` | Archivos SRA y FASTQ originales |
| `01_clean/` | Archivos FASTQ procesados e informes HTML/JSON de fastp |
| `02_bam/` | Alineamientos en formato BAM |
| `03_counts/` | Matriz de conteos por gen y resúmenes de featureCounts |
| `04_multiqc/fastqc_before/` | Informes FastQC de las lecturas originales |
| `04_multiqc/fastqc_after/` | Informes FastQC después del procesamiento inicial con fastp |
| `04_multiqc/fastqc_tail1/` | Informes FastQC después del trimming optimizado |
| `ref/` | Genoma de referencia, anotaciones e índice de alineamiento |

> **Nota:** Los archivos de secuenciación, las referencias y los resultados de gran tamaño se almacenan localmente. Se incluyen explícitamente en Git algunos informes HTML de FastQC para documentar las decisiones de control de calidad.

## 5. Software y entorno de trabajo

El pipeline se desarrolla en **macOS (Apple Silicon)** utilizando un entorno Conda llamado `rnaseq_qc`. FastQC se instaló por separado mediante Homebrew.

| Software | Versión | Función | Estado |
|---|---|---|---|
| SRA Toolkit | 3.4.1 | Descarga y conversión a FASTQ | Utilizado |
| fastp | 1.4.0 | Filtrado de calidad y trimming | Utilizado |
| FastQC | 0.13.0 | Evaluación antes y después del trimming | Tres etapas completadas |
| Subread | 2.1.1 | Preparación del índice y alineamiento | Utilizado |
| featureCounts | Incluido en Subread 2.1.1 | Cuantificación por gen y evaluación de strandedness | Utilizado |
| samtools | 1.24 | Validación e inspección de archivos BAM | Utilizado |
| MultiQC | 1.35 | Integración de informes de calidad | Instalado; informe pendiente |
| Bash | Instalación del sistema | Automatización del flujo de trabajo | Utilizado |

### Uso y validación del software

- SRA Toolkit se utilizó para descargar y convertir tres accesiones SRA en archivos FASTQ.
- fastp se utilizó para el filtrado de calidad y el trimming. La estrategia definitiva incluye `--trim_tail1 1`.
- FastQC se utilizó para evaluar las lecturas originales, las lecturas después del procesamiento inicial con fastp y las lecturas después del trimming optimizado.
- Subread se utilizó para el alineamiento inicial contra el genoma de referencia GRCh38.
- samtools se utilizó para validar los archivos BAM mediante `samtools quickcheck`.
- featureCounts se utilizó para la cuantificación inicial por gen y para comparar las configuraciones de strandedness (`-s 0`, `-s 1` y `-s 2`).
- MultiQC está instalado y se utilizará para integrar los resultados de calidad.

### Instalación de FastQC

FastQC se instaló utilizando Homebrew:

```bash
brew install fastqc
```

La versión instalada es `0.13.0`. El directorio de salida se especifica mediante `--outdir`.

Por ejemplo, las lecturas originales se analizaron utilizando:

```bash
mkdir -p 04_multiqc/fastqc_before

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq
```

FastQC genera informes HTML y ZIP que contienen los resultados y las figuras de control de calidad.

## 6. Ejecución del pipeline

Los siguientes comandos se ejecutan desde el directorio principal del proyecto.

### Paso 1 — Activar el entorno

```bash
conda activate rnaseq_qc
```

### Paso 2 — Verificar la lista de muestras

El archivo `samples.txt` contiene una accesión SRA por línea:

```text
SRR19020908
SRR19020909
SRR19020910
```

### Paso 3 — Descargar los archivos SRA

```bash
bash scripts/01_download.sh
```

El script descarga las accesiones que faltan y almacena los archivos SRA en `00_raw/`.

### Paso 4 — Convertir SRA a FASTQ

```bash
bash scripts/02_fastq.sh
```

Este script convierte los archivos SRA en archivos FASTQ de lecturas single-end.

### Paso 5 — Evaluar la calidad de las lecturas originales

```bash
mkdir -p 04_multiqc/fastqc_before

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq
```

Esta evaluación se completó para las tres muestras.

### Paso 6 — Realizar el filtrado y trimming optimizado

```bash
bash scripts/03_fastp.sh
```

El script actual realiza el filtrado de calidad y elimina el último nucleótido de cada lectura mediante:

```bash
--trim_tail1 1
```

Los archivos FASTQ procesados utilizan el sufijo `_tail1.fastq.gz`.

El script también genera informes de calidad en formatos HTML y JSON.

La ejecución inicial de fastp utilizó la configuración predeterminada y produjo archivos con el sufijo `_clean.fastq.gz`. Estos resultados se conservaron para comparar las estrategias durante la optimización.

### Paso 7 — Evaluar la calidad después del trimming optimizado

```bash
mkdir -p 04_multiqc/fastqc_tail1

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_tail1 \
  01_clean/*_tail1.fastq.gz
```

La evaluación con FastQC se completó para las tres muestras.

La comparación detallada de las tres etapas de control de calidad se encuentra en [FastQC analysis and trimming optimization](docs/fastqc_trimming_analysis.md).

### Paso 8 — Alinear las lecturas

```bash
bash scripts/04_align.sh
```

El alineamiento inicial con Subread se completó para las tres muestras. Los archivos BAM se almacenaron en `02_bam/` y fueron validados con samtools.

**Importante:** El alineamiento inicial se realizó antes de seleccionar la estrategia optimizada de trimming. El script de alineamiento debe actualizarse para utilizar los archivos `_tail1.fastq.gz` antes de repetir esta etapa.

### Paso 9 — Cuantificar las lecturas por gen

```bash
bash scripts/05_featurecounts.sh
```

La cuantificación inicial utilizó la anotación de exones, el atributo `gene_id` y la configuración no específica de cadena (`-s 0`).

La matriz de conteos se guardó en:

`03_counts/gene_counts.txt`

**Importante:** La cuantificación deberá repetirse después del alineamiento con las lecturas procesadas mediante la estrategia definitiva.

### Paso 10 — Integrar los informes de calidad

La generación del informe integrado con MultiQC está pendiente.

El objetivo es combinar los resultados de control de calidad en un solo informe para facilitar la comparación entre muestras y etapas del procesamiento.

## 7. Resultados preliminares

### 7.1. Filtrado y trimming optimizado con fastp

Las tres muestras se procesaron correctamente con la estrategia optimizada de fastp.

| Muestra | Lecturas iniciales | Lecturas después del filtrado | Porcentaje conservado |
|---|---:|---:|---:|
| `SRR19020908` | 3,075,016 | 3,055,200 | 99.36% |
| `SRR19020909` | 3,562,425 | 3,543,079 | 99.46% |
| `SRR19020910` | 3,629,875 | 3,615,115 | 99.59% |

**Observaciones principales:**

- Se conservó más del **99 % de las lecturas** en las tres muestras.
- Aproximadamente el **99 % de las bases alcanzó una calidad Q30** después del filtrado.
- La estrategia optimizada elimina un nucleótido del extremo 3' de cada lectura.
- Los archivos FASTQ procesados se almacenan en formato comprimido.

### 7.2. Evaluación inicial con FastQC

Se generaron y revisaron los informes FastQC de las tres muestras antes del trimming.

| Módulo de FastQC | SRR19020908 | SRR19020909 | SRR19020910 |
|---|---|---|---|
| Basic Statistics | PASS | PASS | PASS |
| Per base sequence quality | FAIL | FAIL | FAIL |
| Per sequence quality scores | PASS | PASS | PASS |
| Per base sequence content | FAIL | FAIL | FAIL |
| Per sequence GC content | PASS | PASS | PASS |
| Per base N content | FAIL | FAIL | FAIL |
| Sequence Length Distribution | PASS | PASS | PASS |
| Sequence Duplication Levels | PASS | PASS | PASS |
| Overrepresented sequences | WARNING | WARNING | WARNING |
| Adapter Content | PASS | PASS | PASS |

Las tres muestras presentaron patrones de calidad similares.

El análisis inicial identificó una disminución de la calidad de secuenciación y un aumento del contenido de N en la posición 51.

También se observó una composición desigual de nucleótidos durante aproximadamente las primeras 13 posiciones. Sin embargo, este patrón no indica automáticamente baja calidad de secuenciación en datos de RNA-seq.

### 7.3. Evaluación de FastQC y optimización del trimming

Se compararon los informes FastQC en tres etapas:

1. Antes del trimming.
2. Después de la ejecución inicial de fastp.
3. Después del trimming optimizado con `--trim_tail1 1`.

Se eliminó el último nucleótido porque los resultados iniciales de FastQC mostraron una disminución de calidad y un aumento del contenido de N en la posición 51.

Se conservaron los primeros 13 nucleótidos porque la composición desigual de bases puede estar relacionada con la preparación de bibliotecas de RNA-seq y no necesariamente justifica eliminar una parte considerable de lecturas cortas.

Después del trimming optimizado, los módulos **Per base sequence quality** y **Per base N content** obtuvieron resultados PASS en las tres muestras.

El módulo **Per base sequence content** continuó con un resultado FAIL, pero este resultado por sí solo no se consideró evidencia suficiente para eliminar las primeras 13 bases.

**Análisis detallado e informes FastQC con figuras:**

- [Análisis completo de FastQC y decisiones de trimming](docs/fastqc_trimming_analysis.md)
- [Informes FastQC antes del trimming](04_multiqc/fastqc_before/)
- [Informes FastQC después del trimming inicial](04_multiqc/fastqc_after/)
- [Informes FastQC después del trimming optimizado](04_multiqc/fastqc_tail1/)

Cada directorio contiene los informes HTML de las tres muestras, incluyendo las figuras de control de calidad.

### 7.4. Alineamiento inicial y cuantificación por gen

El alineamiento inicial se completó para las tres muestras y los archivos BAM superaron la validación con `samtools quickcheck`.

La cuantificación inicial con featureCounts produjo los siguientes resultados:

| Muestra | Lecturas asignadas a genes | Porcentaje asignado |
|---|---:|---:|
| `SRR19020908` | 2,093,495 | 68.5% |
| `SRR19020909` | 2,563,077 | 72.4% |
| `SRR19020910` | 2,572,872 | 71.2% |

Se compararon las configuraciones de strandedness `-s 0`, `-s 1` y `-s 2`. La configuración **no específica de cadena (`-s 0`)** fue seleccionada para la cuantificación inicial.

La matriz de conteos inicial se encuentra en:

`03_counts/gene_counts.txt`

**Estos resultados de alineamiento y cuantificación se obtuvieron antes de seleccionar la estrategia definitiva de trimming y deberán reevaluarse utilizando las lecturas optimizadas.**

Los resultados actuales corresponden a una evaluación técnica del procesamiento y no permiten establecer diferencias biológicas entre los grupos experimentales.

## 8. Reproducibilidad y gestión de datos

- El archivo `samples.txt` controla la selección de muestras.
- Los scripts Bash automatizan las principales etapas del pipeline.
- Las muestras se procesan de forma secuencial cuando es necesario para limitar el consumo de memoria.
- Los archivos originales se conservan durante la validación.
- Los archivos de secuenciación y las referencias de gran tamaño se excluyen del seguimiento de Git.
- Se incluyen informes HTML de FastQC seleccionados para documentar las decisiones de control de calidad.
- Los scripts y la documentación se mantienen bajo control de versiones con Git.
- Las versiones del software y los parámetros de strandedness se documentan para facilitar la reproducibilidad.

### Próximos pasos

1. Actualizar el script de alineamiento para utilizar los archivos definitivos `_tail1.fastq.gz`.
2. Repetir el alineamiento y la cuantificación por gen utilizando las lecturas optimizadas.
3. Generar un informe integrado con MultiQC.
4. Verificar los metadatos de GEO y la asignación de los grupos experimentales.
5. Planificar el análisis de expresión génica.

---

## References / Referencias

### Dataset and metadata / Datos y metadatos

- [NCBI GEO — GSE201955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201955)
- [NCBI BioProject — PRJNA833856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833856)
- [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)

### Bioinformatics software / Herramientas bioinformáticas

- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Subread and featureCounts](https://subread.sourceforge.net/)
- [samtools](https://www.htslib.org/)
- [MultiQC](https://multiqc.info/)

### Reference genome / Genoma de referencia

- [Genome Reference Consortium — Human GRCh38](https://www.ncbi.nlm.nih.gov/grc/human)

---

**Project status:** In progress — SRA download, FASTQ conversion, trimming optimization, three-stage FastQC assessment, and initial alignment and gene counting completed. Reanalysis with final trimmed reads and MultiQC reporting pending.

**Estado del proyecto:** En desarrollo — Descarga de SRA, conversión a FASTQ, optimización del trimming, evaluación FastQC en tres etapas y alineamiento y cuantificación iniciales completados. Pendientes el reprocesamiento con las lecturas definitivas y el informe MultiQC.
