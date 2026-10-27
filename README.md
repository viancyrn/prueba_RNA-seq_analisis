# Bulk RNA-seq Analysis Pipeline | Análisis de RNA-seq

**A reproducible workflow for preprocessing human bulk RNA-seq data.**  
**Un flujo de trabajo reproducible para el preprocesamiento de datos de RNA-seq humano.**

**Dataset:** GSE201955  
**Organism / Organismo:** *Homo sapiens*  
**Sequencing type / Tipo de secuenciación:** Bulk RNA-seq, single-end  
**Status / Estado:** In progress / En desarrollo

---

# 🇬🇧 English

## 1. Project overview

This repository contains a reproducible **bulk RNA-seq analysis pipeline** for publicly available sequencing data from human bronchial epithelial cells.

The dataset **GSE201955** includes samples from healthy controls and individuals with asthma.

The initial objective is to develop, automate, and validate the preprocessing workflow using three SRA accessions. The project will subsequently be extended to read alignment, gene quantification, and quality-report integration.

> **Note:** The current results are for pipeline testing and quality assessment. They are not sufficient for differential gene expression analysis.

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
| `SRR19020908` | Healthy control |
| `SRR19020909` | Asthma |
| `SRR19020910` | Asthma |

The third accession was added to verify that the scripts can automatically process multiple samples after testing the workflow first with one sample.

## 3. Pipeline workflow

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
3. Initial quality assessment before trimming
     FastQC (raw reads)
           |
           v
4. Quality and trimming
          fastp
           |
           v
5. Initial quality assessment after trimming
     FastQC (raw reads)
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
10. Summary quality report
         MultiQC
```

### Pipeline progress

| Step | Tool | Status |
|---|---|---|
| 1. SRA download | SRA Toolkit (`prefetch`) | ✅ Completed |
| 2. SRA-to-FASTQ conversion | SRA Toolkit (`fasterq-dump`) | ✅ Completed |
| 3. Initial quality assessment (raw reads) | FastQC | ✅ Completed |
| 4. Quality control and trimming | `fastp` | ✅ Completed |
| 5. Post-trimming quality assessment |FastQC |⏳ Pending |
| 6. Reference genome preparation | GRCh38 / Subread | ✅ Completed |
| 7. Read alignment | Subread | ✅ Completed |
| 8. Library strandedness assessment | featureCounts (-s 0, -s 1, -s 2) | ✅ Completed |
| 9. Gene-level read counting | featureCounts | ✅ Completed |
| 10. Integrated quality report | MultiQC | ⏳ Planned |

> **Note:** Initial FastQC analysis was performed retrospectively after alignment and read counting. Post-trimming FastQC and MultiQC reporting are the next quality-control steps.

## 4. Repository structure

The GitHub repository contains the scripts and metadata needed to document the workflow.

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
├── 00_raw/
├── 01_clean/
├── 02_bam/
├── 03_counts/
├── 04_multiqc/
│   ├── fastqc_before/
│   └── fastqc_after/
└── ref/
```

During execution, the pipeline also uses the following local directories:

| Directory | Purpose |
|---|---|
| `scripts/` | Bash scripts for downloading, preprocessing, alignment, and read counting |
| `00_raw/` | Original SRA files and raw FASTQ reads |
| `01_clean/` | Trimmed FASTQ files and fastp quality-control reports |
| `02_bam/` | Alignment reads in BAM format |
| `03_counts/` | Gene-level count matrix and featureCounts summaries |
| `04_multiqc/fastqc_before` | FastCQ reports for raw reads |
| `04_multiqc/fastqc_after` | FastCQ reports for trimmed reads |
| `ref/` | Reference genome, annotation files and alignment index |

> **Note:** Large sequencing files, reference files, and generated analysis outputs are stored locally and excluded from Git tracking where configured in .gitignore. GitHub is used primarily to version-control scripts, documentation, and sample metadata.

## 5. Software and environment

The pipeline is developed on **macOS (Apple Silicon)** using a Conda environment named `rnaseq_qc`.

| Software | Version | Purpose | Status | 
|---|---|---|---| 
| SRA Toolkit | 3.4.1 | Download and FASTQ conversion | Used |
| fastp | 1.4.0 | Quality control and trimming | Used |
| FastQC | 0.13.0 | Quality assessment before and after trimming | Raw reads analyzed; clean reads pending |
| Subread | 2.1.1 | Gene indexing and read alignment | Used |
| featureCounts | Included in Subread  2.1.1 | Gene-level quantification and strandedness assessment | Used |
| samtools | 1.24 | BAM file validation and inspection | Used |
| MultiQC | 1.35 | Quality-report integration | installed; reporting pending |
| Bash | — | Workflow automation | Used |

### Software usage and validation
- SRA Toolkit was used to download and convert three SRA accessions into FASTQ files.
- fastp was used to trim and quality-filter sequencing reads, retaining more than 99% of reads for each sample.
- FastQC was used to generate initial quality-control reports for all three raw FASTQ files.
- Subread was used for read alignment against the GRCh38 reference genome.
- samtools was used to validate the resulting BAM files using samtools quickcheck.
- featureCounts was used to quantify gene-level reads and compare strandedness settings (-s 0, -s 1, and -s 2). The unstranded setting (-s 0) was selected.
- MultiQC is installed and will be used to combine quality-control results after completing FastQC analysis of the trimmed reads.

### FastQC installation
FastQC was installed using Homebrew:
```bash
brew install fastqc

The installed version is 0.13.0. On this installation, the output directory is specified using --outdir.

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq

## 6. How to run the preprocessing workflow

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

This script reads the accession list, downloads missing SRA files, and stores them directly in `00_raw/`.

### Step 4 — Convert SRA files to FASTQ

```bash
bash scripts/02_fastq.sh
```

This script converts the downloaded SRA files to single-end FASTQ files.

### Step 5 — Perform quality control and trimming

```bash
bash scripts/03_fastp.sh
```

This script runs `fastp` sequentially for each sample and generates:

- A compressed, cleaned FASTQ file (`.fastq.gz`)
- An HTML quality-control report (`.html`)
- A JSON quality-control report (`.json`)

**The scripts skip samples with existing completed outputs, allowing the workflow to be resumed.**

## 7. Preliminary quality-control results

The three test samples were successfully processed with `fastp`.

| Sample | Reads before filtering | Reads after filtering | Reads retained |
|---|---:|---:|---:|
| `SRR19020908` | 3,075,016 | 3,054,271 | 99.33% |
| `SRR19020909` | 3,562,425 | 3,542,030 | 99.43% |
| `SRR19020910` | 3,629,875 | 3,614,015 | 99.56% |

### Main observations

- More than **99% of reads** were retained in each sample.
- The percentage of bases with **Q30 quality scores** after filtering was approximately 97%.
- `fastp` detected the **Illumina TruSeq Read 1 adapter**.
- Cleaned reads were saved as compressed FASTQ files to reduce storage requirements.

These observations support successful preprocessing, but do not establish biological differences between experimental groups.

## 8. Reproducibility and data management

- The workflow is controlled by the `samples.txt` accession list.
- Bash scripts automate the preprocessing steps.
- Samples are processed sequentially to limit resource usage.
- Raw sequencing data are preserved during testing.
- Large sequencing files are not uploaded to GitHub.
- Scripts and documentation are version-controlled with Git.

---

# 🇨🇴 Español

## 1. Descripción del proyecto

Este repositorio contiene un **pipeline reproducible de análisis de RNA-seq de tipo bulk** para procesar datos públicos de secuenciación de células epiteliales bronquiales humanas.

El conjunto de datos **GSE201955** incluye muestras de controles sanos y personas con asma.

El objetivo es desarrollar, automatizar y validar un flujo de trabajo que abarque la descarga de datos, el control de calidad, el preprocesamiento, el alineamiento y la cuantificación de lecturas por gen.

El pipeline se está probando con tres accesiones SRA. Posteriormente, se ampliará para realizar análisis de expresión génica, una vez verificados los metadatos y los grupos experimentales.

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

El siguiente diagrama representa el orden lógico recomendado del flujo de trabajo.

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
             fastp
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
| 4. Filtrado de calidad y trimming | `fastp` | ✅ Completada |
| 5. Evaluación después del trimming | FastQC | ⏳ Pendiente |
| 6. Preparación del genoma de referencia | GRCh38 / Subread | ✅ Completada |
| 7. Alineamiento de lecturas | Subread | ✅ Completada |
| 8. Evaluación de strandedness | featureCounts (`-s 0`, `-s 1`, `-s 2`) | ✅ Completada |
| 9. Cuantificación de lecturas por gen | featureCounts (`-s 0`) | ✅ Completada |
| 10. Informe integrado de calidad | MultiQC | ⏳ Pendiente |
| 11. Análisis de expresión génica | Por determinar | ⏳ Planificado |

> **Nota:** La evaluación inicial con FastQC se realizó retrospectivamente, después del alineamiento y la cuantificación. Los siguientes pasos de control de calidad son analizar las lecturas procesadas con FastQC y generar el informe integrado con MultiQC.

## 4. Organización del proyecto

El proyecto utiliza directorios separados para los datos originales, las lecturas procesadas, los alineamientos, los conteos por gen, los informes de calidad y los archivos de referencia.

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
├── 00_raw/
├── 01_clean/
├── 02_bam/
├── 03_counts/
├── 04_multiqc/
│   ├── fastqc_before/
│   └── fastqc_after/
└── ref/
```

### Descripción de los directorios

| Directorio | Contenido |
|---|---|
| `scripts/` | Scripts Bash para la descarga, el preprocesamiento, el alineamiento y la cuantificación |
| `00_raw/` | Archivos SRA y FASTQ originales |
| `01_clean/` | Archivos FASTQ procesados e informes de fastp |
| `02_bam/` | Alineamientos en formato BAM |
| `03_counts/` | Matriz de conteos por gen y resúmenes de featureCounts |
| `04_multiqc/fastqc_before/` | Informes FastQC de las lecturas originales |
| `04_multiqc/fastqc_after/` | Informes FastQC de las lecturas procesadas (pendientes) |
| `04_multiqc/` | Directorio destinado a los informes integrados de MultiQC |
| `ref/` | Genoma de referencia, anotaciones e índice de alineamiento |

> **Nota:** Los archivos de secuenciación, las referencias y los resultados generados se almacenan localmente y se excluyen del seguimiento de Git cuando así está configurado en `.gitignore`. GitHub se utiliza principalmente para mantener versiones de los scripts, la documentación y los metadatos de las muestras.

## 5. Software y entorno de trabajo

El pipeline se desarrolla en **macOS (Apple Silicon)** utilizando un entorno Conda llamado `rnaseq_qc`. FastQC se instaló por separado mediante Homebrew.

| Software | Versión | Función | Estado |
|---|---|---|---|
| SRA Toolkit | 3.4.1 | Descarga de datos y conversión a FASTQ | Utilizado |
| fastp | 1.4.0 | Filtrado, trimming e informes de calidad | Utilizado |
| FastQC | 0.13.0 | Evaluación de calidad antes y después del trimming | Lecturas originales analizadas; procesadas pendientes |
| Subread | 2.1.1 | Preparación del índice y alineamiento | Utilizado |
| featureCounts | Incluido en Subread 2.1.1 | Cuantificación por gen y evaluación de strandedness | Utilizado |
| samtools | 1.24 | Validación e inspección de archivos BAM | Utilizado |
| MultiQC | 1.35 | Integración de informes de calidad | Instalado; informe pendiente |
| Bash | Instalación del sistema | Automatización del flujo de trabajo | Utilizado |

### Uso y validación del software

- SRA Toolkit se utilizó para descargar y convertir tres accesiones SRA en archivos FASTQ.
- fastp se utilizó para filtrar y recortar las lecturas, conservando más del 99 % en cada muestra.
- FastQC se utilizó para generar informes de calidad de los tres archivos FASTQ originales.
- Subread se utilizó para alinear las lecturas contra el genoma de referencia GRCh38.
- samtools se utilizó para validar los archivos BAM mediante `samtools quickcheck`.
- featureCounts se utilizó para cuantificar lecturas por gen y comparar las configuraciones de strandedness (`-s 0`, `-s 1` y `-s 2`). Se seleccionó la configuración no específica de cadena (`-s 0`).
- MultiQC está instalado y se utilizará para integrar los resultados de calidad después de completar FastQC en las lecturas procesadas.

### Instalación de FastQC

FastQC se instaló utilizando Homebrew:

```bash
brew install fastqc
```

La versión instalada es `0.13.0`. En esta instalación, el directorio de salida se especifica mediante `--outdir`.

```bash
mkdir -p 04_multiqc/fastqc_before

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_before \
  00_raw/SRR19020908.fastq \
  00_raw/SRR19020909.fastq \
  00_raw/SRR19020910.fastq
```

Los informes HTML y ZIP generados se almacenan en `04_multiqc/fastqc_before/`.

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

El script descarga las accesiones que faltan y organiza los archivos SRA en `00_raw/`.

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

Este análisis ya se completó para las tres muestras.

### Paso 6 — Realizar el filtrado y trimming

```bash
bash scripts/03_fastp.sh
```

El script genera archivos FASTQ comprimidos y procesados, además de informes de calidad en formatos HTML y JSON.

### Paso 7 — Evaluar la calidad después del trimming

```bash
mkdir -p 04_multiqc/fastqc_after

fastqc -t 2 \
  --outdir 04_multiqc/fastqc_after \
  01_clean/*_clean.fastq.gz
```

**Estado:** Pendiente de ejecución.

### Paso 8 — Alinear las lecturas

```bash
bash scripts/04_align.sh
```

El alineamiento con Subread se completó para las tres muestras. Los archivos BAM se almacenan en `02_bam/` y fueron validados con samtools.

### Paso 9 — Cuantificar las lecturas por gen

```bash
bash scripts/05_featurecounts.sh
```

La cuantificación se realizó utilizando la anotación de exones, el identificador `gene_id` y la configuración no específica de cadena (`-s 0`).

Los resultados se almacenan en `03_counts/`.

### Paso 10 — Integrar los informes de calidad

La integración con MultiQC está pendiente y se realizará después de completar FastQC en las lecturas procesadas.

## 7. Resultados preliminares

### 7.1. Filtrado y trimming con fastp

Las tres muestras se procesaron correctamente con fastp.

| Muestra | Lecturas iniciales | Lecturas conservadas | Porcentaje conservado |
|---|---:|---:|---:|
| `SRR19020908` | 3,075,016 | 3,054,271 | 99.33% |
| `SRR19020909` | 3,562,425 | 3,542,030 | 99.43% |
| `SRR19020910` | 3,629,875 | 3,614,015 | 99.56% |

**Observaciones principales:**

- Se conservó más del 99 % de las lecturas en las tres muestras.
- El porcentaje de bases con calidad Q30 después del filtrado fue aproximadamente del 97 %.
- fastp detectó secuencias relacionadas con el adaptador Illumina TruSeq Read 1.
- Los archivos FASTQ procesados se almacenaron en formato comprimido.

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

Las tres muestras presentan patrones similares de calidad.

Los módulos de calidad por base, composición de bases y contenido de bases N requieren una inspección adicional. En datos de RNA-seq, algunas desviaciones en la composición de bases pueden estar relacionadas con características de la preparación de bibliotecas.

Se detectaron secuencias sobrerrepresentadas asociadas con adaptadores TruSeq, aunque el módulo general de contenido de adaptadores obtuvo un resultado PASS.

Estos resultados justifican comparar los informes antes y después del procesamiento con fastp, sin asumir que todos los avisos representan problemas graves.

### 7.3. Alineamiento y cuantificación por gen

El alineamiento se completó para las tres muestras y los archivos BAM superaron la validación con `samtools quickcheck`.

La cuantificación con featureCounts produjo los siguientes resultados:

| Muestra | Lecturas asignadas a genes | Porcentaje asignado |
|---|---:|---:|
| `SRR19020908` | 2,093,495 | 68.5% |
| `SRR19020909` | 2,563,077 | 72.4% |
| `SRR19020910` | 2,572,872 | 71.2% |

Se compararon las configuraciones de strandedness `-s 0`, `-s 1` y `-s 2`. La configuración **no específica de cadena (`-s 0`)** fue seleccionada para la cuantificación final.

La matriz de conteos se encuentra en:

`03_counts/gene_counts.txt`

Estos resultados corresponden a una evaluación técnica del procesamiento y no permiten, por sí solos, establecer diferencias biológicas entre los grupos.

## 8. Reproducibilidad y gestión de datos

- El archivo `samples.txt` controla la selección de muestras.
- Los scripts Bash automatizan las principales etapas del pipeline.
- Las muestras se procesan de forma secuencial cuando es necesario para limitar el consumo de memoria.
- Los archivos originales se conservan durante la validación.
- Los archivos de secuenciación y los resultados de gran tamaño no se suben a GitHub.
- Los scripts y la documentación se mantienen bajo control de versiones con Git.
- Los parámetros de strandedness y las versiones de las herramientas se documentan para facilitar la reproducibilidad.

### Próximos pasos

1. Ejecutar FastQC en las lecturas procesadas.
2. Comparar los informes FastQC antes y después del trimming.
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

**Project status:** In progress — SRA download, FASTQ conversion, trimming, alignment, gene-level counting, strandedness assessment, and initial FastQC completed. Post-trimming FastQC and MultiQC reporting pending.

**Estado del proyecto:** En desarrollo — Descarga de archivos SRA, conversión a FASTQ, trimming, alineamiento, cuantificación por gen, evaluación de strandedness y FastQC inicial completados. Pendientes el análisis FastQC después del trimming y el informe integrado con MultiQC.

**Last updated / Última actualización:** 2026-10-27
