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
| `SRR19020910` | Additional test sample |

The third accession was added to verify that the scripts can automatically process multiple samples.

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
   3. Quality control and trimming
          fastp
           |
           v
   4. Reference genome preparation
         GRCh38
           |
           v
   5. Read alignment
         Subread
           |
           v
   6. Gene-level quantification
       featureCounts
           |
           v
   7. Summary quality report
         MultiQC
```

### Pipeline progress

| Step | Tool | Status |
|---|---|---|
| 1. SRA download | SRA Toolkit (`prefetch`) | ✅ Completed |
| 2. SRA-to-FASTQ conversion | SRA Toolkit (`fasterq-dump`) | ✅ Completed |
| 3. Quality control and trimming | `fastp` | ✅ Completed |
| 4. Reference genome preparation | GRCh38 / Subread | ⏳ Planned |
| 5. Read alignment | Subread | ⏳ Planned |
| 6. Gene-level read counting | featureCounts | ⏳ Planned |
| 7. Integrated quality report | MultiQC | ⏳ Planned |

## 4. Repository structure

The GitHub repository contains the scripts and metadata needed to document the workflow.

```text
prueba_RNA-seq_analisis/
├── README.md
├── .gitignore
├── samples.txt
└── scripts/
    ├── 01_download.sh
    ├── 02_fastq.sh
    └── 03_fastp.sh
```

During execution, the pipeline also uses the following local directories:

| Directory | Purpose |
|---|---|
| `00_raw/` | Downloaded SRA files and raw FASTQ files |
| `01_clean/` | Trimmed FASTQ files and fastp reports |
| `02_bam/` | Alignment files (planned) |
| `03_counts/` | Gene count tables (planned) |
| `04_multiqc/` | Integrated quality reports (planned) |
| `ref/` | Reference genome and index files (planned) |

**Large sequencing files and analysis outputs are excluded from version control using `.gitignore`.**

## 5. Software and environment

The pipeline is developed on **macOS (Apple Silicon)** using a Conda environment named `rnaseq_qc`.

| Software | Version | Purpose |
|---|---|---|
| SRA Toolkit | 3.4.1 | Download and FASTQ conversion |
| fastp | 1.4.0 | Quality control and trimming |
| Subread | 2.1.1 | Read alignment |
| featureCounts | Included in Subread | Gene-level read counting |
| MultiQC | 1.35 | Quality-report integration |
| Bash | — | Workflow automation |

Subread, featureCounts, and MultiQC are available in the environment, but their pipeline stages have not yet been completed.

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

El objetivo inicial es desarrollar, automatizar y validar las etapas de preprocesamiento mediante tres accesiones SRA. Posteriormente, el proyecto incluirá el alineamiento, la cuantificación de genes y la integración de los informes de calidad.

> **Nota:** Los resultados actuales corresponden a pruebas del pipeline y evaluación de calidad. No constituyen un análisis de expresión diferencial.

## 2. Datos y muestras

**Fuente de datos:** NCBI Gene Expression Omnibus (GEO)

| Identificador | Accesión |
|---|---|
| GEO Series | [GSE201955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201955) |
| BioProject | [PRJNA833856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833856) |
| SRA Study | SRP373003 |

### Muestras utilizadas

| Accesión SRA | Descripción |
|---|---|
| `SRR19020908` | Control sano |
| `SRR19020909` | Asma |
| `SRR19020910` | Muestra adicional de prueba |

La tercera muestra se agregó para comprobar el procesamiento automático de múltiples accesiones.

## 3. Etapas del pipeline

| Etapa | Herramienta | Estado |
|---|---|---|
| 1. Descarga de archivos SRA | `prefetch` | ✅ Completada |
| 2. Conversión de SRA a FASTQ | `fasterq-dump` | ✅ Completada |
| 3. Control de calidad y trimming | `fastp` | ✅ Completada |
| 4. Preparación del genoma de referencia | GRCh38 / Subread | ⏳ Pendiente |
| 5. Alineamiento de lecturas | Subread | ⏳ Pendiente |
| 6. Cuantificación de lecturas por gen | featureCounts | ⏳ Pendiente |
| 7. Informe integrado de calidad | MultiQC | ⏳ Pendiente |

## 4. Organización del proyecto

Los scripts y la lista de muestras se almacenan en GitHub. Los archivos grandes de secuenciación y los resultados del análisis permanecen fuera del repositorio.

### Directorios de trabajo

| Directorio | Contenido |
|---|---|
| `00_raw/` | Archivos SRA y FASTQ originales |
| `01_clean/` | FASTQ procesados e informes de fastp |
| `02_bam/` | Archivos de alineamiento (pendiente) |
| `03_counts/` | Tablas de conteos por gen (pendiente) |
| `04_multiqc/` | Informes de calidad integrados (pendiente) |
| `ref/` | Genoma de referencia e índices (pendiente) |

## 5. Ejecución del preprocesamiento

### Activar el entorno

```bash
conda activate rnaseq_qc
```

### Descargar los archivos SRA

```bash
bash scripts/01_download.sh
```

### Convertir los archivos SRA a FASTQ

```bash
bash scripts/02_fastq.sh
```

### Realizar el control de calidad y trimming

```bash
bash scripts/03_fastp.sh
```

Los scripts leen las accesiones de `samples.txt`, procesan las muestras de manera secuencial y omiten los resultados completos que ya existen.

## 6. Resultados preliminares

| Muestra | Lecturas iniciales | Lecturas conservadas | Porcentaje conservado |
|---|---:|---:|---:|
| `SRR19020908` | 3,075,016 | 3,054,271 | 99.33% |
| `SRR19020909` | 3,562,425 | 3,542,030 | 99.43% |
| `SRR19020910` | 3,629,875 | 3,614,015 | 99.56% |

### Observaciones principales

- Se conservó más del **99% de las lecturas** en las tres muestras.
- El porcentaje de bases con calidad **Q30** después del filtrado fue aproximadamente del 97%.
- Se identificaron secuencias del adaptador **Illumina TruSeq Read 1**.
- Los archivos FASTQ procesados se almacenaron en formato comprimido.

Los resultados indican que el preprocesamiento se ejecutó correctamente. No permiten todavía evaluar diferencias de expresión génica entre grupos.

## 7. Reproducibilidad

El pipeline utiliza una lista de muestras y scripts Bash para automatizar las etapas del análisis.

Los archivos de secuenciación no se almacenan en GitHub debido a su tamaño. Los datos pueden volver a descargarse desde las accesiones públicas de SRA.

La documentación y los scripts se actualizan a medida que se completan nuevas etapas.

---

## References / Referencias

- [NCBI GEO — GSE201955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201955)
- [NCBI BioProject — PRJNA833856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA833856)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [fastp](https://github.com/OpenGene/fastp)
- [Subread](https://subread.sourceforge.net/)
- [MultiQC](https://multiqc.info/)

---

**Project status:** Work in progress.  
**Estado del proyecto:** En desarrollo.
