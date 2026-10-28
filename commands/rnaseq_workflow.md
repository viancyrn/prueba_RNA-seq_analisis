# 01 — Download sequencing data

## Purpose

Download the sequencing runs associated with GEO accession GSE201955 from the Sequence Read Archive (SRA).

The sample accession numbers are listed in `samples.txt`, allowing the same script to process multiple runs.

## Input

- `samples.txt`: one SRA run accession per line, such as `SRR19020908`.
- Internet access and the NCBI SRA Toolkit.

## Output

- `data/raw_fastq/<SRR>.sra`: the downloaded SRA archive for each run.

## Main commands and parameters

### `prefetch "$SRR" -O "$TEMP_DIR"`

`prefetch` downloads an SRA run using its accession number.

- `"$SRR"` identifies the sequencing run.
- `-O "$TEMP_DIR"` specifies the temporary download location.

The script first downloads each run into a sample-specific temporary directory, then moves the `.sra` file to `data/raw_fastq/`.

### `while IFS= read -r SRR || [[ -n "$SRR" ]]`

This loop reads the accessions from `samples.txt`, including a final line that does not end with a newline.

`IFS=` and `read -r` help preserve the input text without unintended whitespace processing or backslash interpretation.

### `[[ -s ... ]]`

The `-s` test checks whether a file exists and is nonempty.

The script uses this condition to avoid downloading an existing nonempty SRA file and to check that the expected download was produced.

## Reproducibility considerations

- SRA run accessions identify the source sequencing data.
- The download script skips existing nonempty `.sra` files; it does not independently verify their integrity.
- The script expects the directory layout produced by `prefetch` and the accession list in `samples.txt`.
- The script does not currently ignore comment lines in `samples.txt`; this file should contain only accession numbers and optional blank lines.

## Next step

Run `02_fastq.sh` to convert the downloaded SRA archives into FASTQ files.

---

# 02 — Convert SRA archives to FASTQ

## Purpose

Convert the downloaded SRA archives into single-end FASTQ files containing nucleotide sequences and their associated Phred quality scores.

## Input and output

- Input: `data/raw_fastq/<SRR>.sra`
- Output: `data/raw_fastq/<SRR>.fastq`
- Sample accessions: `samples.txt`

## Main command

`fasterq-dump "$INPUT" -O data/raw_fastq -e 4 -p`

### Parameters

- `"$INPUT"`: the SRA archive to convert.
- `-O data/raw_fastq`: write the FASTQ output into the raw-data directory.
- `-e 4`: use four threads for conversion.
- `-p`: display conversion progress.

The script processes one accession at a time, checks that its SRA archive exists and is nonempty, and skips conversion if the expected FASTQ already exists and is nonempty.

## Scientific considerations

FASTQ preserves the sequencing reads and their per-base quality scores, which are required for subsequent quality control and trimming.

These sequencing runs are treated as single-end data. Therefore, the workflow expects one FASTQ file per accession rather than separate read-pair files.

Conversion does not remove adapters, trim low-quality bases, or change the biological interpretation of the reads.

## Limitations

The script checks that the output file is nonempty but does not independently verify that all expected reads were converted. Read counts and quality metrics should be checked in subsequent steps.

## Next step

Run `03_fastp.sh` to generate the final trimmed FASTQ files.

---

# 03 — Read quality control and trimming with fastp

## Purpose

Generate the final quality-filtered FASTQ files for both downstream analysis routes:

1. Genome alignment with Subread and gene counting with featureCounts.
2. Transcriptome pseudoalignment and quantification with Salmon.

Both routes use the same trimmed reads to support a consistent methodological comparison.

## Input and output

- Input: `data/raw_fastq/<SRR>.fastq`
- Trimmed reads: `data/trimmed_fastq/<SRR>_tail1.fastq.gz`
- HTML reports: `results/quality_control/fastp/<SRR>_tail1_fastp.html`
- JSON reports: `results/quality_control/fastp/<SRR>_tail1_fastp.json`

## Main command

`fastp -i "$INPUT" -o "$OUTPUT" -h "$HTML" -j "$JSON" -w 2 --trim_tail1 1`

### Parameters

- `-i`: input single-end FASTQ file.
- `-o`: output compressed FASTQ file.
- `-h`: HTML quality-control report.
- `-j`: JSON report containing machine-readable metrics.
- `-w 2`: use two worker threads.
- `--trim_tail1 1`: remove one nucleotide from the 3′ end of each read.

Other fastp filtering behavior follows the defaults of the installed version.

## Why was the final nucleotide removed?

The original reads were approximately 51 nucleotides long.

Initial quality-control results showed an unfavorable pattern at the final nucleotide, including a reduction in base quality and an increase in ambiguous nucleotide calls.

An initial fastp run using default settings did not adequately resolve this terminal-base pattern.

After reviewing these results, the workflow was revised to remove exactly one nucleotide from the 3′ end of each read.

The final reads are therefore approximately 50 nucleotides long.

This is a dataset-specific decision based on the observed quality-control results, not a general recommendation to trim one nucleotide from every RNA-seq dataset.

## Results and interpretation

The final fastp processing retained more than 99% of reads across the three samples, with approximately 99% of bases at Q30 or higher.

The final read counts were:

| Sample | Reads after fastp |
|---|---:|
| SRR19020908 | 3,055,200 |
| SRR19020909 | 3,543,079 |
| SRR19020910 | 3,615,115 |

These metrics indicate high retained sequencing quality, but they do not independently establish the absence of technical or biological biases.

## Remaining quality-control considerations

The post-trimming FastQC reports still show a per-base sequence-content warning or failure near the beginning of the reads, particularly within approximately the first 13 nucleotides.

This may reflect library-preparation or sequence-composition bias. Its cause has not been independently established.

Removing the final nucleotide does not address this initial sequence-composition pattern.

## Reproducibility and limitations

The script skips a sample when its trimmed FASTQ and both fastp reports already exist and are nonempty.

This avoids unnecessary reprocessing, but file existence alone does not guarantee that the files were generated using the current parameters.

The final fastp JSON reports provide the metrics used for MultiQC reporting.

The original trimming reports are retained temporarily in `results/quality_control/archive_fastp/` during project cleanup. They are not part of the final reporting workflow.

## Next steps

The final trimmed reads are used for genome alignment with Subread and for transcriptome pseudoalignment with Salmon.

FastQC results before and after trimming are included in both final MultiQC reports.

---

# 04 — Genome alignment with Subread

## Purpose

Align the final trimmed single-end RNA-seq reads to the human reference genome using Subread.

This is the alignment-based branch of the workflow. The resulting BAM files are subsequently used for gene-level counting with featureCounts.

## Reference and input

- Reference genome: Homo sapiens GRCh38, Ensembl release 115.
- Subread genome index: `data/references/indexes/subread/GRCh38_subread`
- Input reads: `data/trimmed_fastq/<SRR>_tail1.fastq.gz`
- Sample accessions: `samples.txt`

## Output

- BAM alignment: `data/aligned_bam/<SRR>_tail1.bam`
- Mapping statistics: `results/alignment_based/mapping_statistics/<SRR>_tail1_flagstat.txt`

## Main alignment command

`subread-align -i "$INDEX" -r "$INPUT" -o "$OUTPUT" -t 0 -T "$THREADS"`

### Parameters

- `-i`: reference genome index generated for Subread.
- `-r`: input single-end reads.
- `-o`: output alignment file.
- `-t 0`: specifies RNA-seq reads in the Subread aligner.
- `-T "$THREADS"`: number of threads, set to two in this workflow.

## BAM validation and mapping statistics

`samtools quickcheck "$OUTPUT"`

Performs a basic structural check of the BAM file. A successful check does not guarantee that every alignment record is correct.

`samtools flagstat "$OUTPUT" > "$STATS"`

Summarizes alignment flags, including mapped reads, and saves the statistics in a text file for subsequent MultiQC reporting.

## Observed mapping results

| Sample | Mapped reads | Mapping rate |
|---|---:|---:|
| SRR19020908 | 2,714,036 | 88.83% |
| SRR19020909 | 3,156,878 | 89.10% |
| SRR19020910 | 3,228,911 | 89.32% |

These values refer to the final reads after trimming one nucleotide from the 3′ end.

## Scientific interpretation

Genome mapping measures whether sequencing reads can be aligned to the reference genome under the selected aligner settings.

A mapped read is not necessarily assigned to a gene by featureCounts. Reads may overlap ambiguous regions, fall outside annotated features, or fail counting-specific criteria.

Consequently, the mapping rate and the gene-assignment rate should be reported separately.

## Reproducibility and limitations

The script skips alignment when a nonempty BAM file already exists. It then runs `samtools quickcheck`.

A missing or empty flagstat file is generated automatically, but an existing nonempty flagstat file is not recalculated.

The use of a reference genome and annotation from the same Ensembl release helps maintain consistency across alignment and counting.

## Next step

Use featureCounts to assign aligned reads to annotated genes.

---

# 05 — Gene-level quantification with featureCounts

## Purpose

Assign genome-aligned RNA-seq reads to annotated genes using featureCounts, part of the Subread package.

This step produces a gene-by-sample count matrix for the alignment-based workflow.

## Input and output

- Input: `data/aligned_bam/SRR*_tail1.bam`
- Annotation: `data/references/annotation/Homo_sapiens.GRCh38.115.gtf.gz`
- Gene counts: `results/alignment_based/featurecounts/gene_counts_tail1.txt`
- Assignment summary: `results/alignment_based/featurecounts/gene_counts_tail1.txt.summary`

## Main command

`featureCounts -T 2 -a "$ANNOTATION" -o "$OUTPUT" -t exon -g gene_id -s 0 data/aligned_bam/SRR*_tail1.bam`

### Parameters

- `-T 2`: use two threads.
- `-a`: provide the Ensembl GRCh38 release 115 GTF annotation.
- `-o`: specify the output count matrix.
- `-t exon`: use annotated exons as the features for read assignment.
- `-g gene_id`: summarize exon-level assignments at gene level.
- `-s 0`: count without imposing a strand-specific orientation.

Other counting behavior follows the defaults of the installed featureCounts version.

## Strandedness comparison

Three strandedness settings were explored:

| Setting | Interpretation |
|---|---|
| `-s 0` | Unstranded counting |
| `-s 1` | Stranded counting |
| `-s 2` | Reverse-stranded counting |

Unstranded counting produced the highest number of assigned reads and was selected for the final count matrix.

However, a higher assignment rate does not independently demonstrate that the original RNA-seq library was unstranded. The library preparation protocol or an independent strandedness assessment would be needed to establish this confidently.

Therefore, `-s 0` is a provisional analytical choice rather than a verified library property.

## Final gene-assignment results

| Sample | Assigned reads | Assignment rate |
|---|---:|---:|
| SRR19020908 | 2,093,666 | 68.5% |
| SRR19020909 | 2,563,479 | 72.4% |
| SRR19020910 | 2,572,873 | 71.2% |

The assignment rate is distinct from the genome-mapping rate reported by SAMtools.

## Comparison with the initial processing

After trimming the final nucleotide, the number of assigned reads increased slightly relative to the initial processing:

| Sample | Additional assigned reads |
|---|---:|
| SRR19020908 | 171 |
| SRR19020909 | 402 |
| SRR19020910 | 1 |

These small differences should not be interpreted as strong evidence of improved gene quantification on their own.

## Limitations

- Reads that map to the genome are not necessarily assigned to annotated genes.
- Gene assignment depends on the reference annotation and counting parameters.
- The script uses the BAM files matching `*_tail1.bam`; any additional BAM files matching this pattern would also be included.
- The script does not independently validate the BAM files before counting.
- Strandedness remains unconfirmed.

## Downstream use

The resulting gene count matrix can be imported into R for exploratory analysis and, if the experimental design supports it, differential expression analysis.

The assignment summary is included in the final alignment-based MultiQC report.

---

# 06 — Transcript quantification with Salmon

## Purpose

Estimate transcript-level abundances from the final trimmed single-end RNA-seq reads using Salmon 2.8.0.

This is the pseudoalignment-based branch of the workflow. It uses the same trimmed FASTQ files as the alignment-based branch, but quantifies transcripts against a transcriptome index instead of first aligning reads to the reference genome.

## Reference and index preparation

- Species: Homo sapiens.
- Reference assembly: GRCh38.
- Source: Ensembl release 115.
- Transcriptome: `data/references/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz`
- Salmon index: `data/references/indexes/salmon/salmon_GRCh38_115_index`

The index was created with:

`salmon index -t <transcriptome.fa.gz> -i <index_directory> -k 21`

The index contains approximately 315,794 transcript references after indexing.

The k-mer length was set to 21 because the trimmed reads are approximately 50 nucleotides long. A smaller k-mer size can be appropriate for short reads, although this choice was not independently optimized or benchmarked for this dataset.

The index is transcriptome-only and does not include decoy sequences. This is a methodological limitation when interpreting mapping specificity.

The Ensembl cDNA reference also does not represent every class of noncoding RNA.

## Input and output

- Input: `data/trimmed_fastq/<SRR>_tail1.fastq.gz`
- Quantification directories: `results/pseudoalignment/salmon_quant/<SRR>/`
- Transcript quantification: `<SRR>/quant.sf`
- Salmon metadata: `<SRR>/aux_info/meta_info.json`

## Main command

`salmon quant -i "$INDEX" -r "$FASTQ" -l U --sketch --fldMean 250 --fldSD 25 -p 8 -o "$OUTPUT"`

In the actual script, the input and output paths are constructed separately for each sample.

### Parameters

- `-i`: Salmon transcriptome index.
- `-r`: input single-end FASTQ file.
- `-l U`: assume an unstranded library.
- `--sketch`: use Salmon's sketch-based pseudoalignment workflow.
- `--fldMean 250`: assume a mean fragment length of 250 nucleotides.
- `--fldSD 25`: assume a fragment length standard deviation of 25 nucleotides.
- `-p 8`: use eight threads.
- `-o`: sample-specific output directory.

## Fragment length assumptions

For single-end sequencing, the original fragment length cannot generally be determined directly from the sequenced read alone.

The values 250 and 25 are therefore model assumptions, not measured fragment lengths.

The fragment length distributions shown in the Salmon MultiQC report reflect the specified prior and must not be interpreted as experimentally measured distributions.

The approximately 50-nucleotide read length is not equivalent to a 50-nucleotide RNA fragment length.

## Library strandedness

The analysis uses `-l U` to remain consistent with the provisional unstranded counting choice made for featureCounts.

However, neither this setting nor the higher assignment rate observed with featureCounts `-s 0` independently verifies the library's true strandedness.

## Observed Salmon mapping rates

| Sample | Salmon mapping rate |
|---|---:|
| SRR19020908 | 83.47% |
| SRR19020909 | Approximately 87.4% |
| SRR19020910 | 85.21% |

These values describe mapping to the indexed transcriptome and should not be directly equated with Subread genome-mapping percentages.

## Quantification output

The `quant.sf` file contains transcript-level estimates, including transcript lengths, effective lengths, TPM values, and estimated read counts.

TPM values and estimated counts serve different purposes. TPM represents length-normalized relative abundance, whereas estimated counts can be summarized for downstream gene-level analysis.

## Transcript-to-gene aggregation

Salmon produces transcript-level estimates, while featureCounts produces gene-level counts.

A subsequent R analysis can use `tximport` with a transcript-to-gene mapping to summarize Salmon estimates at gene level.

Gene-level comparison requires consistent gene identifiers and careful handling of transcript versions.

This aggregation has not yet been completed in the current workflow.

## Validation and limitations

The Salmon metadata confirmed successful quantification for all three samples, convergence of the abundance estimation procedure, and pseudoalignment as the mapping method.

Important limitations include:

- The transcriptome-only index does not contain decoy sequences.
- The assumed fragment length distribution was not independently measured.
- Library strandedness remains unconfirmed.
- The cDNA reference does not cover all RNA biotypes.
- Salmon estimated counts are not directly interchangeable with featureCounts integer counts.

## Next step

Use the Salmon quantification outputs for MultiQC reporting and subsequent transcript-to-gene aggregation in R.

---

# 07–08 — Final quality-control reports with MultiQC

## Purpose

Generate two separate, reproducible MultiQC reports for the completed RNA-seq quantification workflows.

Both reports use the same initial and post-trimming quality-control results, but include different quantification-specific metrics.

## Shared quality-control inputs

- `results/quality_control/fastqc_before/`: FastQC reports for the original reads.
- `results/quality_control/fastqc_after/`: FastQC reports for the final trimmed reads.
- `results/quality_control/fastp/`: final fastp reports.
- `code/config/multiqc_config.yaml`: sample-name normalization configuration.

The configuration removes the `_tail1` and `_tail1_flagstat` suffixes where applicable, allowing results to be grouped under consistent sample identifiers.

## 07 — Alignment-based MultiQC

Script: `code/bash/07_multiqc_alignment.sh`

Additional inputs:

- `results/alignment_based/mapping_statistics/`: SAMtools flagstat results.
- `results/alignment_based/featurecounts/gene_counts_tail1.txt.summary`: final featureCounts assignment summary.

Output:

`results/alignment_based/multiqc/multiqc_final_tail1.html`

This report integrates FastQC, fastp, SAMtools, and featureCounts results.

## 08 — Salmon MultiQC

Script: `code/bash/08_multiqc_salmon.sh`

Additional inputs:

- `results/pseudoalignment/salmon_quant/<SRR>/`: Salmon quantification metadata and auxiliary files.

Output:

`results/pseudoalignment/multiqc/multiqc_salmon_final.html`

This report integrates FastQC, fastp, and Salmon results.

The fragment length distributions reported for Salmon reflect the assumed fragment length prior, not directly measured fragment lengths.

## Important MultiQC parameters

- `-c code/config/multiqc_config.yaml`: apply the shared sample-name normalization configuration.
- `-o`: select the output directory.
- `-n`: specify the final HTML report filename.
- `--force`: overwrite an existing report and its associated MultiQC data directory.

Using fixed output names prevents the accumulation of multiple report versions.

## Validation

Both scripts were successfully executed using MultiQC 1.35 after the project directory reorganization.

For the alignment-based report, MultiQC detected three reports each from featureCounts, SAMtools, fastp, and FastQC.

For the Salmon report, MultiQC detected three Salmon metadata reports, three fastp reports, and three FastQC reports.

The HTML reports should also be inspected to confirm correct sample labeling and interpretation of their plots.

## Reporting strategy

The project retains exactly two final integrated MultiQC HTML reports: one for each quantification route.

Intermediate MultiQC reports generated during method development are not part of the final deliverables.


### Detailed strandedness comparison (original BAM files)

The strandedness comparison was performed using the original BAM files, before the final one-nucleotide 3′ trimming.

| Sample | Unstranded (`-s 0`) | Stranded (`-s 1`) | Reverse-stranded (`-s 2`) |
|---|---:|---:|---:|
| SRR19020908 | 2,093,495 | 1,226,534 | 1,225,694 |
| SRR19020909 | 2,563,077 | 1,457,820 | 1,465,912 |
| SRR19020910 | 2,572,872 | 1,449,836 | 1,461,044 |

Unstranded counting yielded the highest number of assigned reads in all three samples. The two strand-specific settings produced similar assignment counts, without a consistent strong preference for either orientation.

This supports using `-s 0` as a provisional counting strategy, but it does not establish the true library strandedness. Independent verification would be necessary.

These exploratory results must not be confused with the final featureCounts results, which were generated from the `_tail1.bam` files using `-s 0`.

---

# Supplementary QC assessment — Trimming optimization

## Comparison of trimming strategies

Three processing stages were evaluated:

1. Original 51-nucleotide reads.
2. Reads processed with the initial default fastp settings.
3. Reads processed with fastp using `--trim_tail1 1`.

The original FastQC reports showed reduced base quality and increased ambiguous-base content at position 51 across all three samples.

The initial default fastp processing did not adequately resolve this terminal-base pattern.

The final strategy therefore removed one nucleotide from the 3′ end of each read.

## Exact fastp results

| Sample | Reads before filtering | Reads after filtering | Reads retained | Q30 bases after filtering |
|---|---:|---:|---:|---:|
| SRR19020908 | 3,075,016 | 3,055,200 | 99.36% | 98.9929% |
| SRR19020909 | 3,562,425 | 3,543,079 | 99.46% | 98.9646% |
| SRR19020910 | 3,629,875 | 3,615,115 | 99.59% | 98.9824% |

## FastQC results after final trimming

| FastQC module | SRR19020908 | SRR19020909 | SRR19020910 |
|---|---|---|---|
| Per base sequence quality | PASS | PASS | PASS |
| Per sequence quality scores | PASS | PASS | PASS |
| Per base N content | PASS | PASS | PASS |
| Per base sequence content | FAIL | FAIL | FAIL |
| Per sequence GC content | PASS | PASS | PASS |
| Adapter Content | PASS | PASS | PASS |
| Overrepresented sequences | PASS | PASS | PASS |
| Sequence Length Distribution | WARNING | WARNING | WARNING |

## Interpretation of the first 13 nucleotides

The initial approximately 13 positions showed uneven nucleotide composition.

Such patterns can arise from RNA-seq library preparation or priming and do not necessarily indicate low base quality.

Removing 13 nucleotides from approximately 51-nucleotide reads would discard a substantial proportion of the available sequence.

Therefore, the first 13 nucleotides were retained. The persistent FastQC sequence-content failure remains a limitation requiring cautious interpretation.

## Subsequent evaluation

Unlike the initial trimming assessment, the completed workflow also evaluated the final reads using Subread and featureCounts.

Compared with the original processing, the final trimming increased assigned reads by 171, 402, and 1 for SRR19020908, SRR19020909, and SRR19020910, respectively.

These differences are small and should not be interpreted as proof that the trimming strategy improved biological quantification.

The final FastQC reports are stored in `results/quality_control/fastqc_after/`.

