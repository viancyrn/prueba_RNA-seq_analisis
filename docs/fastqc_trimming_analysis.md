# FastQC Analysis and Optimization of fastp Trimming

**Dataset:** GSE201955  
**Samples:** SRR19020908, SRR19020909, SRR19020910  
**Tools:** FastQC and fastp  
**Sequencing:** Single-end RNA-seq; original read length: 51 nucleotides

## 1. Objective

I evaluated the sequencing-read quality across three stages to determine whether the initial `fastp` strategy needed adjustment:

1. FastQC on the original reads, before trimming.
2. FastQC after running `fastp` with its default settings.
3. FastQC after modifying `fastp` to remove the final nucleotide.

My decisions were based on the observed quality profiles, not simply on obtaining a PASS for every FastQC module.

## 2. Initial FastQC observations

### 2.1 Per base sequence quality

In all three samples, I noticed reduced quality at the last position of the original 51-nucleotide reads (position 51). Because the same pattern occurred across samples, I considered removing the final nucleotide.

### 2.2 Per base N content

The **Per base N content** plots showed increased ambiguous-base (`N`) content at position 51. Together with the lower base quality, this suggested that the last nucleotide was less reliable than the preceding positions.

### 2.3 Per base sequence content

The proportions of A, T, G and C were uneven across approximately the first 13 positions. Initially, I considered trimming these positions as well. However, nucleotide-composition bias at the start of RNA-seq reads can arise from library preparation or priming and does not, by itself, indicate poor base quality.

Removing 13 nucleotides from reads that are only 51 nucleotides long would discard a substantial amount of sequence information and could affect alignment. I therefore **decided to retain the first 13 bases**, unless additional evidence indicates that they need to be removed.

## 3. First trimming strategy: default fastp

The initial `scripts/03_fastp.sh` used `fastp` without an explicit option to remove the last nucleotide. Its outputs were named `01_clean/${SRR}_clean.fastq.gz`.

I ran FastQC on these processed reads and used the results to reassess the strategy. The observations motivated a targeted test: removing only the last nucleotide rather than trimming a large portion of each read.

## 4. Modification of the fastp script

I added the following option:

```bash
--trim_tail1 1
```

This removes one nucleotide from the 3′ end of each read. The relevant command became:

```bash
fastp \
    -i "$INPUT" \
    -o "$OUTPUT" \
    -h "$HTML" \
    -j "$JSON" \
    -w 2 \
    --trim_tail1 1
```

To avoid overwriting the initial trimming results, I changed the output names:

```bash
OUTPUT="01_clean/${SRR}_tail1.fastq.gz"
HTML="01_clean/${SRR}_tail1_fastp.html"
JSON="01_clean/${SRR}_tail1_fastp.json"
```

This preserved both sets of results for comparison.

## 5. Results after targeted trimming

The modified script completed successfully for all three samples.

| Sample | Reads before filtering | Reads after filtering | Reads retained | Q30 bases after filtering |
|---|---:|---:|---:|---:|
| SRR19020908 | 3,075,016 | 3,055,200 | 99.36% | 98.9929% |
| SRR19020909 | 3,562,425 | 3,543,079 | 99.46% | 98.9646% |
| SRR19020910 | 3,629,875 | 3,615,115 | 99.59% | 98.9824% |

More than 99% of reads were retained in each sample, with approximately 99% of bases at Q30 or above among the filtered reads.

## 6. FastQC assessment after targeted trimming

I ran FastQC on the `_tail1.fastq.gz` files. The resulting module statuses were:

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

### 6.1 Per base sequence quality and N content

Both modules passed in all three samples after targeted trimming. This is consistent with the original observation that the final nucleotide was problematic and supports removing it.

### 6.2 Per base sequence content

This module still failed in all three samples. I did **not** interpret the failure alone as sufficient justification to remove the first 13 nucleotides, given the possible RNA-seq library-preparation bias and the short read length.

### 6.3 Sequence Length Distribution

This module produced a warning. A distribution of read lengths can arise after trimming and filtering; the warning should be interpreted alongside the other quality metrics rather than in isolation.

## 7. Conclusion and next step

The original FastQC reports showed lower base quality and higher N content at position 51 across all three samples. I therefore changed `fastp` to remove one base from the 3′ end using `--trim_tail1 1`.

The revised processing retained more than 99% of reads in every sample. The subsequent FastQC reports passed the **Per base sequence quality** and **Per base N content** checks. I retained the first 13 bases because their uneven nucleotide composition alone did not establish a quality problem.

**Chosen files for downstream analysis:** `01_clean/${SRR}_tail1.fastq.gz`.

The effect of this decision on alignment and gene quantification has **not yet been evaluated**. The downstream alignment script must be checked to ensure that it uses the `_tail1.fastq.gz` files rather than the earlier `_clean.fastq.gz` outputs.

## 8. Supporting FastQC reports

FastQC HTML reports for all three samples were generated at each stage:

- `04_multiqc/fastqc_before/` — original reads.
- `04_multiqc/fastqc_after/` — reads after default `fastp`.
- `04_multiqc/fastqc_tail1/` — reads after `fastp --trim_tail1 1`.

The figures used to guide the trimming decision were **Per base sequence quality**, **Per base N content**, and **Per base sequence content**. The HTML reports include these plots and the other FastQC modules.
