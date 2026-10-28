#!/usr/bin/env bash
set -euo pipefail

# Ejecutar desde la carpeta principal del proyecto.
cd "$(dirname "$0")/../.."

# Definir las rutas y los parámetros de Salmon.
INDEX="data/references/indexes/salmon/salmon_GRCh38_115_index"
INPUT_DIR="data/trimmed_fastq"
OUTPUT_DIR="results/pseudoalignment/salmon_quant"
THREADS=8

mkdir -p "$OUTPUT_DIR"

# Cuantificar los tres samples con los mismos parámetros.
for sample in SRR19020908 SRR19020909 SRR19020910; do
    echo "=== Salmon quantification: ${sample} ==="

    salmon quant \
        -i "$INDEX" \
        -r "${INPUT_DIR}/${sample}_tail1.fastq.gz" \
        -l U \
        --sketch \
        --fldMean 250 \
        --fldSD 25 \
        -p "$THREADS" \
        -o "${OUTPUT_DIR}/${sample}"

    echo "=== Completed: ${sample} ==="
done

echo "=== All Salmon quantifications completed ==="
