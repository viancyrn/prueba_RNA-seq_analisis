#!/usr/bin/env bash
set -euo pipefail

# Localizar la carpeta principal del proyecto.
cd "$(dirname "$0")/../.."

# Definir las rutas y los parámetros del alineamiento.
INDEX="data/references/indexes/subread/GRCh38_subread"
INPUT_DIR="data/trimmed_fastq"
OUTPUT_DIR="data/aligned_bam"
STATS_DIR="results/alignment_based/mapping_statistics"
THREADS=2

mkdir -p "$OUTPUT_DIR" "$STATS_DIR"

# Procesar cada muestra incluida en samples.txt.
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    [[ -z "$SAMPLE" ]] && continue

    INPUT="${INPUT_DIR}/${SAMPLE}_tail1.fastq.gz"
    OUTPUT="${OUTPUT_DIR}/${SAMPLE}_tail1.bam"
    STATS="${STATS_DIR}/${SAMPLE}_tail1_flagstat.txt"

    # Alinear solamente si el archivo BAM no existe.
    if [[ ! -s "$OUTPUT" ]]; then
        if [[ ! -f "$INPUT" ]]; then
            echo "ERROR: FASTQ not found: $INPUT" >&2
            exit 1
        fi

        echo "=== Aligning ${SAMPLE} ==="
        subread-align \
            -i "$INDEX" \
            -r "$INPUT" \
            -o "$OUTPUT" \
            -t 0 \
            -T "$THREADS"
    else
        echo "=== Existing BAM: ${SAMPLE} ==="
    fi

    # Comprobar la estructura básica del archivo BAM.
    samtools quickcheck "$OUTPUT"

    # Generar estadísticas solamente cuando faltan.
    if [[ ! -s "$STATS" ]]; then
        samtools flagstat "$OUTPUT" > "$STATS"
    fi

    echo "=== BAM validation completed: ${SAMPLE} ==="
done < samples.txt

echo "=== All samples processed ==="
