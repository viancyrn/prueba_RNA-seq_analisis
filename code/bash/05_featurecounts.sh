#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

# Contar lecturas asignadas a genes usando featureCounts.
# Ejecutar desde el directorio principal del proyecto.

ANNOTATION="data/references/annotation/Homo_sapiens.GRCh38.115.gtf.gz"
OUTPUT="results/alignment_based/featurecounts/gene_counts_tail1.txt"

# Crear el directorio de resultados.
mkdir -p results/alignment_based/featurecounts

# Verificar que exista la anotación.
if [[ ! -f "$ANNOTATION" ]]; then
    echo "ERROR: No se encontró la anotación GTF." >&2
    exit 1
fi

# Contar lecturas por gen en las tres muestras.
# -s 0: conteo no específico de cadena (provisional).
featureCounts \
    -T 2 \
    -a "$ANNOTATION" \
    -o "$OUTPUT" \
    -t exon \
    -g gene_id \
    -s 0 \
    data/aligned_bam/SRR*_tail1.bam

echo "Conteo de genes finalizado."

