#!/usr/bin/env bash
set -euo pipefail

# Ejecutar desde la carpeta principal del proyecto.
cd "$(dirname "$0")/../.."

# Integrar los resultados finales del análisis basado en alineamiento.
multiqc \
  results/quality_control/fastqc_before \
  results/quality_control/fastqc_after \
  results/quality_control/fastp \
  results/alignment_based/mapping_statistics \
  results/alignment_based/featurecounts/gene_counts_tail1.txt.summary \
  -c code/config/multiqc_config.yaml \
  -o results/alignment_based/multiqc \
  -n multiqc_final_tail1.html \
  --force

echo "=== Alignment-based MultiQC completed ==="
