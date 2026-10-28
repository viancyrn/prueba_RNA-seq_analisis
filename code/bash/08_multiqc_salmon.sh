#!/usr/bin/env bash
set -euo pipefail

# Ejecutar desde la carpeta principal del proyecto.
cd "$(dirname "$0")/../.."

# Integrar los resultados finales del análisis con Salmon.
multiqc \
  results/quality_control/fastqc_before \
  results/quality_control/fastqc_after \
  results/quality_control/fastp/*_tail1_fastp.json \
  results/pseudoalignment/salmon_quant/SRR19020908 \
  results/pseudoalignment/salmon_quant/SRR19020909 \
  results/pseudoalignment/salmon_quant/SRR19020910 \
  -c code/config/multiqc_config.yaml \
  -o results/pseudoalignment/multiqc \
  -n multiqc_salmon_final.html \
  --force

echo "=== Salmon MultiQC completed ==="
