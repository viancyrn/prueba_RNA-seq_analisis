
#!/bin/bash
set -euo pipefail

# Alineamiento de lecturas RNA-seq contra el genoma humano GRCh38.
# Ejecutar desde la carpeta principal del proyecto.

INDEX="ref/GRCh38_subread"
INPUT_DIR="01_clean"
OUTPUT_DIR="02_bam"
THREADS=2

# Crear el directorio de salida si no existe.
mkdir -p "$OUTPUT_DIR"

# Procesar las muestras indicadas en samples.txt.
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do

    # Ignorar líneas vacías.
    [[ -z "$SAMPLE" ]] && continue

    INPUT="${INPUT_DIR}/${SAMPLE}_clean.fastq.gz"
    OUTPUT="${OUTPUT_DIR}/${SAMPLE}.bam"

    # Evitar repetir alineamientos ya realizados.
    if [[ -s "$OUTPUT" ]]; then
        echo "Omitiendo ${SAMPLE}: archivo BAM existente."
        continue
    fi

    # Verificar que exista el archivo FASTQ.
    if [[ ! -f "$INPUT" ]]; then
        echo "ERROR: no se encontró $INPUT" >&2
        exit 1
    fi

    echo "Alineando ${SAMPLE}..."

    subread-align \
        -i "$INDEX" \
        -r "$INPUT" \
        -o "$OUTPUT" \
        -t 0 \
        -T "$THREADS"

    echo "Alineamiento finalizado: ${SAMPLE}"

done < samples.txt

echo "Procesamiento de muestras finalizado."

