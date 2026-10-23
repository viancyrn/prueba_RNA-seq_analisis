
#!/usr/bin/env bash
set -euo pipefail

# FASE 1: DESCARGA DE DATOS SRA
# Dataset: GSE201955

# Ir al directorio principal del proyecto
cd "$(dirname "$0")/.."

mkdir -p 00_raw

while IFS= read -r SRR || [[ -n "$SRR" ]]; do

    # Ignorar líneas vacías
    [[ -z "$SRR" ]] && continue

    # Saltar muestras ya descargadas
    if [[ -s "00_raw/${SRR}.sra" ]]; then
        echo "$SRR ya existe. Se omite."
        continue
    fi

    echo "Descargando $SRR ..."

    # Descargar en un directorio temporal
    TEMP_DIR="00_raw/.download_${SRR}"
    mkdir -p "$TEMP_DIR"

    prefetch "$SRR" -O "$TEMP_DIR"

    # Localizar el archivo SRA descargado
    SOURCE="$TEMP_DIR/$SRR/$SRR.sra"

    if [[ ! -s "$SOURCE" ]]; then
        echo "ERROR: No se encontró el archivo de $SRR"
        exit 1
    fi

    # Guardar directamente en 00_raw
    mv "$SOURCE" "00_raw/${SRR}.sra"

    # Eliminar directorios temporales vacíos
    rmdir "$TEMP_DIR/$SRR" "$TEMP_DIR"

    echo "$SRR descargado correctamente."

done < samples.txt

echo "Descarga finalizada."

