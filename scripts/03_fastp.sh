
#!/usr/bin/env bash
set -euo pipefail

# FASE 3: CONTROL DE CALIDAD Y TRIMMING CON FASTP
# Dataset: GSE201955

# Ir al directorio principal del proyecto
cd "$(dirname "$0")/.."

# Crear el directorio de resultados
mkdir -p 01_clean

while IFS= read -r SRR || [[ -n "$SRR" ]]; do

    # Ignorar líneas vacías y comentarios
    [[ -z "$SRR" || "$SRR" == \#* ]] && continue

    # Definir archivos de entrada y salida
    INPUT="00_raw/${SRR}.fastq"
    OUTPUT="01_clean/${SRR}_clean.fastq.gz"
    HTML="01_clean/${SRR}_fastp.html"
    JSON="01_clean/${SRR}_fastp.json"

    # Comprobar que existe el archivo FASTQ
    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: No existe $INPUT"
        exit 1
    fi

    # Saltar muestras ya procesadas
    if [[ -s "$OUTPUT" && -s "$HTML" && -s "$JSON" ]]; then
        echo "$SRR ya fue procesado con fastp. Se omite."
        continue
    fi

    echo "Procesando $SRR con fastp ..."

    # Ejecutar control de calidad y trimming
    fastp \
        -i "$INPUT" \
        -o "$OUTPUT" \
        -h "$HTML" \
        -j "$JSON" \
        -w 2

    # Verificar que se generaron los resultados
    if [[ ! -s "$OUTPUT" || ! -s "$HTML" || ! -s "$JSON" ]]; then
        echo "ERROR: Resultados incompletos para $SRR"
        exit 1
    fi

    echo "$SRR procesado correctamente."

done < samples.txt

echo "Control de calidad y trimming finalizados."
